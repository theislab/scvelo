"""Main module."""
import logging
import warnings
from typing import (
    Callable,
    Iterable,
    List,
    Literal,
    NamedTuple,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from joblib import delayed, Parallel

import numpy as np
import pandas as pd
import torch
import torch.nn.functional as F
from scipy.stats import ttest_ind
from torch import nn
from torch.distributions import Categorical, Dirichlet
from torch.distributions import kl_divergence as kl
from torch.distributions import MixtureSameFamily, Normal

from anndata import AnnData
from scvi.data import AnnDataManager
from scvi.data.fields import LayerField
from scvi.dataloaders import DataSplitter
from scvi.model.base import BaseModelClass, UnsupervisedTrainingMixin, VAEMixin
from scvi.module.base import auto_move_data, BaseModuleClass, LossOutput
from scvi.nn import Encoder, FCLayers
from scvi.train import TrainingPlan, TrainRunner
from scvi.utils._docstrings import setup_anndata_dsp

torch.backends.cudnn.benchmark = True


class _REGISTRY_KEYS_NT(NamedTuple):
    X_KEY: str = "X"
    U_KEY: str = "U"


REGISTRY_KEYS = _REGISTRY_KEYS_NT()

DEFAULT_ACTIVATION_FUNCTION = torch.nn.Softplus()


class DecoderVELOVI(nn.Module):
    """Decodes data from latent space of ``n_input`` dimensions ``n_output`` dimensions.

    Uses a fully-connected neural network of ``n_hidden`` layers.

    Parameters
    ----------
    n_input
        The dimensionality of the input (latent space)
    n_output
        The dimensionality of the output (data space)
    n_cat_list
        A list containing the number of categories
        for each category of interest. Each category will be
        included using a one-hot encoding
    n_layers
        The number of fully-connected hidden layers
    n_hidden
        The number of nodes per hidden layer
    dropout_rate
        Dropout rate to apply to each of the hidden layers
    inject_covariates
        Whether to inject covariates in each layer, or just the first (default).
    use_batch_norm
        Whether to use batch norm in layers
    use_layer_norm
        Whether to use layer norm in layers
    linear_decoder
        Whether to use linear decoder for time
    """

    def __init__(
        self,
        n_input: int,
        n_output: int,
        n_cat_list: Iterable[int] = None,
        n_layers: int = 1,
        n_hidden: int = 128,
        inject_covariates: bool = True,
        use_batch_norm: bool = True,
        use_layer_norm: bool = False,
        dropout_rate: float = 0.0,
        linear_decoder: bool = False,
        **kwargs,
    ):
        super().__init__()
        self.n_ouput = n_output
        self.linear_decoder = linear_decoder
        self.rho_first_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden if not linear_decoder else n_output,
            n_cat_list=n_cat_list,
            n_layers=n_layers if not linear_decoder else 1,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm if not linear_decoder else False,
            use_activation=not linear_decoder,
            bias=not linear_decoder,
            **kwargs,
        )

        self.pi_first_decoder = FCLayers(
            n_in=n_input,
            n_out=n_hidden,
            n_cat_list=n_cat_list,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            inject_covariates=inject_covariates,
            use_batch_norm=use_batch_norm,
            use_layer_norm=use_layer_norm,
            **kwargs,
        )

        # categorical pi
        # 4 states
        self.px_pi_decoder = nn.Linear(n_hidden, 4 * n_output)

        # rho for induction
        self.px_rho_decoder = nn.Sequential(nn.Linear(n_hidden, n_output), nn.Sigmoid())

        # tau for repression
        self.px_tau_decoder = nn.Sequential(nn.Linear(n_hidden, n_output), nn.Sigmoid())

        self.linear_scaling_tau = nn.Parameter(torch.zeros(n_output))
        self.linear_scaling_tau_intercept = nn.Parameter(torch.zeros(n_output))

    def forward(self, z: torch.Tensor, latent_dim: int = None):
        """The forward computation for a single sample.

         #. Decodes the data from the latent space using the decoder network
         #. Returns parameters for the ZINB distribution of expression
         #. If ``dispersion != 'gene-cell'`` then value for that param will be ``None``

        Parameters
        ----------
        z :
            tensor with shape ``(n_input,)``
        cat_list
            list of category membership(s) for this sample

        Returns
        -------
        4-tuple of :py:class:`torch.Tensor`
            parameters for the ZINB distribution of expression

        """
        z_in = z
        if latent_dim is not None:
            mask = torch.zeros_like(z)
            mask[..., latent_dim] = 1
            z_in = z * mask
        # The decoder returns values for the parameters of the ZINB distribution
        rho_first = self.rho_first_decoder(z_in)

        if not self.linear_decoder:
            px_rho = self.px_rho_decoder(rho_first)
            px_tau = self.px_tau_decoder(rho_first)
        else:
            px_rho = nn.Sigmoid()(rho_first)
            px_tau = 1 - nn.Sigmoid()(
                rho_first * self.linear_scaling_tau.exp()
                + self.linear_scaling_tau_intercept
            )

        # cells by genes by 4
        pi_first = self.pi_first_decoder(z)
        px_pi = nn.Softplus()(
            torch.reshape(self.px_pi_decoder(pi_first), (z.shape[0], self.n_ouput, 4))
        )

        return px_pi, px_rho, px_tau


# VAE model
class VELOVAE(BaseModuleClass):
    """Variational auto-encoder model.

    This is an implementation of the veloVI model descibed in :cite:p:`GayosoWeiler2022`

    Parameters
    ----------
    n_input
        Number of input genes
    n_hidden
        Number of nodes per hidden layer
    n_latent
        Dimensionality of the latent space
    n_layers
        Number of hidden layers used for encoder and decoder NNs
    dropout_rate
        Dropout rate for neural networks
    log_variational
        Log(data+1) prior to encoding for numerical stability. Not normalization.
    latent_distribution
        One of

        * ``'normal'`` - Isotropic normal
        * ``'ln'`` - Logistic normal with normal params N(0, 1)
    use_layer_norm
        Whether to use layer norm in layers
    use_observed_lib_size
        Use observed library size for RNA as scaling factor in mean of conditional distribution
    var_activation
        Callable used to ensure positivity of the variational distributions' variance.
        When `None`, defaults to `torch.exp`.
    """

    def __init__(
        self,
        n_input: int,
        true_time_switch: Optional[np.ndarray] = None,
        n_hidden: int = 128,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        log_variational: bool = False,
        latent_distribution: str = "normal",
        use_batch_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_layer_norm: Literal["encoder", "decoder", "none", "both"] = "both",
        use_observed_lib_size: bool = True,
        var_activation: Optional[Callable] = DEFAULT_ACTIVATION_FUNCTION,
        model_steady_states: bool = True,
        gamma_unconstr_init: Optional[np.ndarray] = None,
        alpha_unconstr_init: Optional[np.ndarray] = None,
        alpha_1_unconstr_init: Optional[np.ndarray] = None,
        lambda_alpha_unconstr_init: Optional[np.ndarray] = None,
        switch_spliced: Optional[np.ndarray] = None,
        switch_unspliced: Optional[np.ndarray] = None,
        t_max: float = 20,
        penalty_scale: float = 0.2,
        dirichlet_concentration: float = 0.25,
        linear_decoder: bool = False,
        time_dep_transcription_rate: bool = False,
    ):
        super().__init__()
        self.n_latent = n_latent
        self.log_variational = log_variational
        self.latent_distribution = latent_distribution
        self.use_observed_lib_size = use_observed_lib_size
        self.n_input = n_input
        self.model_steady_states = model_steady_states
        self.t_max = t_max
        self.penalty_scale = penalty_scale
        self.dirichlet_concentration = dirichlet_concentration
        self.time_dep_transcription_rate = time_dep_transcription_rate

        if switch_spliced is not None:
            self.register_buffer("switch_spliced", torch.from_numpy(switch_spliced))
        else:
            self.switch_spliced = None
        if switch_unspliced is not None:
            self.register_buffer("switch_unspliced", torch.from_numpy(switch_unspliced))
        else:
            self.switch_unspliced = None

        n_genes = n_input * 2

        # switching time
        self.switch_time_unconstr = torch.nn.Parameter(7 + 0.5 * torch.randn(n_input))
        if true_time_switch is not None:
            self.register_buffer("true_time_switch", torch.from_numpy(true_time_switch))
        else:
            self.true_time_switch = None

        # degradation
        if gamma_unconstr_init is None:
            self.gamma_mean_unconstr = torch.nn.Parameter(-1 * torch.ones(n_input))
        else:
            self.gamma_mean_unconstr = torch.nn.Parameter(
                torch.from_numpy(gamma_unconstr_init)
            )

        # splicing
        # first samples around 1
        self.beta_mean_unconstr = torch.nn.Parameter(0.5 * torch.ones(n_input))

        # transcription
        if alpha_unconstr_init is None:
            self.alpha_unconstr = torch.nn.Parameter(0 * torch.ones(n_input))
        else:
            self.alpha_unconstr = torch.nn.Parameter(
                torch.from_numpy(alpha_unconstr_init)
            )

        # TODO: Add `require_grad`
        if alpha_1_unconstr_init is None:
            self.alpha_1_unconstr = torch.nn.Parameter(0 * torch.ones(n_input))
        else:
            self.alpha_1_unconstr = torch.nn.Parameter(
                torch.from_numpy(alpha_1_unconstr_init)
            )
        self.alpha_1_unconstr.requires_grad = time_dep_transcription_rate

        if lambda_alpha_unconstr_init is None:
            self.lambda_alpha_unconstr = torch.nn.Parameter(0 * torch.ones(n_input))
        else:
            self.lambda_alpha_unconstr = torch.nn.Parameter(
                torch.from_numpy(lambda_alpha_unconstr_init)
            )
        self.lambda_alpha_unconstr.requires_grad = time_dep_transcription_rate

        # likelihood dispersion
        # for now, with normal dist, this is just the variance
        self.scale_unconstr = torch.nn.Parameter(-1 * torch.ones(n_genes, 4))

        use_batch_norm_encoder = use_batch_norm == "encoder" or use_batch_norm == "both"
        use_batch_norm_decoder = use_batch_norm == "decoder" or use_batch_norm == "both"
        use_layer_norm_encoder = use_layer_norm == "encoder" or use_layer_norm == "both"
        use_layer_norm_decoder = use_layer_norm == "decoder" or use_layer_norm == "both"
        self.use_batch_norm_decoder = use_batch_norm_decoder

        # z encoder goes from the n_input-dimensional data to an n_latent-d
        # latent space representation
        n_input_encoder = n_genes
        self.z_encoder = Encoder(
            n_input_encoder,
            n_latent,
            n_layers=n_layers,
            n_hidden=n_hidden,
            dropout_rate=dropout_rate,
            distribution=latent_distribution,
            use_batch_norm=use_batch_norm_encoder,
            use_layer_norm=use_layer_norm_encoder,
            var_activation=var_activation,
            activation_fn=torch.nn.ReLU,
        )
        # decoder goes from n_latent-dimensional space to n_input-d data
        n_input_decoder = n_latent
        self.decoder = DecoderVELOVI(
            n_input_decoder,
            n_input,
            n_layers=n_layers,
            n_hidden=n_hidden,
            use_batch_norm=use_batch_norm_decoder,
            use_layer_norm=use_layer_norm_decoder,
            activation_fn=torch.nn.ReLU,
            linear_decoder=linear_decoder,
        )

    def _get_inference_input(self, tensors):
        spliced = tensors[REGISTRY_KEYS.X_KEY]
        unspliced = tensors[REGISTRY_KEYS.U_KEY]

        input_dict = {
            "spliced": spliced,
            "unspliced": unspliced,
        }
        return input_dict

    def _get_generative_input(self, tensors, inference_outputs):
        z = inference_outputs["z"]
        gamma = inference_outputs["gamma"]
        beta = inference_outputs["beta"]
        alpha = inference_outputs["alpha"]
        alpha_1 = inference_outputs["alpha_1"]
        lambda_alpha = inference_outputs["lambda_alpha"]

        input_dict = {
            "z": z,
            "gamma": gamma,
            "beta": beta,
            "alpha": alpha,
            "alpha_1": alpha_1,
            "lambda_alpha": lambda_alpha,
        }
        return input_dict

    @auto_move_data
    def inference(
        self,
        spliced,
        unspliced,
        n_samples=1,
    ):
        """High level inference method.

        Runs the inference (encoder) model.
        """
        spliced_ = spliced
        unspliced_ = unspliced
        if self.log_variational:
            spliced_ = torch.log(0.01 + spliced)
            unspliced_ = torch.log(0.01 + unspliced)

        encoder_input = torch.cat((spliced_, unspliced_), dim=-1)

        qz_m, qz_v, z = self.z_encoder(encoder_input)

        if n_samples > 1:
            qz_m = qz_m.unsqueeze(0).expand((n_samples, qz_m.size(0), qz_m.size(1)))
            qz_v = qz_v.unsqueeze(0).expand((n_samples, qz_v.size(0), qz_v.size(1)))
            # when z is normal, untran_z == z
            untran_z = Normal(qz_m, qz_v.sqrt()).sample()
            z = self.z_encoder.z_transformation(untran_z)

        gamma, beta, alpha, alpha_1, lambda_alpha = self._get_rates()

        outputs = {
            "z": z,
            "qz_m": qz_m,
            "qz_v": qz_v,
            "gamma": gamma,
            "beta": beta,
            "alpha": alpha,
            "alpha_1": alpha_1,
            "lambda_alpha": lambda_alpha,
        }
        return outputs

    def _get_rates(self):
        # globals
        # degradation
        gamma = torch.clamp(F.softplus(self.gamma_mean_unconstr), 0, 50)
        # splicing
        beta = torch.clamp(F.softplus(self.beta_mean_unconstr), 0, 50)
        # transcription
        alpha = torch.clamp(F.softplus(self.alpha_unconstr), 0, 50)
        if self.time_dep_transcription_rate:
            alpha_1 = torch.clamp(F.softplus(self.alpha_1_unconstr), 0, 50)
            lambda_alpha = torch.clamp(F.softplus(self.lambda_alpha_unconstr), 0, 50)
        else:
            alpha_1 = self.alpha_1_unconstr
            lambda_alpha = self.lambda_alpha_unconstr

        return gamma, beta, alpha, alpha_1, lambda_alpha

    @auto_move_data
    def generative(self, z, gamma, beta, alpha, alpha_1, lambda_alpha, latent_dim=None):
        """Runs the generative model."""
        decoder_input = z
        px_pi_alpha, px_rho, px_tau = self.decoder(decoder_input, latent_dim=latent_dim)
        px_pi = Dirichlet(px_pi_alpha).rsample()

        scale_unconstr = self.scale_unconstr
        scale = F.softplus(scale_unconstr)

        mixture_dist_s, mixture_dist_u, end_penalty = self.get_px(
            px_pi,
            px_rho,
            px_tau,
            scale,
            gamma,
            beta,
            alpha,
            alpha_1,
            lambda_alpha,
        )

        return {
            "px_pi": px_pi,
            "px_rho": px_rho,
            "px_tau": px_tau,
            "scale": scale,
            "px_pi_alpha": px_pi_alpha,
            "mixture_dist_u": mixture_dist_u,
            "mixture_dist_s": mixture_dist_s,
            "end_penalty": end_penalty,
        }

    # TODO: Add docstrings
    def loss(
        self,
        tensors,
        inference_outputs,
        generative_outputs,
        kl_weight: float = 1.0,
        n_obs: float = 1.0,
    ):
        """TODO."""
        spliced = tensors[REGISTRY_KEYS.X_KEY]
        unspliced = tensors[REGISTRY_KEYS.U_KEY]

        qz_m = inference_outputs["qz_m"]
        qz_v = inference_outputs["qz_v"]

        px_pi = generative_outputs["px_pi"]
        px_pi_alpha = generative_outputs["px_pi_alpha"]

        end_penalty = generative_outputs["end_penalty"]
        mixture_dist_s = generative_outputs["mixture_dist_s"]
        mixture_dist_u = generative_outputs["mixture_dist_u"]

        kl_divergence_z = kl(Normal(qz_m, torch.sqrt(qz_v)), Normal(0, 1)).sum(dim=1)

        reconst_loss_s = -mixture_dist_s.log_prob(spliced)
        reconst_loss_u = -mixture_dist_u.log_prob(unspliced)

        reconst_loss = reconst_loss_u.sum(dim=-1) + reconst_loss_s.sum(dim=-1)

        kl_pi = kl(
            Dirichlet(px_pi_alpha),
            Dirichlet(self.dirichlet_concentration * torch.ones_like(px_pi)),
        ).sum(dim=-1)

        # local loss
        kl_local = kl_divergence_z + kl_pi
        weighted_kl_local = kl_weight * (kl_divergence_z) + kl_pi

        local_loss = torch.mean(reconst_loss + weighted_kl_local)

        loss = local_loss + self.penalty_scale * (1 - kl_weight) * end_penalty

        loss_recoder = LossOutput(
            loss=loss, reconstruction_loss=reconst_loss, kl_local=kl_local
        )

        return loss_recoder

    # TODO: Add docstrings
    @auto_move_data
    def get_px(
        self,
        px_pi,
        px_rho,
        px_tau,
        scale,
        gamma,
        beta,
        alpha,
        alpha_1,
        lambda_alpha,
    ) -> torch.Tensor:
        """TODO."""
        t_s = torch.clamp(F.softplus(self.switch_time_unconstr), 0, self.t_max)

        n_cells = px_pi.shape[0]

        # component dist
        comp_dist = Categorical(probs=px_pi)

        # induction
        mean_u_ind, mean_s_ind = self._get_induction_unspliced_spliced(
            alpha, alpha_1, lambda_alpha, beta, gamma, t_s * px_rho
        )

        if self.time_dep_transcription_rate:
            mean_u_ind_steady = (alpha_1 / beta).expand(n_cells, self.n_input)
            mean_s_ind_steady = (alpha_1 / gamma).expand(n_cells, self.n_input)
        else:
            mean_u_ind_steady = (alpha / beta).expand(n_cells, self.n_input)
            mean_s_ind_steady = (alpha / gamma).expand(n_cells, self.n_input)
        scale_u = scale[: self.n_input, :].expand(n_cells, self.n_input, 4).sqrt()

        # repression
        u_0, s_0 = self._get_induction_unspliced_spliced(
            alpha, alpha_1, lambda_alpha, beta, gamma, t_s
        )

        tau = px_tau
        mean_u_rep, mean_s_rep = self._get_repression_unspliced_spliced(
            u_0,
            s_0,
            beta,
            gamma,
            (self.t_max - t_s) * tau,
        )
        mean_u_rep_steady = torch.zeros_like(mean_u_ind)
        mean_s_rep_steady = torch.zeros_like(mean_u_ind)
        scale_s = scale[self.n_input :, :].expand(n_cells, self.n_input, 4).sqrt()

        end_penalty = ((u_0 - self.switch_unspliced).pow(2)).sum() + (
            (s_0 - self.switch_spliced).pow(2)
        ).sum()

        # unspliced
        mean_u = torch.stack(
            (
                mean_u_ind,
                mean_u_ind_steady,
                mean_u_rep,
                mean_u_rep_steady,
            ),
            dim=2,
        )
        scale_u = torch.stack(
            (
                scale_u[..., 0],
                scale_u[..., 0],
                scale_u[..., 0],
                0.1 * scale_u[..., 0],
            ),
            dim=2,
        )
        dist_u = Normal(mean_u, scale_u)
        mixture_dist_u = MixtureSameFamily(comp_dist, dist_u)

        # spliced
        mean_s = torch.stack(
            (mean_s_ind, mean_s_ind_steady, mean_s_rep, mean_s_rep_steady),
            dim=2,
        )
        scale_s = torch.stack(
            (
                scale_s[..., 0],
                scale_s[..., 0],
                scale_s[..., 0],
                0.1 * scale_s[..., 0],
            ),
            dim=2,
        )
        dist_s = Normal(mean_s, scale_s)
        mixture_dist_s = MixtureSameFamily(comp_dist, dist_s)

        return mixture_dist_s, mixture_dist_u, end_penalty

    def _get_induction_unspliced_spliced(
        self, alpha, alpha_1, lambda_alpha, beta, gamma, t, eps=1e-6
    ):
        if self.time_dep_transcription_rate:
            unspliced = alpha_1 / beta * (1 - torch.exp(-beta * t)) - (
                alpha_1 - alpha
            ) / (beta - lambda_alpha) * (
                torch.exp(-lambda_alpha * t) - torch.exp(-beta * t)
            )

            spliced = (
                alpha_1 / gamma * (1 - torch.exp(-gamma * t))
                + alpha_1
                / (gamma - beta + eps)
                * (torch.exp(-gamma * t) - torch.exp(-beta * t))
                - beta
                * (alpha_1 - alpha)
                / (beta - lambda_alpha + eps)
                / (gamma - lambda_alpha + eps)
                * (torch.exp(-lambda_alpha * t) - torch.exp(-gamma * t))
                + beta
                * (alpha_1 - alpha)
                / (beta - lambda_alpha + eps)
                / (gamma - beta + eps)
                * (torch.exp(-beta * t) - torch.exp(-gamma * t))
            )
        else:
            unspliced = (alpha / beta) * (1 - torch.exp(-beta * t))
            spliced = (alpha / gamma) * (1 - torch.exp(-gamma * t)) + (
                alpha / ((gamma - beta) + eps)
            ) * (torch.exp(-gamma * t) - torch.exp(-beta * t))

        return unspliced, spliced

    def _get_repression_unspliced_spliced(self, u_0, s_0, beta, gamma, t, eps=1e-6):
        unspliced = torch.exp(-beta * t) * u_0
        spliced = s_0 * torch.exp(-gamma * t) - (
            beta * u_0 / ((gamma - beta) + eps)
        ) * (torch.exp(-gamma * t) - torch.exp(-beta * t))
        return unspliced, spliced

    def sample(
        self,
    ) -> np.ndarray:
        """Not implemented."""
        raise NotImplementedError

    @torch.no_grad()
    def get_loadings(self) -> np.ndarray:
        """Extract per-gene weights (for each Z, shape is genes by dim(Z)) in the linear decoder."""
        # This is BW, where B is diag(b) batch norm, W is weight matrix
        if self.decoder.linear_decoder is False:
            raise ValueError("Model not trained with linear decoder")
        w = self.decoder.rho_first_decoder.fc_layers[0][0].weight
        if self.use_batch_norm_decoder:
            bn = self.decoder.rho_first_decoder.fc_layers[0][1]
            sigma = torch.sqrt(bn.running_var + bn.eps)
            gamma = bn.weight
            b = gamma / sigma
            b_identity = torch.diag(b)
            loadings = torch.matmul(b_identity, w)
        else:
            loadings = w
        loadings = loadings.detach().cpu().numpy()

        return loadings


logger = logging.getLogger(__name__)


def _softplus_inverse(x: np.ndarray) -> np.ndarray:
    x = torch.from_numpy(x)
    x_inv = torch.where(x > 20, x, x.expm1().log()).numpy()
    return x_inv


class VELOVI(VAEMixin, UnsupervisedTrainingMixin, BaseModelClass):
    """Velocity Variational Inference.

    Parameters
    ----------
    adata
        AnnData object that has been registered via :func:`~velovi.VELOVI.setup_anndata`.
    n_hidden
        Number of nodes per hidden layer.
    n_latent
        Dimensionality of the latent space.
    n_layers
        Number of hidden layers used for encoder and decoder NNs.
    dropout_rate
        Dropout rate for neural networks.
    gamma_init_data
        Initialize gamma using the data-driven technique.
    linear_decoder
        Use a linear decoder from latent space to time.
    **model_kwargs
        Keyword args for :class:`~velovi.VELOVAE`
    """

    def __init__(
        self,
        adata: AnnData,
        n_hidden: int = 256,
        n_latent: int = 10,
        n_layers: int = 1,
        dropout_rate: float = 0.1,
        gamma_init_data: bool = False,
        linear_decoder: bool = False,
        **model_kwargs,
    ):
        super().__init__(adata)
        self.n_latent = n_latent

        spliced = self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        unspliced = self.adata_manager.get_from_registry(REGISTRY_KEYS.U_KEY)

        sorted_unspliced = np.argsort(unspliced, axis=0)
        ind = int(adata.n_obs * 0.99)
        us_upper_ind = sorted_unspliced[ind:, :]

        us_upper = []
        ms_upper = []
        for i in range(len(us_upper_ind)):
            row = us_upper_ind[i]
            us_upper += [unspliced[row, np.arange(adata.n_vars)][np.newaxis, :]]
            ms_upper += [spliced[row, np.arange(adata.n_vars)][np.newaxis, :]]
        us_upper = np.median(np.concatenate(us_upper, axis=0), axis=0)
        ms_upper = np.median(np.concatenate(ms_upper, axis=0), axis=0)

        alpha_unconstr = _softplus_inverse(us_upper)
        alpha_unconstr = np.asarray(alpha_unconstr).ravel()

        alpha_1_unconstr = np.zeros(us_upper.shape).ravel()
        lambda_alpha_unconstr = np.zeros(us_upper.shape).ravel()

        if gamma_init_data:
            gamma_unconstr = np.clip(_softplus_inverse(us_upper / ms_upper), None, 10)
        else:
            gamma_unconstr = None

        self.module = VELOVAE(
            n_input=self.summary_stats["n_vars"],
            n_hidden=n_hidden,
            n_latent=n_latent,
            n_layers=n_layers,
            dropout_rate=dropout_rate,
            gamma_unconstr_init=gamma_unconstr,
            alpha_unconstr_init=alpha_unconstr,
            alpha_1_unconstr_init=alpha_1_unconstr,
            lambda_alpha_unconstr_init=lambda_alpha_unconstr,
            switch_spliced=ms_upper,
            switch_unspliced=us_upper,
            linear_decoder=linear_decoder,
            **model_kwargs,
        )
        self._model_summary_string = (
            "VELOVI Model with the following params: \nn_hidden: {}, n_latent: {}, n_layers: {}, dropout_rate: "
            "{}"
        ).format(
            n_hidden,
            n_latent,
            n_layers,
            dropout_rate,
        )
        self.init_params_ = self._get_init_params(locals())

    def train(
        self,
        max_epochs: Optional[int] = 500,
        lr: float = 1e-2,
        weight_decay: float = 1e-2,
        use_gpu: Optional[Union[str, int, bool]] = None,
        train_size: float = 0.9,
        validation_size: Optional[float] = None,
        batch_size: int = 256,
        early_stopping: bool = True,
        gradient_clip_val: float = 10,
        plan_kwargs: Optional[dict] = None,
        **trainer_kwargs,
    ):
        """Train the model.

        Parameters
        ----------
        max_epochs
            Number of passes through the dataset. If `None`, defaults to
            `np.min([round((20000 / n_cells) * 400), 400])`
        lr
            Learning rate for optimization
        weight_decay
            Weight decay for optimization
        use_gpu
            Use default GPU if available (if None or True), or index of GPU to use (if int),
            or name of GPU (if str, e.g., `'cuda:0'`), or use CPU (if False).
        train_size
            Size of training set in the range [0.0, 1.0].
        validation_size
            Size of the test set. If `None`, defaults to 1 - `train_size`. If
            `train_size + validation_size < 1`, the remaining cells belong to a test set.
        batch_size
            Minibatch size to use during training.
        early_stopping
            Perform early stopping. Additional arguments can be passed in `**kwargs`.
            See :class:`~scvi.train.Trainer` for further options.
        gradient_clip_val
            Val for gradient clipping
        plan_kwargs
            Keyword args for :class:`~scvi.train.TrainingPlan`. Keyword arguments passed to
            `train()` will overwrite values present in `plan_kwargs`, when appropriate.
        **trainer_kwargs
            Other keyword args for :class:`~scvi.train.Trainer`.
        """
        user_plan_kwargs = plan_kwargs.copy() if isinstance(plan_kwargs, dict) else {}
        plan_kwargs = {"lr": lr, "weight_decay": weight_decay, "optimizer": "AdamW"}
        plan_kwargs.update(user_plan_kwargs)

        user_train_kwargs = trainer_kwargs.copy()
        trainer_kwargs = {"gradient_clip_val": gradient_clip_val}
        trainer_kwargs.update(user_train_kwargs)

        data_splitter = DataSplitter(
            self.adata_manager,
            train_size=train_size,
            validation_size=validation_size,
            batch_size=batch_size,
            use_gpu=use_gpu,
        )
        training_plan = TrainingPlan(self.module, **plan_kwargs)

        es = "early_stopping"
        trainer_kwargs[es] = (
            early_stopping if es not in trainer_kwargs.keys() else trainer_kwargs[es]
        )
        runner = TrainRunner(
            self,
            training_plan=training_plan,
            data_splitter=data_splitter,
            max_epochs=max_epochs,
            use_gpu=use_gpu,
            **trainer_kwargs,
        )
        return runner()

    @torch.inference_mode()
    def get_state_assignment(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        gene_list: Optional[Sequence[str]] = None,
        hard_assignment: bool = False,
        n_samples: int = 20,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Tuple[Union[np.ndarray, pd.DataFrame], List[str]]:
        """Returns cells by genes by states probabilities.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        hard_assignment
            Return a hard state assignment
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray",
                    stacklevel=2,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        states = []
        for tensors in scdl:
            minibatch_samples = []
            for _ in range(n_samples):
                _, generative_outputs = self.module.forward(
                    tensors=tensors,
                    compute_loss=False,
                )
                output = generative_outputs["px_pi"]
                output = output[..., gene_mask, :]
                output = output.cpu().numpy()
                minibatch_samples.append(output)
            # samples by cells by genes by four
            states.append(np.stack(minibatch_samples, axis=0))
            if return_mean:
                states[-1] = np.mean(states[-1], axis=0)

        states = np.concatenate(states, axis=0)
        state_cats = [
            "induction",
            "induction_steady",
            "repression",
            "repression_steady",
        ]
        if hard_assignment and return_mean:
            hard_assign = states.argmax(-1)

            hard_assign = pd.DataFrame(
                data=hard_assign, index=adata.obs_names, columns=adata.var_names
            )
            for i, s in enumerate(state_cats):
                hard_assign = hard_assign.replace(i, s)

            states = hard_assign

        return states, state_cats

    @torch.inference_mode()
    def get_latent_time(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        gene_list: Optional[Sequence[str]] = None,
        time_statistic: Literal["mean", "max"] = "mean",
        n_samples: int = 1,
        n_samples_overall: Optional[int] = None,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        """Returns the cells by genes latent time.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        time_statistic
            Whether to compute expected time over states, or maximum a posteriori time over maximal
            probability state.
        n_samples
            Number of posterior samples to use for estimation.
        n_samples_overall
            Number of overall samples to return. Setting this forces n_samples=1.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray",
                    stacklevel=2,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        times = []
        for tensors in scdl:
            minibatch_samples = []
            for _ in range(n_samples):
                _, generative_outputs = self.module.forward(
                    tensors=tensors,
                    compute_loss=False,
                )
                pi = generative_outputs["px_pi"]
                ind_prob = pi[..., 0]
                steady_prob = pi[..., 1]
                rep_prob = pi[..., 2]
                # rep_steady_prob = pi[..., 3]
                switch_time = F.softplus(self.module.switch_time_unconstr)

                ind_time = generative_outputs["px_rho"] * switch_time
                rep_time = switch_time + (
                    generative_outputs["px_tau"] * (self.module.t_max - switch_time)
                )

                if time_statistic == "mean":
                    output = (
                        ind_prob * ind_time
                        + rep_prob * rep_time
                        + steady_prob * switch_time
                        # + rep_steady_prob * self.module.t_max
                    )
                else:
                    t = torch.stack(
                        [
                            ind_time,
                            switch_time.expand(ind_time.shape),
                            rep_time,
                            torch.zeros_like(ind_time),
                        ],
                        dim=2,
                    )
                    max_prob = torch.amax(pi, dim=-1)
                    max_prob = torch.stack([max_prob] * 4, dim=2)
                    max_prob_mask = pi.ge(max_prob)
                    output = (t * max_prob_mask).sum(dim=-1)

                output = output[..., gene_mask]
                output = output.cpu().numpy()
                minibatch_samples.append(output)
            # samples by cells by genes by four
            times.append(np.stack(minibatch_samples, axis=0))
            if return_mean:
                times[-1] = np.mean(times[-1], axis=0)

        if n_samples > 1:
            # The -2 axis correspond to cells.
            times = np.concatenate(times, axis=-2)
        else:
            times = np.concatenate(times, axis=0)

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                times,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return times

    @torch.inference_mode()
    def get_velocity(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        gene_list: Optional[Sequence[str]] = None,
        n_samples: int = 1,
        n_samples_overall: Optional[int] = None,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
        velo_statistic: str = "mean",
        velo_mode: Literal["spliced", "unspliced"] = "spliced",
        clip: bool = True,
    ) -> Union[np.ndarray, pd.DataFrame]:
        """Returns cells by genes velocity estimates.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        gene_list
            Return velocities for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        n_samples
            Number of posterior samples to use for estimation for each cell.
        n_samples_overall
            Number of overall samples to return. Setting this forces n_samples=1.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.
        velo_statistic
            Whether to compute expected velocity over states, or maximum a posteriori velocity over maximal
            probability state.
        velo_mode
            Compute ds/dt or du/dt.
        clip
            Clip to minus spliced value

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        if indices is None:
            indices = np.arange(adata.n_obs)
        if n_samples_overall is not None:
            indices = np.random.choice(indices, n_samples_overall)
            n_samples = 1
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray",
                    stacklevel=2,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        velos = []
        for tensors in scdl:
            minibatch_samples = []
            for _ in range(n_samples):
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    compute_loss=False,
                )
                pi = generative_outputs["px_pi"]
                alpha = inference_outputs["alpha"]
                alpha_1 = inference_outputs["alpha_1"]
                lambda_alpha = inference_outputs["lambda_alpha"]
                beta = inference_outputs["beta"]
                gamma = inference_outputs["gamma"]
                tau = generative_outputs["px_tau"]
                rho = generative_outputs["px_rho"]

                ind_prob = pi[..., 0]
                steady_prob = pi[..., 1]
                rep_prob = pi[..., 2]
                switch_time = F.softplus(self.module.switch_time_unconstr)

                ind_time = switch_time * rho
                u_0, s_0 = self.module._get_induction_unspliced_spliced(
                    alpha, alpha_1, lambda_alpha, beta, gamma, switch_time
                )
                rep_time = (self.module.t_max - switch_time) * tau
                mean_u_rep, mean_s_rep = self.module._get_repression_unspliced_spliced(
                    u_0,
                    s_0,
                    beta,
                    gamma,
                    rep_time,
                )
                if velo_mode == "spliced":
                    velo_rep = beta * mean_u_rep - gamma * mean_s_rep
                else:
                    velo_rep = -beta * mean_u_rep
                mean_u_ind, mean_s_ind = self.module._get_induction_unspliced_spliced(
                    alpha, alpha_1, lambda_alpha, beta, gamma, ind_time
                )
                if velo_mode == "spliced":
                    velo_ind = beta * mean_u_ind - gamma * mean_s_ind
                else:
                    transcription_rate = alpha_1 - (alpha_1 - alpha) * torch.exp(
                        -lambda_alpha * ind_time
                    )
                    velo_ind = transcription_rate - beta * mean_u_ind

                if velo_mode == "spliced":
                    # velo_steady = beta * u_0 - gamma * s_0
                    velo_steady = torch.zeros_like(velo_ind)
                else:
                    # velo_steady = alpha - beta * u_0
                    velo_steady = torch.zeros_like(velo_ind)

                # expectation
                if velo_statistic == "mean":
                    output = (
                        ind_prob * velo_ind
                        + rep_prob * velo_rep
                        + steady_prob * velo_steady
                    )
                # maximum
                else:
                    v = torch.stack(
                        [
                            velo_ind,
                            velo_steady.expand(velo_ind.shape),
                            velo_rep,
                            torch.zeros_like(velo_rep),
                        ],
                        dim=2,
                    )
                    max_prob = torch.amax(pi, dim=-1)
                    max_prob = torch.stack([max_prob] * 4, dim=2)
                    max_prob_mask = pi.ge(max_prob)
                    output = (v * max_prob_mask).sum(dim=-1)

                output = output[..., gene_mask]
                output = output.cpu().numpy()
                minibatch_samples.append(output)
            # samples by cells by genes
            velos.append(np.stack(minibatch_samples, axis=0))
            if return_mean:
                # mean over samples axis
                velos[-1] = np.mean(velos[-1], axis=0)

        if n_samples > 1:
            # The -2 axis correspond to cells.
            velos = np.concatenate(velos, axis=-2)
        else:
            velos = np.concatenate(velos, axis=0)

        spliced = self.adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)

        if clip:
            velos = np.clip(velos, -spliced[indices], None)

        if return_numpy is None or return_numpy is False:
            return pd.DataFrame(
                velos,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
        else:
            return velos

    @torch.inference_mode()
    def get_expression_fit(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        gene_list: Optional[Sequence[str]] = None,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
        restrict_to_latent_dim: Optional[int] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""Returns the fitted spliced and unspliced abundance (s(t) and u(t)).

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)

        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray",
                    stacklevel=2,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        fits_s = []
        fits_u = []
        for tensors in scdl:
            minibatch_samples_s = []
            minibatch_samples_u = []
            for _ in range(n_samples):
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    compute_loss=False,
                    generative_kwargs={"latent_dim": restrict_to_latent_dim},
                )

                gamma = inference_outputs["gamma"]
                beta = inference_outputs["beta"]
                alpha = inference_outputs["alpha"]
                alpha_1 = inference_outputs["alpha_1"]
                lambda_alpha = inference_outputs["lambda_alpha"]
                px_pi = generative_outputs["px_pi"]
                scale = generative_outputs["scale"]
                px_rho = generative_outputs["px_rho"]
                px_tau = generative_outputs["px_tau"]

                (mixture_dist_s, mixture_dist_u, _,) = self.module.get_px(
                    px_pi,
                    px_rho,
                    px_tau,
                    scale,
                    gamma,
                    beta,
                    alpha,
                    alpha_1,
                    lambda_alpha,
                )
                fit_s = mixture_dist_s.mean
                fit_u = mixture_dist_u.mean

                fit_s = fit_s[..., gene_mask]
                fit_s = fit_s.cpu().numpy()
                fit_u = fit_u[..., gene_mask]
                fit_u = fit_u.cpu().numpy()

                minibatch_samples_s.append(fit_s)
                minibatch_samples_u.append(fit_u)

            # samples by cells by genes
            fits_s.append(np.stack(minibatch_samples_s, axis=0))
            if return_mean:
                # mean over samples axis
                fits_s[-1] = np.mean(fits_s[-1], axis=0)
            # samples by cells by genes
            fits_u.append(np.stack(minibatch_samples_u, axis=0))
            if return_mean:
                # mean over samples axis
                fits_u[-1] = np.mean(fits_u[-1], axis=0)

        if n_samples > 1:
            # The -2 axis correspond to cells.
            fits_s = np.concatenate(fits_s, axis=-2)
            fits_u = np.concatenate(fits_u, axis=-2)
        else:
            fits_s = np.concatenate(fits_s, axis=0)
            fits_u = np.concatenate(fits_u, axis=0)

        if return_numpy is None or return_numpy is False:
            df_s = pd.DataFrame(
                fits_s,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
            df_u = pd.DataFrame(
                fits_u,
                columns=adata.var_names[gene_mask],
                index=adata.obs_names[indices],
            )
            return df_s, df_u
        else:
            return fits_s, fits_u

    @torch.inference_mode()
    def get_gene_likelihood(
        self,
        adata: Optional[AnnData] = None,
        indices: Optional[Sequence[int]] = None,
        gene_list: Optional[Sequence[str]] = None,
        n_samples: int = 1,
        batch_size: Optional[int] = None,
        return_mean: bool = True,
        return_numpy: Optional[bool] = None,
    ) -> Union[np.ndarray, pd.DataFrame]:
        r"""Returns the likelihood per gene. Higher is better.

        This is denoted as :math:`\rho_n` in the scVI paper.

        Parameters
        ----------
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.
        indices
            Indices of cells in adata to use. If `None`, all cells are used.
        transform_batch
            Batch to condition on.
            If transform_batch is:

            - None, then real observed batch is used.
            - int, then batch transform_batch is used.
        gene_list
            Return frequencies of expression for a subset of genes.
            This can save memory when working with large datasets and few genes are
            of interest.
        library_size
            Scale the expression frequencies to a common library size.
            This allows gene expression levels to be interpreted on a common scale of relevant
            magnitude. If set to `"latent"`, use the latent libary size.
        n_samples
            Number of posterior samples to use for estimation.
        batch_size
            Minibatch size for data loading into model. Defaults to `scvi.settings.batch_size`.
        return_mean
            Whether to return the mean of the samples.
        return_numpy
            Return a :class:`~numpy.ndarray` instead of a :class:`~pandas.DataFrame`. DataFrame includes
            gene names as columns. If either `n_samples=1` or `return_mean=True`, defaults to `False`.
            Otherwise, it defaults to `True`.

        Returns
        -------
        If `n_samples` > 1 and `return_mean` is False, then the shape is `(samples, cells, genes)`.
        Otherwise, shape is `(cells, genes)`. In this case, return type is :class:`~pandas.DataFrame` unless `return_numpy` is True.
        """
        adata = self._validate_anndata(adata)
        scdl = self._make_data_loader(
            adata=adata, indices=indices, batch_size=batch_size
        )

        if gene_list is None:
            gene_mask = slice(None)
        else:
            all_genes = adata.var_names
            gene_mask = [True if gene in gene_list else False for gene in all_genes]

        if n_samples > 1 and return_mean is False:
            if return_numpy is False:
                warnings.warn(
                    "return_numpy must be True if n_samples > 1 and return_mean is False, returning np.ndarray",
                    stacklevel=2,
                )
            return_numpy = True
        if indices is None:
            indices = np.arange(adata.n_obs)

        rls = []
        for tensors in scdl:
            minibatch_samples = []
            for _ in range(n_samples):
                inference_outputs, generative_outputs = self.module.forward(
                    tensors=tensors,
                    compute_loss=False,
                )
                spliced = tensors[REGISTRY_KEYS.X_KEY]
                unspliced = tensors[REGISTRY_KEYS.U_KEY]

                gamma = inference_outputs["gamma"]
                beta = inference_outputs["beta"]
                alpha = inference_outputs["alpha"]
                alpha_1 = inference_outputs["alpha_1"]
                lambda_alpha = inference_outputs["lambda_alpha"]
                px_pi = generative_outputs["px_pi"]
                scale = generative_outputs["scale"]
                px_rho = generative_outputs["px_rho"]
                px_tau = generative_outputs["px_tau"]

                (mixture_dist_s, mixture_dist_u, _,) = self.module.get_px(
                    px_pi,
                    px_rho,
                    px_tau,
                    scale,
                    gamma,
                    beta,
                    alpha,
                    alpha_1,
                    lambda_alpha,
                )
                device = gamma.device
                reconst_loss_s = -mixture_dist_s.log_prob(spliced.to(device))
                reconst_loss_u = -mixture_dist_u.log_prob(unspliced.to(device))
                output = -(reconst_loss_s + reconst_loss_u)
                output = output[..., gene_mask]
                output = output.cpu().numpy()
                minibatch_samples.append(output)
            # samples by cells by genes by four
            rls.append(np.stack(minibatch_samples, axis=0))
            if return_mean:
                rls[-1] = np.mean(rls[-1], axis=0)

        rls = np.concatenate(rls, axis=0)
        return rls

    # TODO: Add docstrings
    @torch.inference_mode()
    def get_rates(self):
        """TODO."""
        gamma, beta, alpha, alpha_1, lambda_alpha = self.module._get_rates()

        return {
            "beta": beta.cpu().numpy(),
            "gamma": gamma.cpu().numpy(),
            "alpha": alpha.cpu().numpy(),
            "alpha_1": alpha_1.cpu().numpy(),
            "lambda_alpha": lambda_alpha.cpu().numpy(),
        }

    @classmethod
    @setup_anndata_dsp.dedent
    def setup_anndata(
        cls,
        adata: AnnData,
        spliced_layer: str,
        unspliced_layer: str,
        **kwargs,
    ) -> Optional[AnnData]:
        """%(summary)s.

        Parameters
        ----------
        %(param_adata)s
        spliced_layer
            Layer in adata with spliced normalized expression
        unspliced_layer
            Layer in adata with unspliced normalized expression.

        Returns
        -------
        %(returns)s
        """
        setup_method_args = cls._get_setup_method_args(**locals())
        anndata_fields = [
            LayerField(REGISTRY_KEYS.X_KEY, spliced_layer, is_count_data=False),
            LayerField(REGISTRY_KEYS.U_KEY, unspliced_layer, is_count_data=False),
        ]
        adata_manager = AnnDataManager(
            fields=anndata_fields, setup_method_args=setup_method_args
        )
        adata_manager.register_fields(adata, **kwargs)
        cls.register_manager(adata_manager)

    # TODO: Add docstrings
    def get_directional_uncertainty(
        self,
        adata: Optional[AnnData] = None,
        n_samples: int = 50,
        gene_list: Iterable[str] = None,
        n_jobs: int = -1,
    ):
        """TODO."""
        adata = self._validate_anndata(adata)

        logger.info("Sampling from model...")
        velocities_all = self.get_velocity(
            n_samples=n_samples, return_mean=False, gene_list=gene_list
        )  # (n_samples, n_cells, n_genes)

        df, cosine_sims = _compute_directional_statistics_tensor(
            tensor=velocities_all, n_jobs=n_jobs, n_cells=adata.n_obs
        )
        df.index = adata.obs_names

        return df, cosine_sims

    def get_permutation_scores(
        self, labels_key: str, adata: Optional[AnnData] = None
    ) -> Tuple[pd.DataFrame, AnnData]:
        """Compute permutation scores.

        Parameters
        ----------
        labels_key
            Key in adata.obs encoding cell types
        adata
            AnnData object with equivalent structure to initial AnnData. If `None`, defaults to the
            AnnData object used to initialize the model.

        Returns
        -------
        Tuple of DataFrame and AnnData. DataFrame is genes by cell types with score per cell type.
        AnnData is the permutated version of the original AnnData.
        """
        adata = self._validate_anndata(adata)
        adata_manager = self.get_anndata_manager(adata)
        if labels_key not in adata.obs:
            raise ValueError(f"{labels_key} not found in adata.obs")

        # shuffle spliced then unspliced
        bdata = self._shuffle_layer_celltype(
            adata_manager, labels_key, REGISTRY_KEYS.X_KEY
        )
        bdata_manager = self.get_anndata_manager(bdata)
        bdata = self._shuffle_layer_celltype(
            bdata_manager, labels_key, REGISTRY_KEYS.U_KEY
        )
        bdata_manager = self.get_anndata_manager(bdata)

        ms_ = adata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        mu_ = adata_manager.get_from_registry(REGISTRY_KEYS.U_KEY)

        ms_p = bdata_manager.get_from_registry(REGISTRY_KEYS.X_KEY)
        mu_p = bdata_manager.get_from_registry(REGISTRY_KEYS.U_KEY)

        spliced_, unspliced_ = self.get_expression_fit(adata, n_samples=10)
        root_squared_error = np.abs(spliced_ - ms_)
        root_squared_error += np.abs(unspliced_ - mu_)

        spliced_p, unspliced_p = self.get_expression_fit(bdata, n_samples=10)
        root_squared_error_p = np.abs(spliced_p - ms_p)
        root_squared_error_p += np.abs(unspliced_p - mu_p)

        celltypes = np.unique(adata.obs[labels_key])

        dynamical_df = pd.DataFrame(
            index=adata.var_names,
            columns=celltypes,
            data=np.zeros((adata.shape[1], len(celltypes))),
        )
        N = 200
        for ct in celltypes:
            for g in adata.var_names.tolist():
                x = root_squared_error_p[g][adata.obs[labels_key] == ct]
                y = root_squared_error[g][adata.obs[labels_key] == ct]
                ratio = ttest_ind(x[:N], y[:N])[0]
                dynamical_df.loc[g, ct] = ratio

        return dynamical_df, bdata

    def _shuffle_layer_celltype(
        self, adata_manager: AnnDataManager, labels_key: str, registry_key: str
    ) -> AnnData:
        """Shuffle cells within cell types for each gene."""
        from scvi.data._constants import _SCVI_UUID_KEY

        bdata = adata_manager.adata.copy()
        labels = bdata.obs[labels_key]
        del bdata.uns[_SCVI_UUID_KEY]
        self._validate_anndata(bdata)
        bdata_manager = self.get_anndata_manager(bdata)

        # get registry info to later set data back in bdata
        # in a way that doesn't require actual knowledge of location
        unspliced = bdata_manager.get_from_registry(registry_key)
        u_registry = bdata_manager.data_registry[registry_key]
        attr_name = u_registry.attr_name
        attr_key = u_registry.attr_key

        for lab in np.unique(labels):
            mask = np.asarray(labels == lab)
            unspliced_ct = unspliced[mask].copy()
            unspliced_ct = np.apply_along_axis(
                np.random.permutation, axis=0, arr=unspliced_ct
            )
            unspliced[mask] = unspliced_ct
        # e.g., if using adata.X
        if attr_key is None:
            setattr(bdata, attr_name, unspliced)
        # e.g., if using a layer
        elif attr_key is not None:
            attribute = getattr(bdata, attr_name)
            attribute[attr_key] = unspliced
            setattr(bdata, attr_name, attribute)

        return bdata


def _compute_directional_statistics_tensor(
    tensor: np.ndarray, n_jobs: int, n_cells: int
) -> pd.DataFrame:
    df = pd.DataFrame(index=np.arange(n_cells))
    df["directional_variance"] = np.nan
    df["directional_difference"] = np.nan
    df["directional_cosine_sim_variance"] = np.nan
    df["directional_cosine_sim_difference"] = np.nan
    df["directional_cosine_sim_mean"] = np.nan
    logger.info("Computing the uncertainties...")
    results = Parallel(n_jobs=n_jobs, verbose=3)(
        delayed(_directional_statistics_per_cell)(tensor[:, cell_index, :])
        for cell_index in range(n_cells)
    )
    # cells by samples
    cosine_sims = np.stack([results[i][0] for i in range(n_cells)])
    df.loc[:, "directional_cosine_sim_variance"] = [
        results[i][1] for i in range(n_cells)
    ]
    df.loc[:, "directional_cosine_sim_difference"] = [
        results[i][2] for i in range(n_cells)
    ]
    df.loc[:, "directional_variance"] = [results[i][3] for i in range(n_cells)]
    df.loc[:, "directional_difference"] = [results[i][4] for i in range(n_cells)]
    df.loc[:, "directional_cosine_sim_mean"] = [results[i][5] for i in range(n_cells)]

    return df, cosine_sims


def _directional_statistics_per_cell(
    tensor: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Internal function for parallelization.

    Parameters
    ----------
    tensor
        Shape of samples by genes for a given cell.
    """
    n_samples = tensor.shape[0]
    # over samples axis
    mean_velocity_of_cell = tensor.mean(0)
    cosine_sims = [
        _cosine_sim(tensor[i, :], mean_velocity_of_cell) for i in range(n_samples)
    ]
    angle_samples = [np.arccos(el) for el in cosine_sims]
    return (
        cosine_sims,
        np.var(cosine_sims),
        np.percentile(cosine_sims, 95) - np.percentile(cosine_sims, 5),
        np.var(angle_samples),
        np.percentile(angle_samples, 95) - np.percentile(angle_samples, 5),
        np.mean(cosine_sims),
    )


def _centered_unit_vector(vector: np.ndarray) -> np.ndarray:
    """Returns the centered unit vector of the vector."""
    vector = vector - np.mean(vector)
    return vector / np.linalg.norm(vector)


def _cosine_sim(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """Returns cosine similarity of the vectors."""
    v1_u = _centered_unit_vector(v1)
    v2_u = _centered_unit_vector(v2)
    return np.clip(np.dot(v1_u, v2_u), -1.0, 1.0)
