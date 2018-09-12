from ..logging import logg, settings
from .transition_matrix import transition_matrix
from .utils import scale
from scipy.sparse import linalg
import numpy as np


def eigs(T, k=10, eps=1e-3):                            # cred to M. Lange
    eigval, eigvec = linalg.eigs(T.T, k=k, which='LR')  # find k eigs with largest real part
    p = np.argsort(eigval)[::-1]                        # sort in descending order of eigenvalues
    eigval = eigval.real[p]
    eigvec = np.absolute(eigvec.real[:, p])

    eigvec = eigvec[:, (eigval >= 1 - eps)]             # select eigenvectors with eigenvalue of 1
    upper_bound = np.percentile(eigvec.flatten(), 99)   # clip to 99% quantile
    eigvec = np.clip(eigvec, 0, upper_bound)
    return eigval, eigvec


def terminal_states(adata, vkey='velocity', basis=None, weight_diffusion=0, scale_diffusion=1, eps=1e-3, copy=False):
    logg.info('computing root cells', r=True, end=' ')
    T = transition_matrix(adata, vkey=vkey, basis=basis,
                          weight_diffusion=weight_diffusion, scale_diffusion=scale_diffusion, backward=True)
    _, eigvec = eigs(T, eps=eps)
    adata.obs['root'] = scale(eigvec.sum(1))
    logg.info('using ' + str(eigvec.shape[1]) + ' eigenvectors with eigenvalue 1.')

    logg.info('computing end points', end=' ')
    T = transition_matrix(adata, vkey=vkey,
                          weight_diffusion=weight_diffusion, scale_diffusion=scale_diffusion, basis=basis)
    _, eigvec = eigs(T, eps=eps)
    adata.obs['end'] = scale(eigvec.sum(1))
    logg.info('using ' + str(eigvec.shape[1]) + ' eigenvectors with eigenvalue 1.')

    logg.info('    finished', time=True, end=' ' if settings.verbosity > 2 else '\n')
    logg.hint(
        'added to `.obs`\n'
        '    \'root\', root cells of Markov diffusion process\n'
        '    \'end\', end points of Markov diffusion process')
    return adata if copy else None
