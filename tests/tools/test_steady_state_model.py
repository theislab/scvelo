import scvelo as scv
from scvelo.tools import SteadyStateModel


def test_model_run():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    model = SteadyStateModel(adata)
    model.fit()
    model.state_dict()
    model.get_velocity()
