import scvelo as scv
from scvelo.tools import SecondOrderSteadyStateModel, SteadyStateModel


# TODO: Add test to check equivalence with original implementation
def test_first_order_model_run():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    model = SteadyStateModel(adata)
    model.fit()
    model.state_dict()
    model.get_velocity()


# TODO: Add test to check equivalence with original implementation
def test_second_order_model_run():
    adata = scv.datasets.simulation(random_seed=0, n_vars=10)
    scv.pp.filter_and_normalize(adata)
    scv.pp.moments(adata)
    model = SecondOrderSteadyStateModel(adata)
    model.fit()
    model.state_dict()
    model.get_velocity()
