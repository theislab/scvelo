import numpy as np
from scipy.integrate import odeint

from scvelo.tools.dynamical_model_utils import mRNA


def mRNA_rhs(y, t, alpha, beta, gamma):
    """Right-hand side of the mRNA splicing ODE."""
    u, s = y
    du = alpha - beta * u
    ds = beta * u - gamma * s
    return du, ds


def mRNA_num(tau, u0, s0, alpha, beta, gamma):
    """Numerical solution of the mRNA splicing ODE."""
    sol = odeint(
        func=mRNA_rhs,
        y0=np.array([u0, s0]),
        t=np.array([0, tau]),
        args=(alpha, beta, gamma),
    )
    return sol[1, :]


def test_analytical_solution():
    """
    Test whether the analytical solution or the mRNA splicing ODE is close to
    a numerical one.
    """
    # fix parameter values
    tau, u0, s0 = 0.5, 1, 0
    alpha, beta, gamma = 0.5, 0.4, 0.3

    # obtain analytical and numerical solution
    u_ana, s_ana = mRNA(tau, u0, s0, alpha, beta, gamma)
    u_num, s_num = mRNA_num(tau, u0, s0, alpha, beta, gamma)

    # compare
    assert np.allclose([u_ana, s_ana], [u_num, s_num])
