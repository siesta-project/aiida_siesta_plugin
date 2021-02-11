import numpy as np
from aiida.orm import ArrayData

def get_epsilon_from_eps2(eps2_array):

    """
    Computes epsilon(0) from eps2(w) using Kramers-Kronig relations. 
    For this particular case:
       epsilon(0) = 1 + 2/pi int_0^\infinity { eps2(w)/w dw}

    :param: eps2_array: An ArrayData object with a 2D array (e and eps2(e))
                        Note that units do not matter, so e or w can be used.
   
    :return: epsilon(0)
    """

    array = eps2_array.get_array('e_eps2')

    # We assume that the initial energy is 0 in the calculation of eps2,
    # or at least that the start of the range of energies is below the absorption edge,
    # and is non-negative...
    #
    e = array[:,0] + 0.00001   # to regularize at the origin
    delta = e[-1] - e[-2]      # grid separation of energy points
    eps2 = array[:,1]

    epsilon = 1 + (2.0 / np.pi) * delta * sum(eps2/e)

    epsilon = round(epsilon,3)
    return epsilon
