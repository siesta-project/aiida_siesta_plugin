import pytest
import numpy as np
from aiida.common.exceptions import StoringNotAllowed
from aiida_siesta.data.eig import EigData


def test_setters():
    """
    Test setting up of eigs aray and e_fermi property
    """

    eig=EigData()
    
    with pytest.raises(ValueError):
        eig.e_fermi = "ww"
    eig.e_fermi = -1
    assert eig.attributes["e_fermi"]==-1

    with pytest.raises(TypeError):
        eig.set_eigs(np.array(["w","w"]))
    with pytest.raises(IndexError):
        eig.set_eigs(np.array([1,1]))    
    eig.set_eigs(np.array([[[True,True],[False,True]]]))
    eig.set_eigs(np.array([[[1.,1.],[1.,1.]]]))

    assert not eig.is_stored
    eig.store()
    assert eig.is_stored

def test_getter():
    """
    Test get_eigs
    """
    eig=EigData()
    eig.e_fermi = -1
    eig.set_eigs(np.array([[[1.,1.],[1.,1.]]]))

    assert (eig.get_eigs()==np.array([[[1.,1.],[1.,1.]]])).all()
    assert (eig.get_eigs(scale_to_ef=True)==np.array([[[2.,2.],[2.,2.]]])).all()


def test_compute_dos():
    """
    Test compute_dos
    """
    eig=EigData()
    eig.e_fermi = -1
    eig.set_eigs(np.array([[[1.,1.],[2.,4.]]]))
    
    #Even without weights, the dos is calculated with equal weight for every kpoint
    dos = eig.compute_dos()

    #Test error if number of weights is different compared to kpoints in eigs
    kp = np.array([[0.5,0.5,0.5]])
    weights = np.array([0.5])
    eig.set_kpoints(kp, weights=weights)
    with pytest.raises(ValueError):
        dos = eig.compute_dos()

    #normal run with all defaults
    kp = np.array([[0.5,0.5,0.5],[0.0,0.0,0.0]])
    weights = np.array([0.5,0.5])
    eig.set_kpoints(kp, weights=weights)
    dos = eig.compute_dos()
    assert len(dos) == 2
    dos = eig.compute_dos(scale_to_ef = True)

    #Test arguments validator
    with pytest.raises(ValueError):
        dos = eig.compute_dos(e_max="ww")
    with pytest.raises(ValueError):
        dos = eig.compute_dos(e_min="ww")
    with pytest.raises(ValueError):
        dos = eig.compute_dos(d_ene="ww")
