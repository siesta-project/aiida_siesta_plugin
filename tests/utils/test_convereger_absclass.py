# -*- coding: utf-8 -*-
"""
Here I test only the validation checks and the errors of the one classmethod
of SequentialIterator. The rest is tested in ../test_converge.py
"""

from aiida import orm
from aiida.plugins import WorkflowFactory
import pytest

from aiida_siesta.utils.converge_absclass import SequentialConverger


@pytest.fixture
def generate_seq_converger():

    class SubSequentialConverger(SequentialConverger):
        _process_class = WorkflowFactory("siesta.converger")

    return SubSequentialConverger


def test_validation(aiida_profile):
    """Test the validation of subclasses of `SequentialConverger`."""

    with pytest.raises(ValueError):
        class SubSequentialConverger(SequentialConverger):
            _process_class = WorkflowFactory("siesta.base")

def test_iterate_input_serializer(aiida_profile, generate_seq_converger):
    """Test of the classmethod `_iterate_input_serializer` of `BaseIterator`."""

    it_ov_el = generate_seq_converger._iterate_input_serializer([{"www": [1,2]},{"eee":[2,3]}])
    assert isinstance(it_ov_el, orm.List)
    assert isinstance(orm.load_node(it_ov_el[0]), orm.Dict)
