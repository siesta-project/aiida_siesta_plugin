from aiida_siesta.workflows.utils.iterate_absclass import BaseIterator
from aiida.plugins import WorkflowFactory
import pytest
from aiida import orm

@pytest.fixture
def generate_iterator():

    class SubBaseIterator(BaseIterator):
        _process_class = WorkflowFactory("siesta.base")

    return SubBaseIterator


def test_validation(aiida_profile):
    """Test the validation of subclasses of `BaseIterator`."""

    class SubBaseIterator(BaseIterator):
        pass

    import pytest
    with pytest.raises(ValueError):
        SubBaseIterator()

    with pytest.raises(ValueError):      
        class SubBaseIterator(BaseIterator):
            _process_class = "w"


def test_process_input_and_parse_func(aiida_profile, generate_iterator):
    """Test of the classmethod `process_input_and_parse_fun` of `BaseIterator`."""

    import pytest
    with pytest.raises(ValueError):
        generate_iterator.process_input_and_parse_func("w")

    key, pars_func = generate_iterator.process_input_and_parse_func("structure")
    assert key == "structure"
    assert pars_func == None


def test_iterate_input_serializer(aiida_profile, generate_iterator):
    """Test of the classmethod `_iterate_input_serializer` of `BaseIterator`."""

    import pytest
    with pytest.raises(ValueError):
        generate_iterator._iterate_input_serializer({"www": "w"})

    it_over = generate_iterator._iterate_input_serializer({"www": [1,2]})
    it_ov_el = it_over.get_attribute("www")
    assert isinstance(it_ov_el, orm.List)
    assert isinstance(orm.load_node(it_ov_el[0]), orm.Int)


