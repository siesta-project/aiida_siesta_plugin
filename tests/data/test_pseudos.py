def test_pseudos_classmethods(generate_psf_data, generate_psml_data):

    from aiida_siesta.data.psml import PsmlData
    from aiida_siesta.data.psf import PsfData

    assert PsfData.get_psf_groups() == []
    assert PsmlData.get_psml_groups() == []

    #get_psf_group(cls, group_label)

    #from_md5

    #get_or_create

def test_pseudo(generate_psf_data, generate_psml_data):
    
    psf = generate_psf_data('Si')
    psml = generate_psml_data('Si')
    
    assert psf.get_psf_family_names() == []
    assert psml.get_psml_family_names() == []

    #set_file
    #assert 'md5' in psf.attributes
