def test_version():
    import aiida_siesta
    import json

    with open('setup.json') as filee: 
        data = json.load(filee) 

    assert data["version"] == aiida_siesta.__version__

    from aiida_siesta.parsers.siesta import SiestaParser
    assert SiestaParser._version == aiida_siesta.__version__

    from aiida_siesta.parsers.stm import STMParser
    assert STMParser._version == aiida_siesta.__version__

