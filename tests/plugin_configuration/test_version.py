def test_version():
    import aiida_siesta

    from aiida_siesta.parsers.siesta import SiestaParser
    assert SiestaParser._version == aiida_siesta.__version__

    from aiida_siesta.parsers.stm import STMParser
    assert STMParser._version == aiida_siesta.__version__

    with open("aiida_siesta/docs/conf.py") as fil:
        for line in fil.readlines():
            if "version =" in line:
                s=line
    assert str(s.split()[2][2:-1]) == str(aiida_siesta.__version__)

