from __future__ import absolute_import

class ProtocolRegistry:
    _protocol_registy = {
            "standard" : { 
                "description": "blabla",
                "parameters": {
                    'max-scfiterations': 50,
                    'dm-numberpulay': 4,
                    'dm-mixingweight': 0.3,
                    'dm-tolerance': 1.e-3,
                    'solution-method': 'diagon',
                    'electronic-temperature': '25 meV',
                    'write-forces': True,
                    'min_meshcut': 80
                    },
                "atomic_heuristics": {
                    'H': {'cutoff': 100 },
                    'Si': {'cutoff': 101 }
                    },
                "basis": {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZ'
                    },
                "kpoints" : {"distance" : 0.3 ,"offset" : [0., 0., 0.]},
                "pseudo_family" : "sample_psf_family",#nc-fr-04_pbe_standard_psml",
                "threshold_forces" : 0.04, #in ev/ang
                "threshold_stress" : 0.006241509125883258, #in ev/ang**3 = 1 GPa
                },
                   
            "stringent" : { 
                "description": "blublu",
                "parameters": {
                    'max-scfiterations': 50,
                    'dm-numberpulay': 4,
                    'dm-mixingweight': 0.3,
                    'dm-tolerance': 1.e-4,
                    'solution-method': 'diagon',
                    'electronic-temperature': '25 meV',
                    'write-forces': True,
                    'min_meshcut': 100
                    },
                "atomic_heuristics": {
                    'H': {'cutoff': 101 },
                    'Si': {'cutoff': 103 }
                    },
                "basis": {
                    'pao-energy-shift': '100 meV',
                    'pao-basis-size': 'DZP'
                    },
                "kpoints" : {"distance" : 0.2 ,"offset" : [0., 0., 0.]},
                "pseudo_family" : "nc-fr-04_pbe_stringent_psml",
                "threshold_forces" : 0.03, #in ev/ang
                "threshold_stress" : 0.005, #in ev/ang**3
                },
    }


    @classmethod
    def get_protocol_names(cls):
        return cls._protocol_registy.keys()

    @classmethod
    def get_default_protocol_name(cls):
        return "standard"

    @classmethod
    def get_protocol_info(cls, key):
        if key in cls._protocol_registy:
            return cls._protocol_registy[key]["description"]
        else:
            raise ValueError("Wrong protocol: no protocol with name {} implemented".format(key))

    #This maybe should become private!
    @classmethod
    def get_protocol(cls, key):
        if key in cls._protocol_registy:
            return cls._protocol_registy[key]
        else:
            raise ValueError("Wrong protocol: no protocol with name {} implemented".format(key))

