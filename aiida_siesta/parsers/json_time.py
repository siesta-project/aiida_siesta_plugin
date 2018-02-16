def get_timing_info(json_file):
    
    import json

    try:
        data = json.load(open(json_file))
        d1 = data["global_section"]["siesta"]
        t0 = d1["_time"]
        d2 = d1["IterGeom"]
        t1 = d2["state_init"]["_time"]
        t2 = d2["Setup_H0"]["_time"]
        t3 = d2["Setup_H0"]["nlefsm"]["_time"]
        t4 = d2["IterSCF"]["setup_H"]["_time"]
        t5 = d2["IterSCF"]["compute_dm"]["_time"]
        t6 = d2["PostSCF"]["_time"]
        t7 = d2["PostSCF"]["nlefsm"]["_time"]
        t8 = d1["siesta_analysis"]["_time"]

        timing_decomp = {
            "siesta": t0,
            "state_init": t1,
            "setup_H0": t2,
            "nlefsm-1": t3,
            "setup_H": t4,
            "compute_DM": t5,
            "post-SCF": t6,
            "nlefsm-2": t7,
            "siesta_analysis": t8
        }

        return t0, timing_decomp
    except:
        #
        # Either the JSON file is corrupted, or some of the fields
        # are not there. Pending a further polishing, this would happen
        # if the calculation is not complete
        #
        return None, {}




    
