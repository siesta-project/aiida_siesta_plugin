def get_timing_info(json_file):
    
    import json

    timing_decomp = {}
    global_time = None
    
    try:
        data = json.load(open(json_file))
    except:
        #
        # The JSON file is not parseable...
        # Emit message
        return global_time, timing_decomp

    try:
        d1 = data["global_section"]["siesta"]
        global_time = d1["_time"]
        timing_decomp["siesta"] = global_time
    except:
        # wrong structure
        return global_time, timing_decomp
        
    try:
        d2 = d1["IterGeom"]
        timing_decomp["state_init"] = d2["state_init"]["_time"]
    except:
        # This might not be present in a calculation
        pass

    try:
        timing_decomp["setup_H0"] = d2["Setup_H0"]["_time"]
    except:
        # This might not be present in a calculation
        pass

    try:
        timing_decomp["nlefsm-1"] = d2["Setup_H0"]["nlefsm"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    try:
        timing_decomp["setup_H"] = d2["IterSCF"]["setup_H"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    try:
        timing_decomp["compute_DM"] = d2["IterSCF"]["compute_dm"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    try:
        timing_decomp["post-SCF"] = d2["PostSCF"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    try:
        timing_decomp["nlefsm-2"] = d2["PostSCF"]["nlefsm"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    try:
        timing_decomp["siesta_analysis"] = d1["siesta_analysis"]["_time"]
    except:
        # This might not be present in a calculation
        pass
    # Alternate name
    try:
        timing_decomp["siesta_analysis"] = d1["Analysis"]["_time"]
    except:
        # This might not be present in a calculation
        pass

    return global_time, timing_decomp




    
