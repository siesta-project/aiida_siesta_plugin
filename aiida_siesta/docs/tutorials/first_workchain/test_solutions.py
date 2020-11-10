import aiida
from aiida.orm import Code, Str, Dict
import pytest

import deliver
from deliver import deliver_stage, stage_solution

aiida.load_profile("<profile>")

deliver.GENERAL_INPUTS = {
    "code": Code.get_from_string('<code>'),
    "pseudo_family": Str("<pseudo-family>"),
    "options": Dict(
        dict={
            'withmpi': False, 
            'max_wallclock_seconds': 3600 * 2
        }
    ),
    "parameters": Dict(),
}

@pytest.mark.parametrize("stage", [1,2,3,4])
def test_solution(stage):
    deliver_stage(stage, stage_solution(stage).deliverable)
