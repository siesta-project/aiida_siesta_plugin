# Establish symbolic links to plugin files
#
# See README.md
# This is only needed if you are using your own AiiDA installation (i.e., not using Docker)
# It should be done only once, after getting the aiida_core code base.
#
plugin_dir=$(pwd)

ln -sf ${plugin_dir}/siesta_plugin_files/input/siesta \
       ${plugin_dir}/aiida_core/aiida/orm/calculation/job

ln -sf ${plugin_dir}/siesta_plugin_files/parser/siesta \
       ${plugin_dir}/aiida_core/aiida/parsers/plugins/siesta

ln -sf ${plugin_dir}/siesta_plugin_files/examples/submission/siesta \
       ${plugin_dir}/aiida_core/examples/submission/siesta

ln -sf ${plugin_dir}/siesta_plugin_files/orm.data/psf.py \
       ${plugin_dir}/aiida_core/aiida/orm/data/psf.py 

