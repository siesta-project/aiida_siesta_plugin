# Establish symbolic links to plugin files
#
# See README.md
# This is only needed if you are using your own AiiDA installation (i.e., not using Docker)
# It should be done only once, after getting the aiida_core code base.
#
plugin_dir=$(pwd)

echo ${plugin_dir}

ln -sf ${plugin_dir}/siesta_plugin_files/input/siesta \
       ${plugin_dir}/aiida_core/aiida/orm/calculation/job

ln -sf ${plugin_dir}/siesta_plugin_files/parser/siesta \
       ${plugin_dir}/aiida_core/aiida/parsers/plugins

ln -sf ${plugin_dir}/siesta_plugin_files/examples/submission/siesta \
       ${plugin_dir}/aiida_core/examples/submission

ln -sf ${plugin_dir}/siesta_plugin_files/orm.data/psf.py \
       ${plugin_dir}/aiida_core/aiida/orm/data/psf.py 

ln -sf ${plugin_dir}/siesta_plugin_files/docs/siesta \
       ${plugin_dir}/aiida_core/docs/source/plugins

#
# Link the Siesta docs into the standard place in the hierarchy
# (Avoid doing it repeatedly)
#
insert_point=${plugin_dir}/aiida_core/docs/source/plugins/index.rst
grep -n siesta ${insert_point} || echo "   siesta/index" >>  ${insert_point}
