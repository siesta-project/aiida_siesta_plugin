# aiida_siesta_plugin
This repository contains the files that need to be added to the aiida core distribution to
implement the Siesta plugin (input, parser, and tests).
Currently one needs to copy the files to specific places in the aiida tree. When the new
plugin architecture modifications are available, this repository will hold a stand-alone
package.

The docker-compose infrastructure (see README_docker-compose.md) is at the top level. It
is no longer necessary to make an extra copy of the plugin files in the docker context,
as the latter is the top level.

aiida_core commit needed: 7b3b3c3


