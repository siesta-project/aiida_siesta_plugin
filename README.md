# aiida_siesta_plugin
This repository contains the files that need to be added to the aiida core distribution to
implement the Siesta plugin (input, parser, and tests).
Currently one needs to copy the files to specific places in the aiida tree. When the new
plugin architecture modifications are available, this repository will hold a stand-alone
package.

*** DEVELOPMENT AND USAGE guidelines:

You have two options:

A) Use a Docker container stack (with docker-compose) which already includes
a fully configured installation of the aiida_core code, a working database node, and a virtual computer with
torque and a copy of the Siesta code.

The plugin files are automatically linked to the aiida image. They can be changed during development
and will be seen by the aiida machine automatically.

To set things up:

(Install Docker first)

./startup_firsttime.sh

This will create the three images (downloading and compiling the current version of siesta-aiida-4.0 in
a GitHub mirror) and set up the network and configuration needed. The configured computer is 'torquessh', and
the configured code is 'siesta@torquessh'.

To use the machine:

docker-compose exec --user aiida aiida /bin/bash -l

* Rough edges:

1. The scripts in examples/submission/siesta might need to be modified (along the lines of test_siesta_scf_fail.py)
to accomodate the fact that siesta@torquessh does not yet support mpi operation.

2. One needs to specify the code as an extra argument:

./test_siesta_scf_fail.py --send siesta@torquessh

3. There should be support for linking (via virtual filesystems) newer versions of the Siesta code base to make
them available to the development environment of the 'torquessh' image.

4. ... ? please tell me about anything else and we can try to abstract things away.

You can use AiiDA as usual, edit and commit plugin files, etc (in your own branch !)

To get out of the aiida machine, just type 'exit'. The docker-compose environment will keep running.
You can take it down without destroying the virtual filesystems that hold the database, the aiida repo, etc,
by issuing the command (I think !):

docker-compose down

To get it back up again:

docker-compose up -d

I do not know how to make copies of the database and repo if the stack needs to be refreshed (to compile a new
version of Siesta, for example)

* Technical notes

- The docker-compose infrastructure (see README_docker-compose.md) is at the top level. It
is no longer necessary to make an extra copy of the plugin files in the docker context,
as the latter is the top level. 

- If you configure Travis in your branch, it will use the Docker stack infrastructure provided, with a very
simple test as of yet (see .travis.yml)

B) Use your own installation of the database, and your own computers and codes.

- Download the right copy of aiida_core into the plugin directory:

git clone https:/github.com/albgar/aiida_core.git
cd aiida_core
git checkout frozen-pre-0.8

(This will checkout a branch with root at the commit which is currently the 'official' freeze until
0.8 comes out: 7b3b3c3. The branch should contain only very minor modifications over 7b3b3c3)

Go back to the plugin directory and setup symbolic links of the plugin files into the aiida_core tree:

cd ..
. setup_links.sh

Now install aiida with the '-e' option (and optionally install the dependencies for extra utilities such as
pymatgen, cif, etc, and the generation of docs)

cd aiida_core
pip install -e .
pip install -e .[docs,atomic_tools]

The first time you will need to configure AiiDA as explained in the manual.

You can use AiiDA as usual, edit and commit plugin files, etc (in your own branch !)

