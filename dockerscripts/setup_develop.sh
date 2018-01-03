#!/bin/bash
export C_FORCE_ROOT=1  # environment variable to be able to run
                       # aiida daemon under root inside the container

# Check if the container is set up.
# That means that it is re-used by docker-compose.
if [ -e ~/.aiida/config.json ]
then
    echo "Aiida development environment is set up."
else
    # ..if not, let's configure it!
    # First, install siesta plugin from the mounted volume:
    pip install -e /code/siesta_plugin/

    # Setup default aiida profile:
    verdi setup --non-interactive \
          --backend=django \
          --email="aiida@localhost" \
          --db_host=db \
          --db_port=5432 \
          --db_name=aiidadb \
          --db_user=aiida \
          --db_pass=aiidapwd \
          --repo=/root/.aiida_repo/ \
          --first-name="AiiDA" --last-name="In Docker" \
          --institution="None" \
          --no-password

    # Set default profile for `verdi` command and the daemon;
    # Then start the daemon:
    verdi profile setdefault verdi default
    verdi profile setdefault daemon default
    verdi daemon start

    # Setup the container's localhost as the test computer for aiida.
    # It's label will be called `develop`:
    cat /code/siesta_plugin/dockerscripts/computer-setup.txt | verdi computer setup
    verdi computer configure develop
    verdi computer test develop

    # Setup the container's SIESTA code executable as the test code for aiida.
    # It will be labeled `siesta@develop`:
    cat /code/siesta_plugin/dockerscripts/code-setup.txt | verdi code setup

    # Enable verdi auto-completion inside the container:
    echo "eval \"\$(verdi completioncommand)"\" >> ~/.bashrc
fi

find / -name \*.pyc -delete  # remove python bytecode from everywhere
                             # (some of it could point to ghost places)

reentry scan -r aiida  # to make sure that `reentry` finds everything

verdi daemon restart   # restart the daemon
verdi daemon logshow   # eventually, launch the daemon continuous log show
                       # as the `service` to prevent shutdown of the dev container
