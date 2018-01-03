#!/bin/bash
export C_FORCE_ROOT=1

if [ -e ~/.aiida/config.json ]
then
    echo "Aiida development environment is set up."
else
    pip install -e /code/siesta_plugin/

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

    verdi profile setdefault verdi default
    verdi profile setdefault daemon default
    verdi daemon start

    cat /code/siesta_plugin/dockerscripts/computer-setup.txt | verdi computer setup
    verdi computer configure develop
    verdi computer test develop

    cat /code/siesta_plugin/dockerscripts/code-setup.txt | verdi code setup

    echo "eval \"\$(verdi completioncommand)"\" >> ~/.bashrc
fi

find / -name \*.pyc -delete
reentry scan -r aiida

verdi daemon restart
verdi daemon logshow
