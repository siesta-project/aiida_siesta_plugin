#!/bin/bash
export C_FORCE_ROOT=1
find . -name \*.pyc -delete
reentry scan -r aiida

# if [ -e ~/.aiida/config.json ]
# then
#     echo "Aiida test env is set up."
# else
verdi setup --non-interactive \
        --backend=django \
        --email="aiida@localhost" \
        --db_host=db \
        --db_port=5432 \
        --db_name=aiidadb \
        --db_user=aiida \
        --db_pass=example \
        --repo=/home/aiida/.aiida_repo/ \
        --first-name="AiiDA" --last-name="In Docker" \
        --institution="None" \
        --no-password
# fi

verdi profile setdefault verdi default
verdi profile setdefault daemon default
verdi daemon start
verdi daemon logshow

