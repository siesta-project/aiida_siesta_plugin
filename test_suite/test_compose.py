#!/bin/bash

# Stop the script as soon as one step fails
set -e

## I can now setup the computer etc.
docker-compose exec --user aiida aiida /bin/bash -l -c "verdi run /home/aiida/.dockerscripts/plugin/test_plugin.py"



