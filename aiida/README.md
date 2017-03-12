This image starts from aiida_core-base, and adds scripts to setup aiida,
computers, codes, and run some basic tests with quantum espresso.

The set of scripts is divided into two - 'core', for the basic setup of AiiDA
and of the computer 'torquessh' managed by docker-compose, and
'plugin' for anything you want to do that is plugin-specific (in this 
example, quantum-espresso).
