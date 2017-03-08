#TODO: specify a specific version here

FROM aiidateam/aiida_core_base
MAINTAINER AiiDA Team <info@aiida.net>

# Make directory (if COPY in the line below creates it automatically, it then will
# belong to root and therefore you have permission problems)
RUN mkdir /home/aiida/.dockerscripts

## Copy startup scripts, will be used later
COPY ./scripts/ /home/aiida/.dockerscripts/

# Prepare the folders - the next command will assign it to the correct
# user, so it works with proper permissions when I mount 
# an external named volume
RUN mkdir ~/.aiida/ && mkdir ~/.ssh/keys

# Make  sure all things created have the right permissions
USER root
RUN chown -R aiida:aiida /home/aiida/.dockerscripts && \
    chown -R aiida:aiida /home/aiida/.aiida && \
    chown -R aiida:aiida /home/aiida/.ssh 

# TO CHECK: also because we need to start the web service and install Apache
# expose https port (for rest? or maybe instead expose jupyter?)
#EXPOSE 443

### Here we setup the Siesta plugin for AiiDA
USER aiida
WORKDIR /home/aiida/code/
RUN mkdir siesta && cd siesta && \
    git clone https://github.com/albgar/aiida_core.git && \
    cd aiida_core/ && \
    git checkout siesta-0.7 && \
    cp -r aiida/orm/calculation/job/siesta/ ../../aiida_core/aiida/orm/calculation/job/siesta/ && \
    cp -r aiida/parsers/plugins/siesta/ ../../aiida_core/aiida/parsers/plugins/siesta/ && \
    cp -r aiida/orm/data/psf.py ../../aiida_core/aiida/orm/data/psf.py && \
    cp -r aiida/cmdline/commands/data.py ../../aiida_core/aiida/cmdline/commands/data.py

# Important to end as user root (to run the init script)
USER root

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]
