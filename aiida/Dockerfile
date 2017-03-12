#TODO: specify a specific version here

FROM aiidateam/aiida_core_base
MAINTAINER AiiDA Team <info@aiida.net>

# Make directory (if COPY in the line below creates it automatically, it then will
# belong to root and therefore you have permission problems)
RUN mkdir /home/aiida/.dockerscripts

## Copy startup scripts, will be used later
COPY ./aiida/scripts/ /home/aiida/.dockerscripts/

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
COPY ./siesta_plugin_files/input/siesta/ /home/aiida/code/aiida_core/aiida/orm/calculation/job/siesta/
COPY ./siesta_plugin_files/parser/siesta/ /home/aiida/code/aiida_core/aiida/parsers/plugins/siesta/
COPY ./siesta_plugin_files/data/psf.py /home/aiida/code/aiida_core/aiida/orm/data/psf.py 

# Important to end as user root (to run the init script)
# (AG: Also, explicitly set the owner of the plugin files... I am not sure what is different from
#  earlier versions)
#
USER root
RUN chown -R aiida:aiida /home/aiida/code/aiida_core/aiida/orm/calculation/job/siesta && \
    chown -R aiida:aiida /home/aiida/code/aiida_core/aiida/parsers/plugins/siesta && \
    chown    aiida:aiida /home/aiida/code/aiida_core/aiida/orm/data/psf.py
    
# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]
