FROM aiidateam/torquessh_base:1.0
MAINTAINER AiiDA Team <info@aiida.net>

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install required packages
RUN apt-get update \ 
    && apt-get install -y --no-install-recommends \
    gfortran \
    liblapack-dev \
    make \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all

USER app

WORKDIR /home/app

#
# Compile Siesta
#   Get source from a GitHub mirror of the Launchpad sources
#   This avoids having to install bzr and dependencies in the image
#   (temporary hack)
#
RUN curl -L -o master.tar.gz \
    https://github.com/siesta-project/siesta-aiida-4.0/archive/master.tar.gz && \
    tar xzf master.tar.gz && \
    cd siesta-aiida-4.0-master/Obj/ && \
    sh ../Src/configure FC=gfortran && \
    sh ../Src/obj_setup.sh && \
    make && \
    mkdir /home/app/bin/ && \
    cp siesta /home/app/bin/siesta && \
    cd ../../ && \
    rm -rf siesta-aiida-4.0-master master.tar.gz

# I need to go back to root
USER root   


    