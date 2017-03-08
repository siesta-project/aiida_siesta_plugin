FROM aiidateam/torquessh_base:1.0
MAINTAINER AiiDA Team <info@aiida.net>

# Use baseimage-docker's init system.
CMD ["/sbin/my_init"]

# Install required packages
RUN apt-get update \ 
    && apt-get install -y --no-install-recommends \
    bzr \
    gfortran \
    liblapack-dev \
    liblapack-dev \
    make \
    && rm -rf /var/lib/apt/lists/* \
    && apt-get clean all

USER app

WORKDIR /home/app

RUN bzr branch https://code.launchpad.net/~albertog/siesta/4.0-aiida && \
    cd 4.0-aiida/Obj/ && \
    sh ../Src/configure FC=gfortran && \
    sh ../Src/obj_setup.sh && \
    make && \
    mkdir /home/app/bin/ && \
    cp siesta /home/app/bin/siesta && \
    cd ../../ && \
    rm -rf 4.0-aiida

# I need to go back to root
USER root   


    