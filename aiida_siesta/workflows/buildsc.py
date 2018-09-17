#!/usr/bin/env python
# Written by Victor M. Garcia-Suarez. Based on the fcbuild.f program
# of the Siesta/Util/Vibra package. July 2018

def buildsc(scarray,struct):

    try:
        lxmax=scarray.get_array('sca')[0]
    except:
        lxmax=0
    try:
        lymax=scarray.get_array('sca')[1]
    except:
        lymax=0
    try:
        lzmax=scarray.get_array('sca')[2]
    except:
        lzmax=0

    cell=[[0 for x in range(3)] for y in range(3)]
    ia=0
    for vector in struct.cell:
        cell[ia][:]=vector
        ia=ia+1
    scell=[[0 for x in range(3)] for y in range(3)]
    for i in range(3):
        scell[i][0]=(2*lxmax+1)*cell[i][0]
        scell[i][1]=(2*lymax+1)*cell[i][1]
        scell[i][2]=(2*lzmax+1)*cell[i][2]

    ia=0
    nia=len(struct.sites)
    xa=[[0 for x in range(3)] for y in range(nia)]
    spec=[0 for x in range(nia)]
    for site in struct.sites:
        xa[ia][:]=site.position[:]
        spec[ia]=site.kind_name
        ia=ia+1
    nna=nia*(2*lxmax+1)*(2*lymax+1)*(2*lzmax+1)
    xasc=[[0 for x in range(3)] for y in range(nna)]
    specsc=[0 for x in range(nna)]
    iatm=-1
    rr=[0 for x in range(3)]
    for i in range(-lxmax,lxmax+1):
        for j in range(-lymax,lymax+1):
            for k in range(-lzmax,lzmax+1):
                rr[0]=i*cell[0][0]+j*cell[0][1]+k*cell[0][2]
                rr[1]=i*cell[1][0]+j*cell[1][1]+k*cell[1][2]
                rr[2]=i*cell[2][0]+j*cell[2][1]+k*cell[2][2]
                for ia in range(nia):
                    iatm=iatm+1
                    specsc[iatm]=spec[ia]
                    for ix in range(3):
                        xasc[iatm][ix]=xa[ia][ix]+rr[ix]
    return scell, xasc, specsc

