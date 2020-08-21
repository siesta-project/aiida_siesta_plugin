#!/usr/bin/env python
# coding: utf-8

# # Documentation
# Assuming The after supercell relaxation of pure Case :
# with "What=intial" we are going to generate images to take the initial and final image with/without ghost to relaxed with siesta and save the DMs
# 
# After running siesta we will "RERUN" the scipt with "What=relax" two generate the new relaxed images which create a neb folder ready to run.
# 
# NOTE: Have to specify correct paramters like path to the initial and final siesta results....
# 
# NOTE: Be careful of fdfs inputs.
# 
# NOTE: Double check the lua file to be consistent with image labels and iamge DMs files name.
# 
# NOTE: Now the the path defined by fractional coordinates. Have to check for cartesian one!!!.

# # Import

from __future__ import absolute_import

from aiida.engine import WorkChain, calcfunction, ToContext, while_
from aiida import orm


import sisl
import numpy as np
import os, shutil, sys
from linecache import getline
import linecache
from math import sqrt
import ase
from ase.neb import NEB

class ImageGenerationWorkChain(WorkChain):
    """
   Workchain Class for Generating Images

   """
   @classmethod
   def define(cls,spec):
       super(ImageGenerationWorkChain,cls).define(spec)

       spec.input("host_structure",
                   valid_type=orm.StructureData,
                   help="Pristine structure",
                   required=True
               )
       spec.input("initial_defect_site",
                   valid_type=orm.StructureData,
                   help="Initial Defect site position in crystal coordinates"
                    required=True
                   )
       spec.input(
            "final_defect_site",
            valid_type=orm.List,
            help="Final Defect site position in crystal coordinates",
            required=True
        )

       spec.input(
            "image_generation_type",
            valid_type=orm.Str,
            default=orm.Str("idpp")
            help="The Method for image generation Either (LI) Linear or (IDPP) Image dependetn pair potential to Apply for NEB",
               
               )
        spec.input(
            "number_of_image",
            valid_type=orm.Str,
            default=orm.Int(6)
            help="The Number of Image to Apply for Image Generatation",

        )
        spec.input(
                "using_ghost",
                valid_type=orm,Str,
                help="Either using Ghost Atoms or not (Particular to SIESTA)"
                required=True,
                )
        spec.outline(
                cls.generate_vacancy_exchange_images
                )


        #Outputs
        spec.output("images",
                    valid_type=orm.StructureData
                )

        #Checking Atom1 index
        def AtomIndex(Positions,AtomPosition,rtol,atol):
            """
             Positions=FracXYZ
             AtomPosition=InitialAtomPosition
             rtol
             atol
             Read All positions and Specific Postions coordinates to return the index number of array
     
            """   
            #Method 1
            #for i in range(len(Positions)):
            #    if np.allclose(Positions[i][0],AtomPosition[0]) and np.allclose(Positions[i][1],AtomPosition[1]) and np.allclose(Positions[i][2],AtomPosition[2]):
            #        print ("Index is",i)
            #Method 2
            rtol=1e-3
            atol=1e-4

            for i in range(len(Positions)):
                if np.isclose(Positions[i][0],AtomPosition[0],rtol,atol) and np.isclose(Positions[i][1],AtomPosition[1],rtol,atol) and np.isclose(Positions[i][2],AtomPosition[2],rtol,atol):
                #print ("Index and Atomic Position is ",i,Positions[i])
                    index=i
                 else:
                    "Print Couldn't Find the Atom"
            return(index)


        def ASEInitilizer(ASEInitialPositions):
            """
            ASE Initilizer for image
            """
            #initializing Image
            initialized=ase.Atoms()
            #print ("Initialize Image")
            for i in range(len(ASEInitialPositions)):
                #print ("initialize =",i)
                initialized+=ase.Atoms(ASEInitialPositions[i][3],
                                       positions=[(float(ASEInitialPositions[i][0]),
                                                   float(ASEInitialPositions[i][1]),
                                                   float(ASEInitialPositions[i][2]))
                                                             ])
            return(initialized)
        
        
        def generate_vacancy_exchange_images(self):

            import sisl
            
            """
            Generate vacancy exchange images
            """

            self.ase_atoms = self.inputs.host_structure.get_ase()
            self.ase_atoms.write('host.xyz')
            host=sisl.io.get_sile("host.xyz","r") 
            geomtery=host.read_geomtery()
            xyz=geomtery.xyz
            species=np.array([])
            species=np.append(species,self.ase_atoms.get_chemical_symbols())
            species=species.reshape(len(species),1)
            xyzFull=np.hstack([xyz,species])

            #AtomIndex(XYZ,InitialAtomPosition,rtol,atol)
            #AtomIndex(XYZ,FinalAtomPosition,rtol,atol)
            self.report("Removing Vacancies")
            InitialASEXYZ=np.delete(XYZFull,AtomIndex(XYZ,InitialAtomPosition,rtol,atol),0)
            FinalASEXYZ=np.delete(XYZFull,AtomIndex(XYZ,FinalAtomPosition,rtol,atol),0)
            #Saving The Ghost 
            GhostASEXYZ=XYZFull[AtomIndex(XYZ,InitialAtomPosition,rtol,atol)]
            GhostASEXYZ=np.vstack([GhostASEXYZ,XYZFull[AtomIndex(XYZ,FinalAtomPosition,rtol,atol)]])
           
            if AtomIndex(XYZ,FinalAtomPosition,rtol,atol) > AtomIndex(XYZ,InitialAtomPosition,rtol,atol):
                self.report ("Final> Initial")
                InitialASEXYZ=np.delete(InitialASEXYZ,AtomIndex(XYZ,FinalAtomPosition,rtol,atol)-1,0)
                FinalASEXYZ=np.delete(FinalASEXYZ,AtomIndex(XYZ,InitialAtomPosition,rtol,atol),0)
            if AtomIndex(XYZ,FinalAtomPosition,rtol,atol) < AtomIndex(XYZ,InitialAtomPosition,rtol,atol):
                self.report ("Initial > Final")
                InitialASEXYZ=np.delete(InitialASEXYZ,AtomIndex(XYZ,FinalAtomPosition,rtol,atol),0)
                FinalASEXYZ=np.delete(FinalASEXYZ,AtomIndex(XYZ,InitialAtomPosition,rtol,atol)-1,0)
            if (InitialASEXYZ==FinalASEXYZ).all():
                self.report ("Everthing is fine")
            else:
                self.report ("Something go wrong")
                self.report (InitialASEXYZ==FinalASEXYZ)
            self.report("Puting Back the Trace Atoms")
            InitialASEXYZ=np.vstack([InitialASEXYZ,XYZFull[AtomIndex(XYZ,InitialAtomPosition,rtol,atol)]])
            FinalASEXYZ=np.vstack([FinalASEXYZ,XYZFull[AtomIndex(XYZ,FinalAtomPosition,rtol,atol)]])

            initial=ASEInitilizer(InitialASEXYZ)
            final=ASEInitilizer(FinalASEXYZ)
            

            #%%*****************************************************************
            #%% Do the Image Creation with ASE Here
            #%%*****************************************************************
            InterpolationMethod='idpp'
            NumberOfImages=self.inputs.number_of_image
            images=[initial]
            for i in range(NumberOfImages):
                images.append(initial.copy())
                images.append(final)
            #%%
            neb=NEB(images)
            neb.interpolate(self.image_interpolation_method)
            
            # Writing xyz file for lua

            for i in range(NumberOfImages+2):
                neb.images[i].write("image-"+str(i)+".xyz")

            # In The Case of Ghost Atoms switch is on

            if self.using_ghost.upper()=="YES":
                for l in range(NumberOfImages+2):
                    with open("image-"+str(l)+".xyz",'a+') as fp:
                        fp.write("Ghost-"+str(GhostASEXYZ[0][3])+"  "+str(GhostASEXYZ[0][0])+"  "+str(GhostASEXYZ[0][1])+"  "+str(GhostASEXYZ[0][2])+"\n")
                        fp.write("Ghost-"+str(GhostASEXYZ[1][3])+"  "+str(GhostASEXYZ[1][0])+"  "+str(GhostASEXYZ[1][1])+"  "+str(GhostASEXYZ[1][2])+"\n")
        
                #    print ("yes")


            # Writing Animations for checking

            data=""
            for l in range(NumberOfImages+2):
                with open("image-"+str(l)+".xyz",'r') as fp:
                    data += fp.read()
            with open ('Animatiox.xyz', 'w') as fp:
                fp.write(data)








