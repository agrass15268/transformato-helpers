import argparse
import parmed as pm
import logging
import os
import glob

parser=argparse.ArgumentParser(description="This program is designed to do manual Hydrogen Mass Repartitioning instead of CHARMM - GUI, using parmed")
parser.add_argument("inputdir",metavar="inpfold",help="Inputfolder ([...]/openmm/")
parser.add_argument("--inputpsf","-i",metavar="input.psf",help="The input.psf to repartition. Default is step3.1_omm.psf",default="step3.1_omm.psf")
parser.add_argument("--output","-o",metavar="output.psf",help="The repartitioned output .psf. Default is step3_input.psf",default="step3_input.psf")

args=parser.parse_args()

print(args)
os.chdir(args.inputdir)
print(f"Changed Working directory to {os.getcwd()}")
print(f"Reading psf file: {args.inputpsf}")
psf=pm.charmm.CharmmPsfFile(f"{args.inputpsf}")
print("Reading toppar.str")
#""" with open(f"{args.inputdir}/toppar.str","r") as topparfile:
#    topparlines=topparfile.readlines()
#    topparlines=[line.strip("\n") for line in topparlines]
#    topparlines=[f"{args.inputdir}/{line}" for line in topparlines]
#    topparlines[:]=[x for x in topparlines if "." in x] # remove empty lines
#
#print(topparlines) """

parfiles=[x for x in glob.glob("../**/*") if ".rtf" in x or ".par" in x or ".str" in x or ".prm" in x]
print(parfiles)
parmset=pm.charmm.CharmmParameterSet(*parfiles)

psf.load_parameters(parmset)
print(f"PSF read successfully. Performing HMR.")
hmr=pm.tools.actions.HMassRepartition(psf)
hmr.execute()
print(f"HMR successfull. Writing to {args.output}")
psf.save(f"{args.output}",overwrite=True)