import argparse
import MDAnalysis
import os
import logging
import sys
import glob
import yaml
import transformato.restraints as tfrs
import numpy as np
import matplotlib.pyplot as plt

logging.basicConfig(level=logging.DEBUG,stream=sys.stdout)
logger=logging.getLogger(__name__)



parser=argparse.ArgumentParser(description="Select which run you want to analyse. By defualt, visualizes distances as experienced by the restraints.")
parser.add_argument("rundir",nargs="+",help="The run (replicate) directories containing the trajectories you want to analyze")
parser.add_argument("--altrestraints",help="For systems without restraints, you may specify an alternative restraints.yaml")


args=parser.parse_args()

print(args)
rundirs=[]
for rundir in args.rundir:
    if os.path.exists(rundir):
        
        logger.info(f"Found rundir: {rundir}")
        rundirs.append(rundir)
    else:
        logger.warning(f"Could not find rundir: {rundir}. Skipping.")


    

class Intstate():
    def __init__(self,name):
        self.name=name
        self.results={"time":[]}


class ReplicateToAnalyze():
    def __init__(self,rundir):
        
        self.rundir=rundir
        self.logger=logging.getLogger(f"{rundir}")
       

        

        self.dirs_intstates=glob.glob(f"{rundir}/**/**/intst*/")
        self.num_intstates=len(self.dirs_intstates)
        
        self.logger.info(f"Processing run {self.rundir}. Found {self.num_intstates} intermediary states to process.")
        
        self.results={"time":[]}

        self.pdbpath=glob.glob(f"{rundir}/**/**/intst4/lig_in_complex.pdb")[0]
        if args.altrestraints==None:
            self.restraint_path=glob.glob(f"{rundir}/**/**/intst4/restraints.yaml")[0]
            
        else:
            self.restraint_path=glob.glob(f"{args.altrestraints}/**/**/intst4/restraints.yaml")[0]
        with open(self.restraint_path) as stream:
                
            self.config=yaml.safe_load(stream)
            self.logger.info(f"Loaded restraints configuration at {self.restraint_path}")

        self.logger.info(f"Restraints for this run were set as: {self.config['simulation']['restraints']}")
        self.restraints_to_analyze=tfrs.CreateRestraintsFromConfig(self.config,self.pdbpath)

        r_id=1
        for restraint in self.restraints_to_analyze:
            restraint.createForce(self.config["system"]["structure"]["ccs"])
            self.results[f"b_{r_id}"]=[]
            r_id+=1
        
        self.intstates=[]
        for intstatepath in self.dirs_intstates:
            self.intstates.append(self.analyze_intstate(intstatepath))
        
        subfigs=[]
        c=1
        for intstate in self.intstates:
            currentfig=plt.subplot(round(len(self.intstates)),2,c)
            c+=1    
            print (intstate.results)
            [currentfig.plot(intstate.results["time"],intstate.results[f"b_{i}"]) for i in range(1,1+len(self.restraints_to_analyze))]
            subfigs.append(currentfig)

        plt.show()
    def analyze_intstate(self,path_to_intstate,scale=100):
        
        self.universe=MDAnalysis.Universe(self.pdbpath,f"{path_to_intstate}/lig_in_complex.dcd")

        self.current_int_name=path_to_intstate.split("/")[-2]
        self.current_int_state=Intstate(self.current_int_name)
        
        self.logger.info(f"analyzing intstate {self.current_int_name}")

        r_id=1
        for restraint in self.restraints_to_analyze:
            
            
            self.current_int_state.results[f"b_{r_id}"]=[]
            r_id+=1
        for frame in self.universe.trajectory[::scale]:
            self.logger.info(self.universe.trajectory.time)
            self.current_int_state.results["time"].append(self.universe.trajectory.time)

            r_id=1
            for restraint in self.restraints_to_analyze:
                
                self.current_int_state.results[f"b_{r_id}"].append((self.get_distance(restraint.g1_in_cc,restraint.g2)))
            
                r_id+=1
        
        
        return self.current_int_state

    def get_distance(self,g1,g2):
        #self.logger.info(f"{g1}\n{g2}")
        tempg1=self.universe.atoms[g1.indices]
        tempg2=self.universe.atoms[g2.indices]
        p1=tempg1.center_of_mass()
        p2=tempg2.center_of_mass()
        vec=p1-p2
        dist=np.linalg.norm(vec)
        return dist

if __name__=="__main__":

    all_runs=[]

    for rundir in rundirs:
        all_runs.append(ReplicateToAnalyze(rundir))