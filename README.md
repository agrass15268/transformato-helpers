# transformato-helpers

A variety of helper scripts/tools to use with [transformato](https://github.com/JohannesKarwou/transformato). Currently includes the following tools:

## hmr_via_parmed

Python script that uses parmed to do manual hydrogren mass repartitioning. Use if you need to manually modify the output of CHARMM-GUI after solvation.

## makereplicates

Bash script that takes transformato replicate folders suffixed with -1/ and creates a desired number of copies. Really just a fancy wrapper for cp -rf and a loop.



## metaanalysis

Bash script that reads transformato analysis output files along with information from the folder names and puts them in a handy .csv

## trajectory_analysis

Python script that analyses the distances of specified restraints during transformato trajectories. Requires a restraint.yaml, but the system does not necessarily have been run with restraints.

Able to both output a .csv with the data, as well as use matplotlib to display plots of the data. May
output either absolute or relative distances.

Planned features: RMSD data output.




For usage, all scripts include a -h option to display help.