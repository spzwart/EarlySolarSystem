# EarlySolarSystem
AMUSE production script for simulating the effect of a nearby supernova explosion on the Sun's circumstellar disk.

## Usage

Get and build specific version of AMUSE and Rebound:
> get_amuse.sh

The script is run in 3 steps:
Step 1: Calculate the effect of the radiation of the supernova on the progoplanetary disk.
Step 2: Calculate the effect of the hydrodynamical blastwave on the (in step 1) irradiated disk.
Step 3: Process the data


Step 1:
%> mpiexec python -u src/supernova_irradiation.py --SNname 11aof --inc 45 --Ndisk 100000 --Nray 7 --dt_diag 3 --dt 0.5 --ff 0.01 --Rsn 0.4 > ID_Nd5Nr7Ms23Rs04pcff-2.data

Step 2:
%> mpiexec python -u src/supernova_impact_noplanets.py -f ID_Nd5Nr7Ms23Rs04pcff-2_i0244  > ID_Nd5Nr7Ms23Rs04pcff-2_i0244.data

Step 3:
%> mpiexec python -u rdc/process_last_impact_files.py -n 5 > processed.data


