# getPotentialEnergy
Quick and dirty python code to get the potential energy of a near-simulation 
ready biomolecule structure, energy minimize it, and run a brief simulation of 
it.

The majority of this code was pulled directly from: http://docs.openmm.org/development/userguide/application/02_running_sims.html#a-first-example

# requirements: 
OpenMM and pdbfixer. 

```
conda install -c conda-forge openmm
conda install pdbfixer
```

# running: 
```
python3 getPotentialEnergy.py strct.pdb
```

where strct.pdb is the near-simulation ready structure of a biomolecule. 
By near-simulation ready, I mean that the user has removed any non-standard 
residues that aren't expected in the molecular mechanics force field parameter
set. In the example provided, there is a GOL residue in the original 2jlq.pdb 
file; this residue needs to be removed since Amber ff14SB doesn't inherently 
have parameters defined for this molecule. So always do a pre-flight check of 
your models. 

This code will create a 'fixed' version that is simulation ready, print an 
initial energy value to stdout, perform an energy minimization calculation, 
print a post-minimization energy value to stdout, and then perform 10,000 steps
of MD simulation. Potential energy and temperate values will be printed to 
stdout at an interval of 1,000 steps. A trajectory file will also be written to
visualize the MD trajectory. 

## example:
The 2JLQ structure is an Apo enzyme but has a resolved glycerol molecule. Remove
that molecule and then the structure is ready for simulation!

```
python3 ../getPotentialEnergy.py 2jlq_edited.pdb > stdout
```

The stdout file should look something like this: 

```
Initial Energy: 31673.81337060561 kJ/mol
Final Energy: -88596.38975439439 kJ/mol
#"Step","Potential Energy (kJ/mole)","Temperature (K)"
1000,-65079.03037939439,295.95934187451127
2000,-65883.40537939439,304.78781083833167
3000,-66303.29600439439,303.1338326496441
4000,-66417.84287939439,302.834927493111
5000,-67124.43662939439,301.86161891683923
6000,-67130.24912939439,300.2661094687141
7000,-67010.60850439439,301.7395257499953
8000,-67604.43662939439,302.89986410250305
9000,-67491.96787939439,306.11505700176184
10000,-67815.65537939439,296.79914299225317
```

IMPORTANT NOTE: your results will not match exactly! This is because numerous
lines in the code are non-deterministic; given the same input, they will not 
report the same value over multiple trials. Simulation engines use random number
generators to provide initial velocities, placement of atoms by the "fixer" 
function may vary slightly, etc etc. This is why you should always run multiple
trails of any MD simulations/analyses rather than taking a single trial at face
value. 

