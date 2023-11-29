
from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout, argv

input_pdb_file = argv[1]

# -----------------------
# added this python module to the preamble so that we can check our structure
# files for any simulation-ending issues
import pdbfixer

# added this function to do the issue checking and fixing
def fix_protein(pdb_file):
    """ Check the structure model for any fixes before sending it to be 
        parameterized and minimized.
    INPUT:
        pdb_file: string of the pdb structure file that is to be "fixed". 
    
    RETURNS:
            new_pdb_file: string for the "fixed" pdb file to be energy 
                          minimized; format: "fixed_{pdb_file}"
    """

    try:
        with open(pdb_file,'r') as pdb, open('fixed_'+pdb_file,'w') as save_file:
            fixer = pdbfixer.PDBFixer(pdbfile=pdb)
            fixer.findMissingResidues()
            fixer.findMissingAtoms()
            fixer.addMissingAtoms(seed=0)
            fixer.addMissingHydrogens()
            PDBFile.writeFile(fixer.topology,
                              fixer.positions,
                              file=save_file)
    except Exception as e:
        print(f'For {pdb_file}: fix_protein function failed with the following exception:\n{e}\n')

    return 'fixed_'+pdb_file
# -----------------------


# check to make sure your input file isn't missing any atoms or residues.
# this function should throw an error for most drastic issues but haven't tested
# for all possible issues
SimReadyStructure = fix_protein(input_pdb_file)

# read in the simulation ready structure
pdb = PDBFile(SimReadyStructure)
# use the Amber ff14SB force field and tip3p-fb water model parameters
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')
# create the OpenMM-ready system object: 
# this contains information about the bonding (topology) of your structure and
# important methodological details about how to treat the long-range 
# interactions. Also, constraining any bonds that include a hydrogen since 
# those fluctuations are too quick for our simulation time step to accurately
# capture. 
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME,
        nonbondedCutoff=1*nanometer, constraints=HBonds)
# specifying which MD integrator to use for the simulation as well as important
# details of how that integrator works. Simulation temperature, friction 
# coefficient, and MD time step size. 
integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)
# create the simulation object that brings all these things together
simulation = Simulation(pdb.topology, system, integrator)
# setting the atomic positions of the model to the positions defined in the pdb
# file
simulation.context.setPositions(pdb.positions)

# -----------------------
# grab initial energy
state = simulation.context.getState(getEnergy=True)
einit = state.getPotentialEnergy() # check which units are used by default 
print(f'Initial Energy: {einit}')
# -----------------------

# run a energy minimization calculation to remove any horrendous clashes or 
# high energy interactions
simulation.minimizeEnergy()

# -----------------------
# grab final minimization energy
state = simulation.context.getState(getEnergy=True, getPositions=True)
efinal = state.getPotentialEnergy()
print(f'Final Energy: {efinal}')
# -----------------------

# tell the simulation object that you want to write a trajectory to 'output.pdb'
simulation.reporters.append(PDBReporter('output.pdb', 1000))
# tell the simulation object that you want to get energy and temperate values 
# from the simulation on every 1000th step. 
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True,
        potentialEnergy=True, temperature=True))
# run a MD simulation for 10,000 steps; total simulation time: 40 picoseconds
simulation.step(10000)
