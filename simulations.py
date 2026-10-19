

import openmm as mm
from openmm.unit import *
from openmm import app
from pathlib import Path
from sys import stdout





def prepare_system(
        pdb_file: str
):
    """
    Load a protein structure and prepare it for simulation.

    Removes crystallographic water, adds hydrogens (protonation states at
    pH 7), and solvates the protein in a TIP3P-FB water box with 1.0 nm
    padding and 0.15 M NaCl (also neutralizing the net charge).

    Args:
        pdb_file: Path to the input PDB file, e.g. "1UBQ.pdb".

    Returns:
        tuple (modeller, forcefield):
            modeller: app.Modeller holding the solvated topology and positions.
            forcefield: app.ForceField (Amber ff19SB + TIP3P-FB) used to build
                the system, needed later by forcefield.createSystem().
    """
    #retrieving structure
    pdb = app.PDBFile(pdb_file)
    forcefield = app.ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, padding=1.0*nanometer, ionicStrength=0.15*molar)

    return modeller, forcefield


def setup_simulation(
        modeller: app.Modeller,
        forcefield: app.ForceField,

):
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,nonbondedCutoff=1.0*nanometer,constraints=app.HBonds) #Setup the actuall system
    system.addForce(mm.MonteCarloBarostat(1*bar,300*kelvin)) #Gives the NPT-ensemple
    integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/(picosecond),0.002*picosecond)
    simulation = app.Simulation(modeller.topology, system, integrator)
    simulation.context.setPositions(modeller.positions)
    print(f"Platform: {simulation.context.getPlatform().getName()}")
    return simulation


def run_simulation(
        simulation: app.Simulation,
        n_steps: int,
        frame_interval: int = 1000,
        output_dir: Path | None = None
):
    if output_dir is None:
        output_dir = Path("output")
    output_dir.mkdir(parents=True, exist_ok=True)
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    state = simulation.context.getState(getPositions=True)
    with open(output_dir / "start.pdb","w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    simulation.reporters.append(app.DCDReporter(str(output_dir / "traj.dcd"), frame_interval))
    simulation.reporters.append(app.StateDataReporter(
        str(output_dir / "log.csv"), frame_interval, step=True, time=True, potentialEnergy=True, temperature=True, density=True, speed=True))
    simulation.reporters.append(app.StateDataReporter(
            stdout, 5000, step=True, potentialEnergy=True, temperature=True, speed=True, progress=True, remainingTime=True, totalSteps=n_steps))
    simulation.reporters.append(app.CheckpointReporter(str(output_dir / "checkpoint.chk"),10_000))
    simulation.step(n_steps)



def main():
    test_protein = "1UBQ" #ubiquitin
    modeller, forcefield = prepare_system(f"{test_protein}.pdb")
    simulation = setup_simulation(modeller,forcefield)
    run_simulation(simulation,25_000, 500)

if __name__=="__main__":
    main()






