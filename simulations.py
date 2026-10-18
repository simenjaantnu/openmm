import openmm as mm
from openmm.unit import *
from openmm import app
from dataclasses import dataclass, field
from pathlib import Path
from sys import stdout

@dataclass
class SimulationConfig:
    """Parameters for a short, geometry-preserving metalloprotein simulation."""

    coordination_cutoff: float = 2.8  # Angstrom
    metal_ligand_k: float = 200.0  # kcal/mol/A^2
    ligand_ligand_k: float = 50.0  # kcal/mol/A^2
    temperature: float = 300.0  # K
    friction: float = 1.0  # 1/ps
    timestep: float = 0.002  # ps
    equilibration_steps: int = 5_000  # 10 ps
    production_steps: int = 25_000  # 50 ps
    frame_interval: int = 50  # 1 ps -> 100 frames
    solvent_padding: float = 1.0  # nm
    ionic_strength: float = 0.15  # molar
    ph: float = 7.0
    explicit_solvent: bool = True
    platform: str | None = None
    output_dir: Path 


def prepare_system(
        pdb_file: str
):
    #retrieving structure
    pdb = app.PDBFile(pdb_file)
    forcefield = app.ForceField('amber19-all.xml', 'amber19/tip3pfb.xml')

    modeller = app.Modeller(pdb.topology, pdb.positions)
    modeller.deleteWater()
    modeller.addHydrogens(forcefield)
    modeller.addSolvent(forcefield, padding=1.0*nano*meter, ionicStrength=0.15*molar)

    return modeller


def setup_simulation(
        modeller: app.Modeller,
        forcefield: app.ForceField,

):
    system = forcefield.createSystem(modeller.topology, nonbondedMethod=app.PME,nonbondedCutoff=1.0*nano*meter,constraints=app.HBonds) #Setup the actuall system
    system.addForce(mm.MonteCarloBarostat(1*unit.bar,300*unit.kelvin)) #Gives the NPT-ensemple
    integrator = mm.LangevinMiddleIntegrator(300*kelvin, 1/(pico*second),0.004*pico*second)
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
    simulation.minimizeEnergy()
    simulation.context.setVelocitiesToTemperature(300*kelvin)
    state = simulation.context.getState(getPositions=True)
    with open(output_dir / "start.pdb","w") as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    simulation.reporters.append(app.DCDReporter("traj.dcd", frame_interval))
    simulation.reporters.append(app.StateDataReporter(
        stdout, frame_interval, step=True, time=True, potentialEnergy=True, temperature=True, density=True, speed=True))
    
    







