import openmm as mm
from dataclasses import dataclass, field
from pathlib import Path

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
    output_dir: Path = field(
        default_factory=lambda: config.directory.structures.parent / "simulations"
    )