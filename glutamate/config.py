"""Central parameter definitions for the whole project.

Python port of ``InputParametersSR.m``. One :class:`Config` object is the single
source of parameters for the diffusion Monte Carlo and both receptor drivers,
mirroring the MATLAB groups ``BasicParameters, Balls, ComputationSet, Cleft,
SizeCube, Adhesion, AMPA, NMDA, Analysis``.

Units: spatial quantities in nm, time in ms, unless a field name says otherwise
(``_um``, ``_uM`` ...). Values match the current MATLAB defaults exactly.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import List


@dataclass
class SizeCube:
    Size: float = 4000.0            # nm, side length of the cubic domain
    AxisSize: float = field(init=False)  # nm, half-width (graphics + open-boundary)

    def __post_init__(self) -> None:
        self.AxisSize = self.Size / 2.0


@dataclass
class BasicParameters:
    SaveFiles: str = "y"
    SaveReceptorInputData: bool = True
    ShowGraphics: str = "n"          # 'n' by default for headless runs
    ShowObstacleBalls: bool = False
    ObstacleBallDisplayMode: str = "adhesive"
    AdhesiveBallColor: tuple = (0.90, 0.15, 0.10)
    NonAdhesiveBallColor: tuple = (0.20, 0.45, 0.85)
    GraphicsUpdateEveryNSteps: int = 1000
    ReleaseSource: str = "synapse"   # 'astrocyte' or 'synapse'
    ReleaseTimes: List[float] = field(default_factory=lambda: [0.0, 4.0])       # ms
    ReleaseFractions: List[float] = field(default_factory=lambda: [0.5, 0.5])
    AstrocyteRadius: float = 300.0
    AstrocyteCenterDistance: float = 300.0
    ReleasePatchRadius: float = 100.0
    NumberOfParticles: int = 500
    charge: float = 0.0
    TotalTime: float = 8.0           # ms
    Numbersphere: int = 1500
    UnboundProb: float = 0.0
    BeginTime: float = 0.8           # fraction of TotalTime
    EndTime: float = 0.9             # fraction of TotalTime
    Trials: int = 1
    Control: int = 2                 # 1 -> save FinalDistribution.txt
    UseExponentialDecay: bool = True
    TimeDecay: float = 8.0           # ms, exponential removal time constant


@dataclass
class Balls:
    name: str = "Balls"
    MaxBallRadius: float = 400.0
    MinBallRadius: float = 50.0
    BettaControl: float = 0.8
    PlacementMode: str = "overlapping"   # 'overlapping' or 'nonoverlapping'
    MinimumBallSeparation: float = 0.1
    MaximumPlacementAttempts: int = 100000
    OccupancySampleCount: int = 100000
    UseAdhesionDecay: bool = True
    AdhesionDecayTime: float = 1000.0    # ms
    ProbabilityofAstrocytes: float = 0.1


@dataclass
class ComputationSet:
    DifCoefINIT: float = 0.4e6           # nm^2/ms (= 0.4 um^2/ms)
    GTimeStep: float = 0.0001            # ms
    StepLocalMovement: float = 0.1
    DistanceBins: int = 100
    TimeBins: int = 100
    SnapshotTimes: List[float] = field(default_factory=lambda: [0.1, 0.3, 1.0, 3.0])
    StepRad: float = field(init=False)   # nm, 3D RMS Brownian step

    def __post_init__(self) -> None:
        self.StepRad = math.sqrt(6.0 * self.DifCoefINIT * self.GTimeStep)


@dataclass
class Cleft:
    RadiusofCleft: float = 100.0
    CleftWidth: float = 10.0
    DistanceAstrocyte: float = 10.0
    SynapticEnvelopeRadius: float = 110.0
    SpilloverRadius: float = 300.0


@dataclass
class Adhesion:
    TimeOfDeAdhesion: float = 4.0        # ms
    MaxProbAdhesive: float = 1.0
    TimeINsideAdhesiveZone: float = 1.0  # ms


@dataclass
class AMPAConfig:
    """AMPA spatial-analysis parameters (Unified_AMPA_SpaceSR)."""
    TimeGlu_ms: float = 8.0              # linked to BasicParameters.TotalTime
    BinGlu: int = 100                    # linked to ComputationSet.DistanceBins
    TotalTime_ms: float = 200.0          # ODE integration horizon
    OutputTimeStep_ms: float = 0.25
    MaxDistance_um: float = 2.0


@dataclass
class NMDAConfig:
    """NMDA spatial-analysis parameters (NMDA_SpaceSR_Optimized)."""
    TotalTime_ms: float = 200.0
    LegacyMasterTotalTime_ms: float = 25.0
    LegacyMaximumDistance_um: float = 2.0


@dataclass
class Analysis:
    """Shared receptor-analysis constants."""
    AvogadroNumber: float = 6.02214076e23
    # No default data folder: the path to a results folder is always passed
    # explicitly (CLI --data, GUI, or the driver argument), so nothing in the
    # package depends on a particular machine or directory layout.
    DataFolderName: str = ""
    ODE_RelTol: float = 1e-8
    ODE_AbsTol: float = 1e-7             # non-stiff models (AMPA, AMPA1, NMDA)
    ODE_AbsTol_Stiff: float = 1e-12      # stiff AMPA2 (BDF / Radau)


@dataclass
class Config:
    basic: BasicParameters = field(default_factory=BasicParameters)
    balls: Balls = field(default_factory=Balls)
    computation: ComputationSet = field(default_factory=ComputationSet)
    cleft: Cleft = field(default_factory=Cleft)
    cube: SizeCube = field(default_factory=SizeCube)
    adhesion: Adhesion = field(default_factory=Adhesion)
    ampa: AMPAConfig = field(default_factory=AMPAConfig)
    nmda: NMDAConfig = field(default_factory=NMDAConfig)
    analysis: Analysis = field(default_factory=Analysis)

    def __post_init__(self) -> None:
        # Link duplicated values to the diffusion parameters (as in InputParametersSR.m).
        self.ampa.TimeGlu_ms = self.basic.TotalTime
        self.ampa.BinGlu = self.computation.DistanceBins


def default_config() -> Config:
    """Return the full-scale parameter set used for production runs."""
    return Config()


def quick_config() -> Config:
    """A small configuration that runs in seconds.

    Same physics, far fewer particles, obstacles and time steps: meant for
    trying the workflow, teaching and tests, not for quantitative results.
    Matches the GUI's "Quick preset" and the CLI's --quick flag.
    """
    cfg = Config()
    cfg.basic.NumberOfParticles = 80
    cfg.basic.Numbersphere = 60
    cfg.basic.TotalTime = 0.5
    cfg.balls.OccupancySampleCount = 3000
    cfg.computation.DistanceBins = 25
    cfg.computation.TimeBins = 15
    cfg.computation.SnapshotTimes = [0.1, 0.3]
    cfg.computation.__post_init__()
    cfg.__post_init__()
    return cfg
