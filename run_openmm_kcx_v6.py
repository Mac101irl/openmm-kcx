#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX (Carboxylated Lysine) Support - v6
Tamarind Platform Compatible Version

Conventions:
- File inputs: Read from inputs/settingName.ext
- Job Settings: Access via os.getenv('settingName')
- Outputs: Save ALL results to out/ directory
- Job name: Available as os.getenv('JobName')

Author: Generated for BMS Hydantoinase Project
Version: 6.0 (Tamarind Compatible)
"""

import os
import sys
import json
import subprocess
import shutil
import glob
from pathlib import Path

# =============================================================================
# SETTINGS - Read from Environment Variables (Tamarind Convention)
# =============================================================================

def get_settings():
    """Get all settings from environment variables"""
    
    def get_env(name, default=None, type_fn=str):
        """Get environment variable with type conversion"""
        val = os.getenv(name, default)
        if val is None:
            return None
        if type_fn == bool:
            return str(val).lower() in ('true', 'yes', '1')
        try:
            return type_fn(val)
        except (ValueError, TypeError):
            return default
    
    settings = {
        # Job info
        'jobName': get_env('JobName', 'openmm_simulation'),
        
        # Input files (these come from inputs/ directory)
        'pdbFile': 'inputs/pdbFile.pdb',  # Tamarind convention
        'ligandFile': 'inputs/ligandFile.sdf' if os.path.exists('inputs/ligandFile.sdf') else None,
        'ligandCharge': get_env('ligandCharge', 0, int),
        
        # Force field
        'forceField': get_env('forceField', 'ff19SB'),
        'waterModel': get_env('waterModel', 'tip3p'),
        
        # Solvation
        'boxSize': get_env('boxSize', 12.0, float),
        'boxShape': get_env('boxShape', 'cube'),
        'ionicStrength': get_env('ionicStrength', 0.15, float),
        'positiveIon': get_env('positiveIon', 'Na+'),
        'negativeIon': get_env('negativeIon', 'Cl-'),
        'pH': get_env('pH', 7.0, float),
        'removeWaters': get_env('removeWaters', 'no'),
        
        # Nonbonded
        'nonbondedMethod': get_env('nonbondedMethod', 'PME'),
        'nonbondedCutoff': get_env('nonbondedCutoff', 1.0, float),
        'ewaldErrorTolerance': get_env('ewaldErrorTolerance', 0.0005, float),
        'switchDistance': get_env('switchDistance', None, float),
        
        # Constraints
        'constraints': get_env('constraints', 'HBonds'),
        'rigidWater': get_env('rigidWater', 'yes'),
        'hydrogenMass': get_env('hydrogenMass', None, float),
        
        # Simulation
        'minimizationSteps': get_env('minimizationSteps', 10000, int),
        'minimizationTolerance': get_env('minimizationTolerance', 10.0, float),
        'equilibrationTime': get_env('equilibrationTime', 0.2, float),
        'productionTime': get_env('productionTime', 5.0, float),
        'timestep': get_env('timestep', 2.0, float),
        
        # Temperature & Pressure
        'temperature': get_env('temperature', 310.0, float),
        'pressure': get_env('pressure', 1.0, float),
        'frictionCoeff': get_env('frictionCoeff', 1.0, float),
        'barostatInterval': get_env('barostatInterval', 25, int),
        
        # Integrator
        'integrator': get_env('integrator', 'LangevinMiddle'),
        
        # Output
        'equilTrajFreq': get_env('equilTrajFreq', 1000, int),
        'prodTrajFreq': get_env('prodTrajFreq', 1000, int),
        'checkpointFreq': get_env('checkpointFreq', 10000, int),
        'trajectoryFormat': get_env('trajectoryFormat', 'DCD'),
        
        # Analysis
        'stepSize': get_env('stepSize', 5, int),
        'rmsdMask': get_env('rmsdMask', 'backbone'),
    }
    
    # Check for ligand file with different extensions
    for ext in ['.sdf', '.mol2', '.pdb']:
        ligand_path = f'inputs/ligandFile{ext}'
        if os.path.exists(ligand_path):
            settings['ligandFile'] = ligand_path
            break
    
    return settings


class Settings:
    """Settings object with attribute access"""
    def __init__(self, settings_dict):
        for key, value in settings_dict.items():
            setattr(self, key, value)


# =============================================================================
# KCX PARAMETER FILES
# =============================================================================

def create_kcx_frcmod():
    """Create AMBER frcmod file for KCX (carboxylated lysine) parameters"""
    return """KCX (N6-carboxylysine) parameters for AMBER ff19SB
Remark: Parameters for carboxylated lysine found in metalloenzymes
MASS
c3 12.01         0.878   sp3 carbon (GAFF2)
c  12.01         0.616   sp2 carbonyl carbon
o  16.00         0.434   carbonyl oxygen
n  14.01         0.530   sp2 nitrogen with 1 H

BOND
c3-n   330.6    1.456   from GAFF2
c -n   427.6    1.379   from GAFF2
c -o   637.7    1.218   from GAFF2
n -hn   403.2    1.013   from GAFF2

ANGLE
c3-c3-n    66.0      111.04  from GAFF2
c3-n -c    63.9      121.35  from GAFF2
c3-n -hn   46.0      116.78  from GAFF2
c -n -hn   48.3      117.55  from GAFF2
n -c -o    75.8      122.03  from GAFF2
o -c -o    78.2      129.52  carboxylate

DIHE
c3-c3-n -c     1    0.650       180.0     2.    from GAFF2
c3-c3-n -hn    1    0.000         0.0     2.    from GAFF2
c3-n -c -o     1    2.500       180.0     2.    from GAFF2
hn-n -c -o     1    2.500       180.0    -2.    from GAFF2
hn-n -c -o     1    2.000         0.0     1.    from GAFF2
X -c -n -X     4   10.000       180.0     2.    general

IMPROPER
c -n -o -o     1.1         180.         2.    carboxylate planarity
n -c3-c -hn   1.1         180.         2.    nitrogen planarity

NONBON
c3    1.9080  0.1094    sp3 carbon GAFF2
c     1.9080  0.0860    sp2 carbon
o     1.6612  0.2100    carbonyl oxygen
n     1.8240  0.1700    sp2 nitrogen
hn    0.6000  0.0157    H on nitrogen
"""


def create_kcx_lib():
    """Create AMBER library entry for KCX residue"""
    return """!!index array str
 "KCX"
!entry.KCX.unit.atoms table  str name  str type  int typex  int reession  int flags  int seq  int elmnt  dbl chg
 "N" "N" 0 1 131072 1 7 -0.4157
 "H" "H" 0 1 131072 2 1 0.2719
 "CA" "CX" 0 1 131072 3 6 -0.0581
 "HA" "H1" 0 1 131072 4 1 0.1360
 "CB" "2C" 0 1 131072 5 6 -0.0070
 "HB2" "HC" 0 1 131072 6 1 0.0367
 "HB3" "HC" 0 1 131072 7 1 0.0367
 "CG" "2C" 0 1 131072 8 6 0.0104
 "HG2" "HC" 0 1 131072 9 1 0.0103
 "HG3" "HC" 0 1 131072 10 1 0.0103
 "CD" "2C" 0 1 131072 11 6 -0.0377
 "HD2" "HC" 0 1 131072 12 1 0.0621
 "HD3" "HC" 0 1 131072 13 1 0.0621
 "CE" "2C" 0 1 131072 14 6 0.1313
 "HE2" "HP" 0 1 131072 15 1 0.1135
 "HE3" "HP" 0 1 131072 16 1 0.1135
 "NZ" "n" 0 1 131072 17 7 -0.5163
 "HZ" "hn" 0 1 131072 18 1 0.3339
 "CX" "c" 0 1 131072 19 6 0.7231
 "OQ1" "o" 0 1 131072 20 8 -0.7855
 "OQ2" "o" 0 1 131072 21 8 -0.7855
 "C" "C" 0 1 131072 22 6 0.5973
 "O" "O" 0 1 131072 23 8 -0.5679
!entry.KCX.unit.connectivity table  int atom1x  int atom2x  int flags
 1 2 1
 1 3 1
 3 4 1
 3 5 1
 3 22 1
 5 6 1
 5 7 1
 5 8 1
 8 9 1
 8 10 1
 8 11 1
 11 12 1
 11 13 1
 11 14 1
 14 15 1
 14 16 1
 14 17 1
 17 18 1
 17 19 1
 19 20 2
 19 21 1
 22 23 2
!entry.KCX.unit.residueconnect table  int c1x  int c2x  int c3x  int c4x  int c5x  int c6x
 1 22 0 0 0 0
!entry.KCX.unit.name single str
 "KCX"
"""


def write_kcx_parameters(output_dir):
    """Write KCX parameter files to output directory"""
    os.makedirs(output_dir, exist_ok=True)
    
    frcmod_path = os.path.join(output_dir, 'kcx.frcmod')
    lib_path = os.path.join(output_dir, 'kcx.lib')
    
    with open(frcmod_path, 'w') as f:
        f.write(create_kcx_frcmod())
    
    with open(lib_path, 'w') as f:
        f.write(create_kcx_lib())
    
    return frcmod_path, lib_path


# =============================================================================
# LIGAND PREPARATION
# =============================================================================

def prepare_ligand(ligand_file, ligand_charge, output_dir='prep'):
    """Prepare ligand with antechamber using GAFF2 and AM1-BCC charges"""
    print(f"\n{'='*60}")
    print(f"LIGAND PREPARATION")
    print(f"{'='*60}")
    print(f"  Ligand file: {ligand_file}")
    print(f"  Net charge: {ligand_charge}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Determine input file type
    ext = os.path.splitext(ligand_file)[1].lower()
    if ext == '.sdf':
        input_type = 'sdf'
    elif ext in ['.mol2']:
        input_type = 'mol2'
    elif ext == '.pdb':
        input_type = 'pdb'
    else:
        input_type = 'sdf'
    
    mol2_out = os.path.join(output_dir, 'ligand.mol2')
    frcmod_out = os.path.join(output_dir, 'ligand.frcmod')
    
    # Run antechamber for AM1-BCC charges
    print("  Running antechamber for AM1-BCC charges...")
    cmd_antechamber = [
        'antechamber',
        '-i', ligand_file,
        '-fi', input_type,
        '-o', mol2_out,
        '-fo', 'mol2',
        '-c', 'bcc',
        '-at', 'gaff2',
        '-nc', str(ligand_charge),
        '-rn', 'LIG',
        '-pf', 'y'
    ]
    
    result = subprocess.run(cmd_antechamber, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(mol2_out):
        print(f"  WARNING: Antechamber failed: {result.stderr}")
        print("  Trying with sqm charges...")
        cmd_antechamber[cmd_antechamber.index('bcc')] = 'gas'
        result = subprocess.run(cmd_antechamber, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: Antechamber failed completely")
            return None, None
    
    # Run parmchk2 for missing parameters
    print("  Running parmchk2 for missing parameters...")
    cmd_parmchk = [
        'parmchk2',
        '-i', mol2_out,
        '-f', 'mol2',
        '-o', frcmod_out,
        '-s', 'gaff2'
    ]
    
    result = subprocess.run(cmd_parmchk, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  WARNING: parmchk2 warning: {result.stderr}")
    
    print(f"  Ligand prepared: {mol2_out}")
    return mol2_out, frcmod_out


# =============================================================================
# PROTEIN PREPARATION
# =============================================================================

def prepare_protein(pdb_file, remove_waters, output_dir='prep'):
    """Prepare protein with pdb4amber, handling KCX residues and Zn ions"""
    print(f"\n{'='*60}")
    print(f"PROTEIN PREPARATION")
    print(f"{'='*60}")
    print(f"  PDB file: {pdb_file}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Read PDB and check for special residues
    with open(pdb_file, 'r') as f:
        pdb_content = f.read()
    
    has_kcx = 'KCX' in pdb_content
    has_zn = ' ZN ' in pdb_content or 'ZN2' in pdb_content
    
    print(f"  Contains KCX (carboxylated lysine): {has_kcx}")
    print(f"  Contains Zn2+ ions: {has_zn}")
    
    protein_out = os.path.join(output_dir, 'protein.pdb')
    
    # Run pdb4amber
    cmd = ['pdb4amber', '-i', pdb_file, '-o', protein_out]
    if remove_waters == 'yes':
        cmd.append('--dry')
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  WARNING: pdb4amber warning: {result.stderr}")
        # Copy original file if pdb4amber fails
        shutil.copy(pdb_file, protein_out)
    
    print(f"  Protein prepared: {protein_out}")
    return protein_out, has_kcx, has_zn


# =============================================================================
# TLEAP SYSTEM BUILDING
# =============================================================================

def create_tleap_input(settings, protein_pdb, ligand_mol2, ligand_frcmod, 
                       kcx_frcmod, kcx_lib, has_kcx, has_zn):
    """Create tleap input script for system building"""
    
    # Select force field files
    if settings.forceField == 'ff19SB':
        ff_protein = 'leaprc.protein.ff19SB'
    else:
        ff_protein = 'leaprc.protein.ff14SB'
    
    # Select water model
    water_models = {
        'tip3p': 'leaprc.water.tip3p',
        'tip3pfb': 'leaprc.water.tip3p',
        'tip4pew': 'leaprc.water.tip4pew',
        'tip4pfb': 'leaprc.water.tip4pew',
        'spce': 'leaprc.water.spce',
        'opc': 'leaprc.water.opc',
        'opc3': 'leaprc.water.opc3'
    }
    ff_water = water_models.get(settings.waterModel, 'leaprc.water.tip3p')
    
    # Map ion names for tleap
    ion_map = {
        'Na+': 'Na+', 'K+': 'K+', 'Li+': 'Li+', 'Cs+': 'Cs+', 'Rb+': 'Rb+',
        'Cl-': 'Cl-', 'Br-': 'Br-', 'F-': 'F-', 'I-': 'I-'
    }
    pos_ion = ion_map.get(settings.positiveIon, 'Na+')
    neg_ion = ion_map.get(settings.negativeIon, 'Cl-')
    
    # Calculate number of ions for ionic strength
    n_ions = max(1, int(settings.ionicStrength * 60))
    
    # Box size in Angstroms
    box_size = settings.boxSize
    
    script = f"""# tleap input for OpenMM-KCX system
# Force fields
source {ff_protein}
source {ff_water}
source leaprc.gaff2
"""
    
    # Load KCX parameters if needed
    if has_kcx:
        script += f"""
# Load KCX (carboxylated lysine) parameters
loadamberparams {kcx_frcmod}
loadoff {kcx_lib}
"""
    
    # Load Zn parameters if needed
    if has_zn:
        script += """
# Load Zn2+ parameters
loadamberparams frcmod.ions234lm_126_tip3p
"""
    
    # Load ligand if present
    if ligand_mol2 and os.path.exists(ligand_mol2):
        script += f"""
# Load ligand
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}
"""
    
    # Load protein
    script += f"""
# Load protein
mol = loadpdb {protein_pdb}
check mol
"""
    
    # Combine protein and ligand
    if ligand_mol2 and os.path.exists(ligand_mol2):
        script += """
# Combine protein and ligand
complex = combine {mol LIG}
"""
    else:
        script += """
# Use protein only
complex = mol
"""
    
    # Solvate based on box shape
    if settings.boxShape == 'octahedron':
        script += f"""
# Solvate with truncated octahedron
solvateOct complex TIP3PBOX {box_size}
"""
    else:
        script += f"""
# Solvate with cubic box
solvatebox complex TIP3PBOX {box_size}
"""
    
    # Add ions
    script += f"""
# Neutralize and add ions
addions complex {pos_ion} 0
addions complex {neg_ion} 0
addions complex {pos_ion} {n_ions}
addions complex {neg_ion} {n_ions}

# Check and save
check complex

# Save AMBER files
saveamberparm complex system.prmtop system.inpcrd

# Save PDB for reference
savepdb complex system.pdb

quit
"""
    
    # Write tleap input file
    tleap_file = 'tleap.in'
    with open(tleap_file, 'w') as f:
        f.write(script)
    
    return tleap_file


def run_tleap(tleap_file):
    """Run tleap to build the system"""
    print(f"\n{'='*60}")
    print(f"RUNNING TLEAP")
    print(f"{'='*60}")
    
    result = subprocess.run(['tleap', '-f', tleap_file], capture_output=True, text=True)
    
    print(result.stdout)
    if result.stderr:
        print("STDERR:", result.stderr)
    
    # Check if output files were created
    if not os.path.exists('system.prmtop') or not os.path.exists('system.inpcrd'):
        print("ERROR: tleap failed to create topology files!")
        print("Check tleap.log for details")
        sys.exit(1)
    
    print("  System built successfully!")
    return True


# =============================================================================
# OPENMM SIMULATION
# =============================================================================

def run_simulation(settings):
    """Run OpenMM MD simulation with all specified parameters"""
    print(f"\n{'='*60}")
    print(f"OPENMM SIMULATION")
    print(f"{'='*60}")
    
    # Import OpenMM
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    
    # Create output directory
    os.makedirs('out', exist_ok=True)
    
    # Load AMBER files
    print("Loading AMBER topology and coordinates...")
    prmtop = app.AmberPrmtopFile('system.prmtop')
    inpcrd = app.AmberInpcrdFile('system.inpcrd')
    
    # Create System
    print("Creating OpenMM system...")
    
    # Nonbonded method
    nb_methods = {
        'PME': app.PME,
        'NoCutoff': app.NoCutoff,
        'CutoffPeriodic': app.CutoffPeriodic,
        'LJPME': app.LJPME
    }
    nonbonded_method = nb_methods.get(settings.nonbondedMethod, app.PME)
    
    # Constraints
    constraint_types = {
        'None': None,
        'HBonds': app.HBonds,
        'AllBonds': app.AllBonds,
        'HAngles': app.HAngles
    }
    constraints = constraint_types.get(settings.constraints, app.HBonds)
    
    # Create system with specified parameters
    system_kwargs = {
        'nonbondedMethod': nonbonded_method,
        'nonbondedCutoff': settings.nonbondedCutoff * unit.nanometer,
        'constraints': constraints,
        'rigidWater': settings.rigidWater == 'yes'
    }
    
    # Add hydrogen mass repartitioning if specified
    if settings.hydrogenMass is not None:
        system_kwargs['hydrogenMass'] = settings.hydrogenMass * unit.amu
    
    system = prmtop.createSystem(**system_kwargs)
    
    # Set PME parameters if using PME
    if settings.nonbondedMethod in ['PME', 'LJPME']:
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setEwaldErrorTolerance(settings.ewaldErrorTolerance)
                if settings.switchDistance is not None:
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(settings.switchDistance * unit.nanometer)
    
    # Add barostat for NPT
    system.addForce(mm.MonteCarloBarostat(
        settings.pressure * unit.bar,
        settings.temperature * unit.kelvin,
        settings.barostatInterval
    ))
    
    # Create Integrator
    print(f"Creating {settings.integrator} integrator...")
    
    timestep = settings.timestep * unit.femtosecond
    temperature = settings.temperature * unit.kelvin
    friction = settings.frictionCoeff / unit.picosecond
    
    if settings.integrator == 'LangevinMiddle':
        integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
    elif settings.integrator == 'Langevin':
        integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    elif settings.integrator == 'NoseHoover':
        integrator = mm.NoseHooverIntegrator(temperature, friction, timestep)
    elif settings.integrator == 'Verlet':
        integrator = mm.VerletIntegrator(timestep)
        system.addForce(mm.AndersenThermostat(temperature, friction))
    elif settings.integrator == 'Brownian':
        integrator = mm.BrownianIntegrator(temperature, friction, timestep)
    else:
        integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # Create Simulation
    print("Setting up simulation...")
    
    # Try CUDA, fall back to CPU
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        print("  Using CUDA platform")
    except Exception:
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
            properties = {}
            print("  Using OpenCL platform")
        except Exception:
            platform = mm.Platform.getPlatformByName('CPU')
            properties = {}
            print("  Using CPU platform")
    
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    # Energy Minimization
    print(f"\nMinimizing energy ({settings.minimizationSteps} steps max)...")
    
    initial_state = simulation.context.getState(getEnergy=True)
    print(f"  Initial energy: {initial_state.getPotentialEnergy()}")
    
    simulation.minimizeEnergy(
        tolerance=settings.minimizationTolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=settings.minimizationSteps
    )
    
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"  Final energy: {final_state.getPotentialEnergy()}")
    
    # Save minimized structure
    with open('out/minimized.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)
    
    # Equilibration
    equil_steps = int(settings.equilibrationTime * 1e6 / settings.timestep)
    print(f"\nEquilibrating for {settings.equilibrationTime} ns ({equil_steps} steps)...")
    
    # Add reporters for equilibration
    simulation.reporters.append(app.StateDataReporter(
        'out/equilibration.log', settings.equilTrajFreq,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=equil_steps
    ))
    
    if settings.trajectoryFormat == 'DCD':
        simulation.reporters.append(app.DCDReporter('out/equilibration.dcd', settings.equilTrajFreq))
    elif settings.trajectoryFormat == 'XTC':
        simulation.reporters.append(app.XTCReporter('out/equilibration.xtc', settings.equilTrajFreq))
    
    simulation.step(equil_steps)
    
    # Clear reporters
    simulation.reporters.clear()
    
    # Production
    prod_steps = int(settings.productionTime * 1e6 / settings.timestep)
    print(f"\nRunning production for {settings.productionTime} ns ({prod_steps} steps)...")
    
    # Add reporters for production
    simulation.reporters.append(app.StateDataReporter(
        'out/production.log', settings.prodTrajFreq,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=prod_steps
    ))
    
    if settings.trajectoryFormat == 'DCD':
        simulation.reporters.append(app.DCDReporter('out/production.dcd', settings.prodTrajFreq))
    elif settings.trajectoryFormat == 'XTC':
        simulation.reporters.append(app.XTCReporter('out/production.xtc', settings.prodTrajFreq))
    else:
        simulation.reporters.append(app.PDBReporter('out/production.pdb', settings.prodTrajFreq))
    
    simulation.reporters.append(app.CheckpointReporter('out/checkpoint.chk', settings.checkpointFreq))
    
    simulation.step(prod_steps)
    
    # Save final structure
    state = simulation.context.getState(getPositions=True, getVelocities=True)
    with open('out/final.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    simulation.saveState('out/final_state.xml')
    
    print("\nSimulation complete!")


# =============================================================================
# TRAJECTORY ANALYSIS
# =============================================================================

def run_analysis(settings):
    """Run trajectory analysis with MDTraj"""
    print(f"\n{'='*60}")
    print(f"TRAJECTORY ANALYSIS")
    print(f"{'='*60}")
    
    try:
        import mdtraj as md
        import numpy as np
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except ImportError as e:
        print(f"  WARNING: Analysis skipped - missing package: {e}")
        return
    
    # Determine trajectory file
    if settings.trajectoryFormat == 'DCD':
        traj_file = 'out/production.dcd'
    elif settings.trajectoryFormat == 'XTC':
        traj_file = 'out/production.xtc'
    else:
        traj_file = 'out/production.pdb'
    
    if not os.path.exists(traj_file):
        print(f"  Trajectory file not found: {traj_file}")
        return
    
    print(f"  Loading trajectory: {traj_file}")
    traj = md.load(traj_file, top='system.pdb', stride=settings.stepSize)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    # RMSD
    print("  Calculating RMSD...")
    
    if settings.rmsdMask == 'backbone':
        rmsd_atoms = traj.topology.select('backbone')
    elif settings.rmsdMask == 'heavy':
        rmsd_atoms = traj.topology.select('mass > 2')
    elif settings.rmsdMask == 'protein':
        rmsd_atoms = traj.topology.select('protein')
    else:
        rmsd_atoms = traj.topology.select('all')
    
    rmsd = md.rmsd(traj, traj, 0, atom_indices=rmsd_atoms) * 10  # nm to Å
    time_ns = traj.time / 1000  # ps to ns
    
    np.savetxt('out/rmsd.csv', 
               np.column_stack([time_ns, rmsd]),
               delimiter=',', header='time_ns,rmsd_angstrom', comments='')
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, rmsd, 'b-', linewidth=0.5)
    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('RMSD (Å)', fontsize=12)
    plt.title(f'Backbone RMSD (mean: {np.mean(rmsd):.2f} Å)', fontsize=14)
    plt.tight_layout()
    plt.savefig('out/rmsd.png', dpi=150)
    plt.close()
    
    # RMSF
    print("  Calculating RMSF...")
    
    ca_atoms = traj.topology.select('name CA')
    if len(ca_atoms) > 0:
        rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_atoms) * 10  # nm to Å
        residues = [traj.topology.atom(i).residue.resSeq for i in ca_atoms]
        
        np.savetxt('out/rmsf.csv',
                   np.column_stack([residues, rmsf]),
                   delimiter=',', header='residue,rmsf_angstrom', comments='')
        
        plt.figure(figsize=(12, 6))
        plt.plot(residues, rmsf, 'b-', linewidth=1)
        plt.xlabel('Residue Number', fontsize=12)
        plt.ylabel('RMSF (Å)', fontsize=12)
        plt.title(f'Per-Residue RMSF (mean: {np.mean(rmsf):.2f} Å)', fontsize=14)
        plt.tight_layout()
        plt.savefig('out/rmsf.png', dpi=150)
        plt.close()
    
    print("  Analysis complete!")


# =============================================================================
# MAIN
# =============================================================================

def main():
    # Get settings from environment variables (Tamarind convention)
    settings_dict = get_settings()
    settings = Settings(settings_dict)
    
    print("="*70)
    print("  OpenMM MD Simulation with KCX Support - v6 (Tamarind)")
    print("="*70)
    print(f"\nJob Name: {settings.jobName}")
    print("\nSIMULATION PARAMETERS:")
    print("-"*70)
    print(f"  Input PDB:           {settings.pdbFile}")
    print(f"  Ligand file:         {settings.ligandFile}")
    print(f"  Ligand charge:       {settings.ligandCharge}")
    print("-"*70)
    print(f"  Force field:         {settings.forceField}")
    print(f"  Water model:         {settings.waterModel}")
    print(f"  Box size:            {settings.boxSize} Å")
    print(f"  Box shape:           {settings.boxShape}")
    print(f"  Ionic strength:      {settings.ionicStrength} M")
    print(f"  pH:                  {settings.pH}")
    print("-"*70)
    print(f"  Minimization steps:  {settings.minimizationSteps}")
    print(f"  Equilibration:       {settings.equilibrationTime} ns")
    print(f"  Production:          {settings.productionTime} ns")
    print(f"  Timestep:            {settings.timestep} fs")
    print(f"  Temperature:         {settings.temperature} K")
    print(f"  Pressure:            {settings.pressure} bar")
    print("-"*70)
    print(f"  Integrator:          {settings.integrator}")
    print(f"  Constraints:         {settings.constraints}")
    print(f"  Nonbonded method:    {settings.nonbondedMethod}")
    print(f"  Cutoff:              {settings.nonbondedCutoff} nm")
    print("="*70)
    
    # Check input file exists
    if not os.path.exists(settings.pdbFile):
        print(f"ERROR: PDB file not found: {settings.pdbFile}")
        print("Available files in inputs/:")
        for f in glob.glob('inputs/*'):
            print(f"  {f}")
        sys.exit(1)
    
    # Write KCX parameters
    kcx_frcmod, kcx_lib = write_kcx_parameters('params')
    
    # Prepare ligand if provided
    ligand_mol2, ligand_frcmod = None, None
    if settings.ligandFile and os.path.exists(settings.ligandFile):
        ligand_mol2, ligand_frcmod = prepare_ligand(
            settings.ligandFile, settings.ligandCharge, 'prep'
        )
    
    # Prepare protein
    protein_pdb, has_kcx, has_zn = prepare_protein(
        settings.pdbFile, settings.removeWaters, 'prep'
    )
    
    # Create tleap input
    tleap_file = create_tleap_input(
        settings, protein_pdb, ligand_mol2, ligand_frcmod,
        kcx_frcmod, kcx_lib, has_kcx, has_zn
    )
    
    # Run tleap
    run_tleap(tleap_file)
    
    # Run simulation
    run_simulation(settings)
    
    # Run analysis
    run_analysis(settings)
    
    print("\n" + "="*70)
    print("  ALL DONE!")
    print("="*70)
    print("\nOutput files in out/:")
    print("  - minimized.pdb      : Minimized structure")
    print("  - equilibration.*    : Equilibration trajectory")
    print("  - production.*       : Production trajectory")
    print("  - final.pdb          : Final structure")
    print("  - rmsd.csv/png       : RMSD analysis")
    print("  - rmsf.csv/png       : RMSF analysis")
    print("="*70)


if __name__ == '__main__':
    main()
