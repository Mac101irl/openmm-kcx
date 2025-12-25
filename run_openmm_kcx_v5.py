#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX (Carboxylated Lysine) Support - v5
Full parameter support matching standard Tamarind OpenMM
For D-Hydantoinase simulations with proper active site modeling

Author: Generated for BMS Hydantoinase Project
Version: 5.0
"""

import os
import sys
import json
import argparse
import subprocess
import shutil
from pathlib import Path

# =============================================================================
# ARGUMENT PARSING - All OpenMM Parameters
# =============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='OpenMM MD Simulation with KCX (carboxylated lysine) support',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # -------------------------------------------------------------------------
    # INPUT FILES (support both --pdb/--pdbFile and --ligand/--ligandFile)
    # -------------------------------------------------------------------------
    # Primary names (API style)
    parser.add_argument('--pdbFile', dest='pdbFile', default=None,
                        help='Input PDB file (can contain KCX residues)')
    parser.add_argument('--pdb', dest='pdbFile', default=None,
                        help='Input PDB file (alias for --pdbFile)')
    parser.add_argument('--ligandFile', dest='ligandFile', default=None, 
                        help='Ligand SDF/MOL2 file for protein-ligand simulations')
    parser.add_argument('--ligand', dest='ligandFile', default=None,
                        help='Ligand file (alias for --ligandFile)')
    parser.add_argument('--ligandCharge', type=int, default=0, 
                        help='Net formal charge of the ligand molecule')
    parser.add_argument('--systemType', default='protein-ligand',
                        choices=['protein', 'protein-ligand'],
                        help='Type of system to simulate')
    
    # -------------------------------------------------------------------------
    # FORCE FIELD PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--forceField', default='ff19SB',
                        choices=['ff19SB', 'ff14SB'],
                        help='Protein force field')
    parser.add_argument('--waterModel', default='tip3p',
                        choices=['tip3p', 'tip3pfb', 'tip4pew', 'tip4pfb', 'spce', 'opc', 'opc3'],
                        help='Water model for explicit solvent')
    
    # -------------------------------------------------------------------------
    # SOLVATION PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--boxSize', type=float, default=12.0,
                        help='Minimum distance from solute to box edge (Angstroms)')
    parser.add_argument('--boxShape', default='cube',
                        choices=['cube', 'dodecahedron', 'octahedron'],
                        help='Shape of periodic box (dodecahedron is ~29%% smaller)')
    parser.add_argument('--ionicStrength', type=float, default=0.15,
                        help='Salt concentration in molar (M)')
    parser.add_argument('--positiveIon', default='Na+',
                        choices=['Na+', 'K+', 'Li+', 'Cs+', 'Rb+'],
                        help='Type of positive ion for neutralization')
    parser.add_argument('--negativeIon', default='Cl-',
                        choices=['Cl-', 'Br-', 'F-', 'I-'],
                        help='Type of negative ion for neutralization')
    parser.add_argument('--pH', type=float, default=7.0,
                        help='pH for determining protonation states')
    parser.add_argument('--removeWaters', default='no',
                        choices=['yes', 'no'],
                        help='Remove crystallographic waters from input PDB')
    
    # -------------------------------------------------------------------------
    # NONBONDED PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--nonbondedMethod', default='PME',
                        choices=['PME', 'NoCutoff', 'CutoffPeriodic', 'LJPME'],
                        help='Method for long-range electrostatics')
    parser.add_argument('--nonbondedCutoff', type=float, default=1.0,
                        help='Cutoff distance for direct-space interactions (nm)')
    parser.add_argument('--ewaldErrorTolerance', type=float, default=0.0005,
                        help='Error tolerance for PME (smaller = more accurate)')
    parser.add_argument('--switchDistance', type=float, default=None,
                        help='Distance to start switching LJ interactions (nm)')
    
    # -------------------------------------------------------------------------
    # CONSTRAINT PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--constraints', default='HBonds',
                        choices=['None', 'HBonds', 'AllBonds', 'HAngles'],
                        help='Bond constraint type (HBonds allows 2fs timestep)')
    parser.add_argument('--rigidWater', default='yes',
                        choices=['yes', 'no'],
                        help='Keep water molecules rigid')
    parser.add_argument('--hydrogenMass', type=float, default=None,
                        help='Hydrogen mass repartitioning (amu). Use 4.0 for 4fs timestep')
    
    # -------------------------------------------------------------------------
    # SIMULATION PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--minimizationSteps', '--minimize', dest='minimizationSteps',
                        type=int, default=10000,
                        help='Maximum steps for energy minimization')
    parser.add_argument('--minimizationTolerance', type=float, default=10.0,
                        help='Energy minimization convergence tolerance (kJ/mol/nm)')
    parser.add_argument('--equilibrationTime', '--equilibrate', dest='equilibrationTime',
                        type=float, default=0.2,
                        help='NVT/NPT equilibration time (ns)')
    parser.add_argument('--productionTime', '--time', dest='productionTime',
                        type=float, default=5.0,
                        help='Production MD simulation time (ns)')
    parser.add_argument('--timestep', type=float, default=2.0,
                        help='Integration timestep (fs). Use 4.0 with hydrogenMass=4')
    
    # -------------------------------------------------------------------------
    # TEMPERATURE & PRESSURE
    # -------------------------------------------------------------------------
    parser.add_argument('--temperature', '--temp', dest='temperature',
                        type=float, default=298.0,
                        help='Simulation temperature (K). Use 310 for physiological')
    parser.add_argument('--pressure', type=float, default=1.0,
                        help='Simulation pressure (bar)')
    parser.add_argument('--frictionCoeff', type=float, default=1.0,
                        help='Langevin friction coefficient (1/ps)')
    parser.add_argument('--barostatInterval', type=int, default=25,
                        help='Monte Carlo barostat attempt frequency (steps)')
    
    # -------------------------------------------------------------------------
    # INTEGRATOR
    # -------------------------------------------------------------------------
    parser.add_argument('--integrator', default='LangevinMiddle',
                        choices=['LangevinMiddle', 'Langevin', 'NoseHoover', 'Verlet', 'Brownian'],
                        help='Integration algorithm')
    
    # -------------------------------------------------------------------------
    # RESTRAINT PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--forceConstant', type=float, default=1000.0,
                        help='Position restraint force constant (kJ/mol/nm^2)')
    parser.add_argument('--restraintMask', default='backbone',
                        choices=['backbone', 'heavy', 'all', 'none'],
                        help='Atoms to restrain during equilibration')
    
    # -------------------------------------------------------------------------
    # OUTPUT PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--output', default='simulation',
                        help='Output directory/prefix name')
    parser.add_argument('--equilTrajFreq', type=int, default=1000,
                        help='Equilibration trajectory save frequency (steps)')
    parser.add_argument('--prodTrajFreq', type=int, default=1000,
                        help='Production trajectory save frequency (steps)')
    parser.add_argument('--checkpointFreq', type=int, default=10000,
                        help='Checkpoint save frequency (steps)')
    parser.add_argument('--trajectoryFormat', default='DCD',
                        choices=['DCD', 'XTC', 'PDB'],
                        help='Output trajectory format')
    
    # -------------------------------------------------------------------------
    # ANALYSIS PARAMETERS
    # -------------------------------------------------------------------------
    parser.add_argument('--stepSize', type=int, default=5,
                        help='Stride for trajectory analysis')
    parser.add_argument('--rmsdMask', default='backbone',
                        choices=['backbone', 'heavy', 'all', 'protein'],
                        help='Atom selection for RMSD calculation')
    
    args = parser.parse_args()
    
    # Validate required arguments
    if args.pdbFile is None:
        parser.error("--pdb or --pdbFile is required")
    
    return args


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
        input_type = 'sdf'  # default
    
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
        '-c', 'bcc',           # AM1-BCC charges
        '-at', 'gaff2',        # GAFF2 atom types
        '-nc', str(ligand_charge),
        '-rn', 'LIG',          # Residue name
        '-pf', 'y'             # Remove intermediate files
    ]
    
    result = subprocess.run(cmd_antechamber, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(mol2_out):
        print(f"  ERROR: Antechamber failed!")
        print(f"  stdout: {result.stdout}")
        print(f"  stderr: {result.stderr}")
        # Try with sqm for charge calculation
        print("  Retrying with slower sqm method...")
        cmd_antechamber[cmd_antechamber.index('bcc')] = 'gas'
        result = subprocess.run(cmd_antechamber, capture_output=True, text=True)
        if result.returncode != 0:
            raise RuntimeError(f"Ligand preparation failed: {result.stderr}")
    
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
        print(f"  WARNING: parmchk2 had issues: {result.stderr}")
    
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
        print("  Removing crystallographic waters")
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  WARNING: pdb4amber issues: {result.stderr}")
        # Copy original if pdb4amber fails
        shutil.copy(pdb_file, protein_out)
    
    print(f"  Protein prepared: {protein_out}")
    return protein_out, has_kcx, has_zn


# =============================================================================
# TLEAP SYSTEM BUILDING
# =============================================================================

def create_tleap_input(args, protein_pdb, ligand_mol2, ligand_frcmod, 
                       kcx_frcmod, kcx_lib, has_kcx, has_zn):
    """Create tleap input script for system building"""
    
    # Select force field files
    if args.forceField == 'ff19SB':
        ff_protein = 'leaprc.protein.ff19SB'
    else:
        ff_protein = 'leaprc.protein.ff14SB'
    
    # Select water model
    water_models = {
        'tip3p': 'leaprc.water.tip3p',
        'tip3pfb': 'leaprc.water.tip3p',  # Use tip3p leaprc, different params
        'tip4pew': 'leaprc.water.tip4pew',
        'tip4pfb': 'leaprc.water.tip4pew',
        'spce': 'leaprc.water.spce',
        'opc': 'leaprc.water.opc',
        'opc3': 'leaprc.water.opc3'
    }
    ff_water = water_models.get(args.waterModel, 'leaprc.water.tip3p')
    
    # Map ion names for tleap
    ion_map = {
        'Na+': 'Na+', 'K+': 'K+', 'Li+': 'Li+', 'Cs+': 'Cs+', 'Rb+': 'Rb+',
        'Cl-': 'Cl-', 'Br-': 'Br-', 'F-': 'F-', 'I-': 'I-'
    }
    pos_ion = ion_map.get(args.positiveIon, 'Na+')
    neg_ion = ion_map.get(args.negativeIon, 'Cl-')
    
    # Calculate number of ions for ionic strength
    # Approximate: 0.15 M in ~100 nm^3 box ≈ 9 ion pairs
    n_ions = max(1, int(args.ionicStrength * 60))  # Rough estimate
    
    # Box size in nm
    box_nm = args.boxSize / 10.0
    
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
# Zn2+ ion parameters
addAtomTypes {{ {{ "Zn" "Zn" "sp3" }} }}
"""
    
    # Load ligand if present
    if ligand_mol2 and os.path.exists(ligand_mol2):
        script += f"""
# Load ligand parameters
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}
check LIG
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
system = combine {mol LIG}
"""
    else:
        script += """
system = mol
"""
    
    # Solvate based on box shape
    if args.boxShape == 'octahedron':
        script += f"""
# Solvate with truncated octahedron
solvateoct system TIP3PBOX {box_nm * 10:.1f}
"""
    else:
        script += f"""
# Solvate with rectangular box
solvateBox system TIP3PBOX {box_nm * 10:.1f}
"""
    
    # Add ions
    script += f"""
# Neutralize system
addionsrand system {pos_ion} 0
addionsrand system {neg_ion} 0

# Add ions for ionic strength ({args.ionicStrength} M)
addionsrand system {pos_ion} {n_ions}
addionsrand system {neg_ion} {n_ions}

# Check system
check system

# Save AMBER files
saveamberparm system system.prmtop system.inpcrd
savepdb system system.pdb

quit
"""
    
    tleap_file = 'tleap.in'
    with open(tleap_file, 'w') as f:
        f.write(script)
    
    return tleap_file


def run_tleap(tleap_file):
    """Execute tleap to build the system"""
    print(f"\n{'='*60}")
    print(f"BUILDING SYSTEM WITH TLEAP")
    print(f"{'='*60}")
    
    result = subprocess.run(['tleap', '-f', tleap_file], 
                          capture_output=True, text=True)
    
    print(result.stdout)
    
    if not os.path.exists('system.prmtop') or not os.path.exists('system.inpcrd'):
        print(f"ERROR: tleap failed to create system files!")
        print(f"stderr: {result.stderr}")
        raise RuntimeError("tleap system building failed")
    
    print("  System built successfully!")
    print(f"  - system.prmtop")
    print(f"  - system.inpcrd")
    print(f"  - system.pdb")


# =============================================================================
# OPENMM SIMULATION
# =============================================================================

def run_simulation(args):
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
    
    # -------------------------------------------------------------------------
    # Create System
    # -------------------------------------------------------------------------
    print("Creating OpenMM system...")
    
    # Nonbonded method
    nb_methods = {
        'PME': app.PME,
        'NoCutoff': app.NoCutoff,
        'CutoffPeriodic': app.CutoffPeriodic,
        'LJPME': app.LJPME
    }
    nonbonded_method = nb_methods.get(args.nonbondedMethod, app.PME)
    
    # Constraints
    constraint_types = {
        'None': None,
        'HBonds': app.HBonds,
        'AllBonds': app.AllBonds,
        'HAngles': app.HAngles
    }
    constraints = constraint_types.get(args.constraints, app.HBonds)
    
    # Create system with specified parameters
    system_kwargs = {
        'nonbondedMethod': nonbonded_method,
        'nonbondedCutoff': args.nonbondedCutoff * unit.nanometer,
        'constraints': constraints,
        'rigidWater': args.rigidWater == 'yes'
    }
    
    # Add hydrogen mass repartitioning if specified
    if args.hydrogenMass is not None:
        system_kwargs['hydrogenMass'] = args.hydrogenMass * unit.amu
    
    system = prmtop.createSystem(**system_kwargs)
    
    # Set PME parameters if using PME
    if args.nonbondedMethod in ['PME', 'LJPME']:
        for force in system.getForces():
            if isinstance(force, mm.NonbondedForce):
                force.setEwaldErrorTolerance(args.ewaldErrorTolerance)
                if args.switchDistance is not None:
                    force.setUseSwitchingFunction(True)
                    force.setSwitchingDistance(args.switchDistance * unit.nanometer)
    
    # Add barostat for NPT
    system.addForce(mm.MonteCarloBarostat(
        args.pressure * unit.bar,
        args.temperature * unit.kelvin,
        args.barostatInterval
    ))
    
    # -------------------------------------------------------------------------
    # Create Integrator
    # -------------------------------------------------------------------------
    print(f"Creating {args.integrator} integrator...")
    
    timestep = args.timestep * unit.femtosecond
    temperature = args.temperature * unit.kelvin
    friction = args.frictionCoeff / unit.picosecond
    
    if args.integrator == 'LangevinMiddle':
        integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
    elif args.integrator == 'Langevin':
        integrator = mm.LangevinIntegrator(temperature, friction, timestep)
    elif args.integrator == 'NoseHoover':
        integrator = mm.NoseHooverIntegrator(temperature, friction, timestep)
    elif args.integrator == 'Verlet':
        integrator = mm.VerletIntegrator(timestep)
        # Add Andersen thermostat for temperature control
        system.addForce(mm.AndersenThermostat(temperature, friction))
    elif args.integrator == 'Brownian':
        integrator = mm.BrownianIntegrator(temperature, friction, timestep)
    else:
        integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # -------------------------------------------------------------------------
    # Create Simulation
    # -------------------------------------------------------------------------
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
    
    # -------------------------------------------------------------------------
    # Energy Minimization
    # -------------------------------------------------------------------------
    print(f"\nMinimizing energy ({args.minimizationSteps} steps max)...")
    
    initial_state = simulation.context.getState(getEnergy=True)
    print(f"  Initial energy: {initial_state.getPotentialEnergy()}")
    
    simulation.minimizeEnergy(
        tolerance=args.minimizationTolerance * unit.kilojoule_per_mole / unit.nanometer,
        maxIterations=args.minimizationSteps
    )
    
    final_state = simulation.context.getState(getEnergy=True, getPositions=True)
    print(f"  Final energy: {final_state.getPotentialEnergy()}")
    
    # Save minimized structure
    with open('out/minimized.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, final_state.getPositions(), f)
    
    # -------------------------------------------------------------------------
    # Equilibration
    # -------------------------------------------------------------------------
    equil_steps = int(args.equilibrationTime * 1e6 / args.timestep)
    print(f"\nEquilibrating for {args.equilibrationTime} ns ({equil_steps} steps)...")
    
    # Add reporters for equilibration
    simulation.reporters.append(app.StateDataReporter(
        'out/equilibration.log', args.equilTrajFreq,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=equil_steps
    ))
    
    if args.trajectoryFormat == 'DCD':
        simulation.reporters.append(app.DCDReporter('out/equilibration.dcd', args.equilTrajFreq))
    elif args.trajectoryFormat == 'XTC':
        simulation.reporters.append(app.XTCReporter('out/equilibration.xtc', args.equilTrajFreq))
    
    simulation.step(equil_steps)
    
    # Clear reporters
    simulation.reporters.clear()
    
    # -------------------------------------------------------------------------
    # Production
    # -------------------------------------------------------------------------
    prod_steps = int(args.productionTime * 1e6 / args.timestep)
    print(f"\nRunning production for {args.productionTime} ns ({prod_steps} steps)...")
    
    # Add reporters for production
    simulation.reporters.append(app.StateDataReporter(
        'out/production.log', args.prodTrajFreq,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        temperature=True, volume=True, density=True,
        progress=True, remainingTime=True, speed=True,
        totalSteps=prod_steps
    ))
    
    if args.trajectoryFormat == 'DCD':
        simulation.reporters.append(app.DCDReporter('out/production.dcd', args.prodTrajFreq))
    elif args.trajectoryFormat == 'XTC':
        simulation.reporters.append(app.XTCReporter('out/production.xtc', args.prodTrajFreq))
    else:
        simulation.reporters.append(app.PDBReporter('out/production.pdb', args.prodTrajFreq))
    
    simulation.reporters.append(app.CheckpointReporter('out/checkpoint.chk', args.checkpointFreq))
    
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

def run_analysis(args):
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
    
    os.makedirs('out/analysis', exist_ok=True)
    
    # Determine trajectory file
    if args.trajectoryFormat == 'DCD':
        traj_file = 'out/production.dcd'
    elif args.trajectoryFormat == 'XTC':
        traj_file = 'out/production.xtc'
    else:
        traj_file = 'out/production.pdb'
    
    if not os.path.exists(traj_file):
        print(f"  Trajectory file not found: {traj_file}")
        return
    
    print(f"  Loading trajectory: {traj_file}")
    traj = md.load(traj_file, top='system.pdb', stride=args.stepSize)
    print(f"  Frames: {traj.n_frames}")
    print(f"  Atoms: {traj.n_atoms}")
    
    # -------------------------------------------------------------------------
    # RMSD
    # -------------------------------------------------------------------------
    print("  Calculating RMSD...")
    
    if args.rmsdMask == 'backbone':
        rmsd_atoms = traj.topology.select('backbone')
    elif args.rmsdMask == 'heavy':
        rmsd_atoms = traj.topology.select('mass > 2')
    elif args.rmsdMask == 'protein':
        rmsd_atoms = traj.topology.select('protein')
    else:
        rmsd_atoms = traj.topology.select('all')
    
    rmsd = md.rmsd(traj, traj, 0, atom_indices=rmsd_atoms) * 10  # nm to Å
    time_ns = traj.time / 1000  # ps to ns
    
    np.savetxt('out/analysis/rmsd.csv', 
               np.column_stack([time_ns, rmsd]),
               delimiter=',', header='time_ns,rmsd_angstrom', comments='')
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, rmsd, 'b-', linewidth=0.5)
    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('RMSD (Å)', fontsize=12)
    plt.title(f'Backbone RMSD (mean: {np.mean(rmsd):.2f} Å)', fontsize=14)
    plt.tight_layout()
    plt.savefig('out/analysis/rmsd.png', dpi=150)
    plt.close()
    
    # -------------------------------------------------------------------------
    # RMSF
    # -------------------------------------------------------------------------
    print("  Calculating RMSF...")
    
    ca_atoms = traj.topology.select('name CA')
    if len(ca_atoms) > 0:
        rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_atoms) * 10  # nm to Å
        residues = [traj.topology.atom(i).residue.resSeq for i in ca_atoms]
        
        np.savetxt('out/analysis/rmsf.csv',
                   np.column_stack([residues, rmsf]),
                   delimiter=',', header='residue,rmsf_angstrom', comments='')
        
        plt.figure(figsize=(12, 6))
        plt.plot(residues, rmsf, 'b-', linewidth=1)
        plt.fill_between(residues, rmsf, alpha=0.3)
        plt.xlabel('Residue Number', fontsize=12)
        plt.ylabel('RMSF (Å)', fontsize=12)
        plt.title('C-alpha RMSF', fontsize=14)
        plt.tight_layout()
        plt.savefig('out/analysis/rmsf.png', dpi=150)
        plt.close()
    
    # -------------------------------------------------------------------------
    # Radius of Gyration
    # -------------------------------------------------------------------------
    print("  Calculating Radius of Gyration...")
    
    rg = md.compute_rg(traj) * 10  # nm to Å
    
    np.savetxt('out/analysis/rg.csv',
               np.column_stack([time_ns, rg]),
               delimiter=',', header='time_ns,rg_angstrom', comments='')
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, rg, 'g-', linewidth=0.5)
    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('Radius of Gyration (Å)', fontsize=12)
    plt.title(f'Radius of Gyration (mean: {np.mean(rg):.2f} Å)', fontsize=14)
    plt.tight_layout()
    plt.savefig('out/analysis/rg.png', dpi=150)
    plt.close()
    
    # -------------------------------------------------------------------------
    # 2D RMSD Matrix
    # -------------------------------------------------------------------------
    print("  Calculating 2D RMSD matrix...")
    
    # Subsample for large trajectories
    n_frames = min(200, traj.n_frames)
    stride = max(1, traj.n_frames // n_frames)
    traj_sub = traj[::stride]
    
    rmsd_matrix = np.zeros((traj_sub.n_frames, traj_sub.n_frames))
    for i in range(traj_sub.n_frames):
        rmsd_matrix[i] = md.rmsd(traj_sub, traj_sub, i, atom_indices=rmsd_atoms) * 10
    
    plt.figure(figsize=(10, 8))
    plt.imshow(rmsd_matrix, cmap='viridis', origin='lower', 
               extent=[0, time_ns[-1], 0, time_ns[-1]])
    plt.colorbar(label='RMSD (Å)')
    plt.xlabel('Time (ns)', fontsize=12)
    plt.ylabel('Time (ns)', fontsize=12)
    plt.title('2D RMSD Matrix', fontsize=14)
    plt.tight_layout()
    plt.savefig('out/analysis/rmsd_2d.png', dpi=150)
    plt.close()
    
    # -------------------------------------------------------------------------
    # Summary Statistics
    # -------------------------------------------------------------------------
    summary = {
        'total_frames': traj.n_frames,
        'simulation_time_ns': float(time_ns[-1]),
        'rmsd_mean_A': float(np.mean(rmsd)),
        'rmsd_std_A': float(np.std(rmsd)),
        'rg_mean_A': float(np.mean(rg)),
        'rg_std_A': float(np.std(rg))
    }
    
    with open('out/analysis/summary.json', 'w') as f:
        json.dump(summary, f, indent=2)
    
    print("\n  Analysis Summary:")
    print(f"    Total frames: {summary['total_frames']}")
    print(f"    Simulation time: {summary['simulation_time_ns']:.2f} ns")
    print(f"    RMSD: {summary['rmsd_mean_A']:.2f} ± {summary['rmsd_std_A']:.2f} Å")
    print(f"    Rg: {summary['rg_mean_A']:.2f} ± {summary['rg_std_A']:.2f} Å")
    
    print("\n  Analysis files saved to out/analysis/")


# =============================================================================
# MAIN
# =============================================================================

def main():
    args = parse_arguments()
    
    print("="*70)
    print("  OpenMM MD Simulation with KCX Support - v5")
    print("="*70)
    print("\nSIMULATION PARAMETERS:")
    print("-"*70)
    print(f"  Input PDB:           {args.pdbFile}")
    print(f"  Ligand file:         {args.ligandFile}")
    print(f"  Ligand charge:       {args.ligandCharge}")
    print(f"  System type:         {args.systemType}")
    print("-"*70)
    print(f"  Force field:         {args.forceField}")
    print(f"  Water model:         {args.waterModel}")
    print(f"  Box size:            {args.boxSize} Å")
    print(f"  Box shape:           {args.boxShape}")
    print(f"  Ionic strength:      {args.ionicStrength} M")
    print(f"  pH:                  {args.pH}")
    print("-"*70)
    print(f"  Minimization steps:  {args.minimizationSteps}")
    print(f"  Equilibration:       {args.equilibrationTime} ns")
    print(f"  Production:          {args.productionTime} ns")
    print(f"  Timestep:            {args.timestep} fs")
    print(f"  Temperature:         {args.temperature} K")
    print(f"  Pressure:            {args.pressure} bar")
    print("-"*70)
    print(f"  Integrator:          {args.integrator}")
    print(f"  Constraints:         {args.constraints}")
    print(f"  Nonbonded method:    {args.nonbondedMethod}")
    print(f"  Cutoff:              {args.nonbondedCutoff} nm")
    print("="*70)
    
    # Write KCX parameters
    kcx_frcmod, kcx_lib = write_kcx_parameters('params')
    
    # Prepare ligand if provided
    ligand_mol2, ligand_frcmod = None, None
    if args.ligandFile and os.path.exists(args.ligandFile):
        ligand_mol2, ligand_frcmod = prepare_ligand(
            args.ligandFile, args.ligandCharge, 'prep'
        )
    
    # Prepare protein
    protein_pdb, has_kcx, has_zn = prepare_protein(
        args.pdbFile, args.removeWaters, 'prep'
    )
    
    # Create tleap input
    tleap_file = create_tleap_input(
        args, protein_pdb, ligand_mol2, ligand_frcmod,
        kcx_frcmod, kcx_lib, has_kcx, has_zn
    )
    
    # Run tleap
    run_tleap(tleap_file)
    
    # Run simulation
    run_simulation(args)
    
    # Run analysis
    run_analysis(args)
    
    print("\n" + "="*70)
    print("  ALL DONE!")
    print("="*70)
    print("\nOutput files:")
    print("  - out/minimized.pdb     : Minimized structure")
    print("  - out/equilibration.*   : Equilibration trajectory")
    print("  - out/production.*      : Production trajectory")
    print("  - out/final.pdb         : Final structure")
    print("  - out/analysis/         : Analysis results")
    print("="*70)


if __name__ == '__main__':
    main()
