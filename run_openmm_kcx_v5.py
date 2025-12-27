#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX (Carboxylated Lysine) Support
Tamarind Platform Compatible Version

Tamarind Conventions:
- File inputs: inputs/{fieldName} (no extension)
- Parameters: os.getenv('paramName', 'default')
- Outputs: Save ALL results to out/ directory
- Job name: os.getenv('JobName', 'job')

Author: Generated for BMS Hydantoinase Project
"""

import os
import sys
import subprocess
import shutil

# =============================================================================
# FILE INPUTS (Tamarind convention: inputs/{fieldName})
# =============================================================================
pdb_file = "inputs/pdbFile"
ligand_file = "inputs/ligandFile"

# =============================================================================
# PARAMETERS FROM ENVIRONMENT VARIABLES
# =============================================================================
job_name = os.getenv('JobName', 'openmm_simulation')
ligand_charge = int(os.getenv('ligandCharge', '0'))
force_field = os.getenv('forceField', 'ff19SB')
water_model = os.getenv('waterModel', 'tip3p')
box_size = float(os.getenv('boxSize', '12.0'))
ionic_strength = float(os.getenv('ionicStrength', '0.15'))
pH = float(os.getenv('pH', '7.0'))
minimization_steps = int(os.getenv('minimizationSteps', '10000'))
equilibration_time = float(os.getenv('equilibrationTime', '0.2'))
production_time = float(os.getenv('productionTime', '1.0'))
timestep = float(os.getenv('timestep', '2.0'))
temperature = float(os.getenv('temperature', '310.0'))
pressure = float(os.getenv('pressure', '1.0'))
constraints = os.getenv('constraints', 'HBonds')
integrator = os.getenv('integrator', 'LangevinMiddle')
prod_traj_freq = int(os.getenv('prodTrajFreq', '1000'))
remove_waters = os.getenv('removeWaters', 'no')

# Additional settings with defaults
nonbonded_method = os.getenv('nonbondedMethod', 'PME')
nonbonded_cutoff = float(os.getenv('nonbondedCutoff', '1.0'))
rigid_water = os.getenv('rigidWater', 'yes')
friction_coeff = float(os.getenv('frictionCoeff', '1.0'))
barostat_interval = int(os.getenv('barostatInterval', '25'))
equil_traj_freq = int(os.getenv('equilTrajFreq', '1000'))
checkpoint_freq = int(os.getenv('checkpointFreq', '10000'))
step_size = int(os.getenv('stepSize', '5'))

# =============================================================================
# KCX PARAMETER FILES
# =============================================================================

def create_kcx_frcmod():
    """Create AMBER frcmod file for KCX (carboxylated lysine) parameters"""
    return """KCX (N6-carboxylysine) parameters for AMBER ff19SB
MASS
c3 12.01         0.878
c  12.01         0.616
o  16.00         0.434
n  14.01         0.530

BOND
c3-n   330.6    1.456
c -n   427.6    1.379
c -o   637.7    1.218
n -hn  403.2    1.013

ANGLE
c3-c3-n    66.0      111.04
c3-n -c    63.9      121.35
c3-n -hn   46.0      116.78
c -n -hn   48.3      117.55
n -c -o    75.8      122.03
o -c -o    78.2      129.52

DIHE
c3-c3-n -c     1    0.650       180.0     2.
c3-c3-n -hn    1    0.000         0.0     2.
c3-n -c -o     1    2.500       180.0     2.
hn-n -c -o     1    2.500       180.0    -2.
hn-n -c -o     1    2.000         0.0     1.
X -c -n -X     4   10.000       180.0     2.

IMPROPER
c -n -o -o     1.1         180.         2.
n -c3-c -hn   1.1         180.         2.

NONBON
c3    1.9080  0.1094
c     1.9080  0.0860
o     1.6612  0.2100
n     1.8240  0.1700
hn    0.6000  0.0157
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

def prepare_ligand(lig_file, lig_charge, output_dir='prep'):
    """Prepare ligand with antechamber using GAFF2 and AM1-BCC charges"""
    print(f"\n{'='*60}")
    print(f"LIGAND PREPARATION")
    print(f"{'='*60}")
    print(f"  Ligand file: {lig_file}")
    print(f"  Net charge: {lig_charge}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # Detect input type from content or default to sdf
    input_type = 'sdf'
    
    mol2_out = os.path.join(output_dir, 'ligand.mol2')
    frcmod_out = os.path.join(output_dir, 'ligand.frcmod')
    
    print("  Running antechamber...")
    cmd = [
        'antechamber',
        '-i', lig_file,
        '-fi', input_type,
        '-o', mol2_out,
        '-fo', 'mol2',
        '-c', 'bcc',
        '-at', 'gaff2',
        '-nc', str(lig_charge),
        '-rn', 'LIG',
        '-pf', 'y'
    ]
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not os.path.exists(mol2_out):
        print(f"  WARNING: Antechamber BCC failed, trying gas charges...")
        cmd[cmd.index('bcc')] = 'gas'
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            print(f"  ERROR: Antechamber failed: {result.stderr}")
            return None, None
    
    print("  Running parmchk2...")
    cmd_parmchk = [
        'parmchk2',
        '-i', mol2_out,
        '-f', 'mol2',
        '-o', frcmod_out,
        '-s', 'gaff2'
    ]
    subprocess.run(cmd_parmchk, capture_output=True, text=True)
    
    print(f"  Ligand prepared: {mol2_out}")
    return mol2_out, frcmod_out


# =============================================================================
# PROTEIN PREPARATION
# =============================================================================

def prepare_protein(pdb_path, remove_wat, output_dir='prep'):
    """Prepare protein with pdb4amber"""
    print(f"\n{'='*60}")
    print(f"PROTEIN PREPARATION")
    print(f"{'='*60}")
    print(f"  PDB file: {pdb_path}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open(pdb_path, 'r') as f:
        pdb_content = f.read()
    
    has_kcx = 'KCX' in pdb_content
    has_zn = ' ZN ' in pdb_content or 'ZN2' in pdb_content
    
    print(f"  Contains KCX: {has_kcx}")
    print(f"  Contains Zn: {has_zn}")
    
    protein_out = os.path.join(output_dir, 'protein.pdb')
    
    cmd = ['pdb4amber', '-i', pdb_path, '-o', protein_out]
    if remove_wat == 'yes':
        cmd.append('--dry')
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  WARNING: pdb4amber warning: {result.stderr}")
        shutil.copy(pdb_path, protein_out)
    
    print(f"  Protein prepared: {protein_out}")
    return protein_out, has_kcx, has_zn


# =============================================================================
# TLEAP SYSTEM BUILDING
# =============================================================================

def create_tleap_input(protein_pdb, ligand_mol2, ligand_frcmod, 
                       kcx_frcmod, kcx_lib, has_kcx, has_zn):
    """Create tleap input script"""
    
    ff_protein = 'leaprc.protein.ff19SB' if force_field == 'ff19SB' else 'leaprc.protein.ff14SB'
    
    water_models_map = {
        'tip3p': 'leaprc.water.tip3p',
        'opc': 'leaprc.water.opc',
        'spce': 'leaprc.water.spce',
    }
    ff_water = water_models_map.get(water_model, 'leaprc.water.tip3p')
    
    n_ions = max(1, int(ionic_strength * 60))
    
    script = f"""# tleap input for OpenMM-KCX system
source {ff_protein}
source {ff_water}
source leaprc.gaff2
"""
    
    if has_kcx:
        script += f"""
loadamberparams {kcx_frcmod}
loadoff {kcx_lib}
"""
    
    if has_zn:
        script += """
loadamberparams frcmod.ions234lm_126_tip3p
"""
    
    if ligand_mol2 and os.path.exists(ligand_mol2):
        script += f"""
loadamberparams {ligand_frcmod}
LIG = loadmol2 {ligand_mol2}
"""
    
    script += f"""
mol = loadpdb {protein_pdb}
"""
    
    if ligand_mol2 and os.path.exists(ligand_mol2):
        script += "complex = combine {mol LIG}\n"
    else:
        script += "complex = mol\n"
    
    script += f"solvatebox complex TIP3PBOX {box_size}\n"
    
    script += f"""
addions complex Na+ 0
addions complex Cl- 0
addions complex Na+ {n_ions}
addions complex Cl- {n_ions}

check complex
saveamberparm complex system.prmtop system.inpcrd
savepdb complex system.pdb
quit
"""
    
    with open('tleap.in', 'w') as f:
        f.write(script)
    
    return 'tleap.in'


def run_tleap(tleap_file):
    """Run tleap"""
    print(f"\n{'='*60}")
    print(f"RUNNING TLEAP")
    print(f"{'='*60}")
    
    result = subprocess.run(['tleap', '-f', tleap_file], capture_output=True, text=True)
    print(result.stdout)
    
    if not os.path.exists('system.prmtop') or not os.path.exists('system.inpcrd'):
        print("ERROR: tleap failed!")
        print(result.stderr)
        sys.exit(1)
    
    print("  System built successfully!")
    return True


# =============================================================================
# OPENMM SIMULATION
# =============================================================================

def run_simulation():
    """Run OpenMM MD simulation"""
    print(f"\n{'='*60}")
    print(f"OPENMM SIMULATION")
    print(f"{'='*60}")
    
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit
    
    os.makedirs('out', exist_ok=True)
    
    print("Loading topology...")
    prmtop = app.AmberPrmtopFile('system.prmtop')
    inpcrd = app.AmberInpcrdFile('system.inpcrd')
    
    print("Creating system...")
    nb_methods = {'PME': app.PME, 'NoCutoff': app.NoCutoff, 'CutoffPeriodic': app.CutoffPeriodic}
    constraints_map = {'HBonds': app.HBonds, 'AllBonds': app.AllBonds, 'None': None}
    
    system = prmtop.createSystem(
        nonbondedMethod=nb_methods.get(nonbonded_method, app.PME),
        nonbondedCutoff=nonbonded_cutoff * unit.nanometer,
        constraints=constraints_map.get(constraints, app.HBonds),
        rigidWater=rigid_water == 'yes'
    )
    
    system.addForce(mm.MonteCarloBarostat(
        pressure * unit.bar,
        temperature * unit.kelvin,
        barostat_interval
    ))
    
    print(f"Creating {integrator} integrator...")
    ts = timestep * unit.femtosecond
    temp = temperature * unit.kelvin
    friction = friction_coeff / unit.picosecond
    
    integ = mm.LangevinMiddleIntegrator(temp, friction, ts)
    
    print("Setting up simulation...")
    try:
        platform = mm.Platform.getPlatformByName('CUDA')
        properties = {'CudaPrecision': 'mixed'}
        print("  Using CUDA")
    except:
        try:
            platform = mm.Platform.getPlatformByName('OpenCL')
            properties = {}
            print("  Using OpenCL")
        except:
            platform = mm.Platform.getPlatformByName('CPU')
            properties = {}
            print("  Using CPU")
    
    simulation = app.Simulation(prmtop.topology, system, integ, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    # Minimization
    print(f"\nMinimizing ({minimization_steps} steps)...")
    simulation.minimizeEnergy(maxIterations=minimization_steps)
    
    state = simulation.context.getState(getPositions=True)
    with open('out/minimized.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    # Equilibration
    equil_steps = int(equilibration_time * 1e6 / timestep)
    print(f"\nEquilibrating ({equilibration_time} ns)...")
    
    simulation.reporters.append(app.StateDataReporter(
        'out/equilibration.log', equil_traj_freq,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True
    ))
    simulation.step(equil_steps)
    simulation.reporters.clear()
    
    # Production
    prod_steps = int(production_time * 1e6 / timestep)
    print(f"\nProduction ({production_time} ns)...")
    
    simulation.reporters.append(app.StateDataReporter(
        'out/production.log', prod_traj_freq,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True
    ))
    simulation.reporters.append(app.DCDReporter('out/production.dcd', prod_traj_freq))
    simulation.reporters.append(app.CheckpointReporter('out/checkpoint.chk', checkpoint_freq))
    
    simulation.step(prod_steps)
    
    # Save final
    state = simulation.context.getState(getPositions=True)
    with open('out/final.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    print("\nSimulation complete!")


# =============================================================================
# ANALYSIS
# =============================================================================

def run_analysis():
    """Run trajectory analysis"""
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
        print(f"  Skipping analysis: {e}")
        return
    
    if not os.path.exists('out/production.dcd'):
        print("  No trajectory found")
        return
    
    print("  Loading trajectory...")
    traj = md.load('out/production.dcd', top='system.pdb', stride=step_size)
    print(f"  Frames: {traj.n_frames}")
    
    # RMSD
    print("  Calculating RMSD...")
    rmsd_atoms = traj.topology.select('backbone')
    rmsd = md.rmsd(traj, traj, 0, atom_indices=rmsd_atoms) * 10
    time_ns = traj.time / 1000
    
    np.savetxt('out/rmsd.csv', np.column_stack([time_ns, rmsd]), 
               delimiter=',', header='time_ns,rmsd_A', comments='')
    
    plt.figure(figsize=(10, 6))
    plt.plot(time_ns, rmsd)
    plt.xlabel('Time (ns)')
    plt.ylabel('RMSD (Å)')
    plt.title(f'RMSD - {job_name}')
    plt.savefig('out/rmsd.png', dpi=150)
    plt.close()
    
    # RMSF
    print("  Calculating RMSF...")
    ca_atoms = traj.topology.select('name CA')
    if len(ca_atoms) > 0:
        rmsf = md.rmsf(traj, traj, 0, atom_indices=ca_atoms) * 10
        residues = [traj.topology.atom(i).residue.resSeq for i in ca_atoms]
        
        np.savetxt('out/rmsf.csv', np.column_stack([residues, rmsf]),
                   delimiter=',', header='residue,rmsf_A', comments='')
        
        plt.figure(figsize=(12, 6))
        plt.plot(residues, rmsf)
        plt.xlabel('Residue')
        plt.ylabel('RMSF (Å)')
        plt.title(f'RMSF - {job_name}')
        plt.savefig('out/rmsf.png', dpi=150)
        plt.close()
    
    print("  Analysis complete!")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("  OpenMM MD Simulation with KCX Support")
    print("="*70)
    
    # Print job info
    print(f"\nJob Name: {job_name}")
    print(f"PDB File: {pdb_file}")
    print(f"Ligand File: {ligand_file}")
    print(f"Production Time: {production_time} ns")
    print(f"Temperature: {temperature} K")
    
    # Check inputs directory
    print("\nFiles in inputs/:")
    if os.path.exists('inputs'):
        for f in os.listdir('inputs'):
            fpath = os.path.join('inputs', f)
            size = os.path.getsize(fpath) if os.path.isfile(fpath) else 'DIR'
            print(f"  {f} ({size} bytes)" if isinstance(size, int) else f"  {f}/")
    else:
        print("  inputs/ directory does not exist!")
    
    # Check PDB file exists
    if not os.path.exists(pdb_file):
        print(f"\nERROR: PDB file not found at {pdb_file}")
        print("Tamarind should provide file inputs at inputs/{fieldName}")
        sys.exit(1)
    
    print(f"\n✓ PDB file found: {pdb_file}")
    
    # Check ligand file (optional)
    has_ligand = os.path.exists(ligand_file)
    if has_ligand:
        print(f"✓ Ligand file found: {ligand_file}")
    else:
        print(f"  Ligand file not provided (protein-only simulation)")
    
    # Create output directory
    os.makedirs('out', exist_ok=True)
    
    # Setup KCX parameters
    kcx_frcmod, kcx_lib = write_kcx_parameters('params')
    
    # Prepare ligand if present
    ligand_mol2, ligand_frcmod = None, None
    if has_ligand:
        ligand_mol2, ligand_frcmod = prepare_ligand(ligand_file, ligand_charge, 'prep')
    
    # Prepare protein
    protein_pdb, has_kcx, has_zn = prepare_protein(pdb_file, remove_waters, 'prep')
    
    # Build system with tleap
    tleap_file = create_tleap_input(
        protein_pdb, ligand_mol2, ligand_frcmod,
        kcx_frcmod, kcx_lib, has_kcx, has_zn
    )
    run_tleap(tleap_file)
    
    # Run simulation
    run_simulation()
    
    # Run analysis
    run_analysis()
    
    print("\n" + "="*70)
    print("  ALL DONE!")
    print("="*70)
    print("\nOutputs saved to out/ directory:")
    if os.path.exists('out'):
        for f in os.listdir('out'):
            print(f"  {f}")


if __name__ == '__main__':
    main()
