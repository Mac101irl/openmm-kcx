#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX Support
Optimized for NVIDIA A100 on Tamarind Platform
"""

import os
import sys
import subprocess
import shutil
import traceback

import matplotlib
matplotlib.use('Agg')

# =============================================================================
# TAMARIND PATH CONFIGURATION
# =============================================================================
PDB_FILE_INPUT = "inputs/pdbFile.pdb"
LIGAND_FILE_INPUT = "inputs/ligandFile.sdf"
OUTPUT_DIR = "out"
PREP_DIR = os.path.join(OUTPUT_DIR, "prep")
PARAMS_DIR = os.path.join(OUTPUT_DIR, "params")
KCX_PARAMS_DIR = "/app/kcx_params"

# =============================================================================
# PARAMETERS FROM ENVIRONMENT VARIABLES
# =============================================================================
job_name = os.getenv('JobName', 'openmm_simulation')
ligand_charge = int(os.getenv('ligandCharge', '0'))
force_field = os.getenv('forceField', 'ff19sb').lower()
water_model = os.getenv('waterModel', 'tip3p').lower()
box_size = float(os.getenv('boxSize', '12.0'))
ionic_strength = float(os.getenv('ionicStrength', '0.15'))
minimization_steps = int(os.getenv('minimizationSteps', '10000'))
equilibration_time = float(os.getenv('equilibrationTime', '0.2'))
production_time = float(os.getenv('productionTime', '1.0'))
timestep = float(os.getenv('timestep', '2.0'))
temperature = float(os.getenv('temperature', '310.0'))
pressure = float(os.getenv('pressure', '1.0'))
pH = float(os.getenv('pH', '7.0'))
constraints = os.getenv('constraints', 'HBonds')
prod_traj_freq = int(os.getenv('prodTrajFreq', '5000'))
step_size = int(os.getenv('stepSize', '5'))
cuda_precision = os.getenv('CUDA_PRECISION', 'mixed')
default_platform = os.getenv('OPENMM_DEFAULT_PLATFORM', 'CUDA')

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def run_command(cmd, description):
    """Run subprocess command with error handling"""
    print(f"--> Running: {description}")
    print(f"    Command: {' '.join(cmd)}")
    process = subprocess.Popen(
        cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
    )
    stdout, stderr = process.communicate()
    
    if process.returncode != 0:
        print(f"ERROR in {description}:")
        print(f"STDOUT: {stdout}")
        print(f"STDERR: {stderr}")
        sys.exit(1)
    else: 
        if stdout.strip():
            print(f"    Output: {stdout.strip()[:200]}")
    return process

def print_gpu_info():
    """Print GPU information for debugging"""
    print("\n--> GPU Environment Configuration:")
    gpu_vars = [
        'NVIDIA_VISIBLE_DEVICES', 'CUDA_VISIBLE_DEVICES',
        'OPENMM_DEFAULT_PLATFORM', 'CUDA_CACHE_DISABLE',
        'CUDA_MPS_PIPE_DIRECTORY'
    ]
    for var in gpu_vars:
        value = os.getenv(var, 'not set')
        print(f"    {var}: {value}")

# =============================================================================
# KCX & SYSTEM PREP UTILITIES
# =============================================================================

def get_kcx_parameters():
    """Use pre-bundled KCX parameters from Docker image"""
    frcmod_path = os.path.join(KCX_PARAMS_DIR, 'kcx.frcmod')
    lib_path = os.path.join(KCX_PARAMS_DIR, 'kcx.lib')
    
    if not os.path.exists(frcmod_path):
        print(f"WARNING: KCX frcmod not found at {frcmod_path}")
        print("    Continuing without KCX parameters (protein-only mode)")
        return None, None
    if not os.path.exists(lib_path):
        print(f"WARNING: KCX lib not found at {lib_path}")
        print("    Continuing without KCX parameters (protein-only mode)")
        return None, None
    
    print(f"--> Using KCX parameters from {KCX_PARAMS_DIR}")
    return frcmod_path, lib_path

def validate_kcx_residues(pdb_path):
    """
    Validate KCX (carboxylated lysine) residues in PDB file.
    
    Checks:
    1. Presence of KCX residues
    2. Required carbamate atoms (NZ, HZ, CY, OQ1, OQ2)
    3. Basic geometry validation
    4. pH-dependent warnings
    
    Returns:
        dict: Validation results with KCX info
    """
    print(f"--> Validating KCX residues in {pdb_path}")
    
    # Expected carbamate atoms in KCX (from kcx.lib)
    REQUIRED_CARBAMATE_ATOMS = {'NZ', 'CY', 'OQ1', 'OQ2'}  # HZ may be missing in some PDBs
    OPTIONAL_ATOMS = {'HZ'}  # Hydrogen may not be in crystal structures
    
    # Parse PDB and collect KCX residues
    kcx_residues = {}  # {(chain, resnum): {atom_name: (x, y, z)}}
    
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                res_name = line[17:20].strip()
                if res_name == 'KCX':
                    chain = line[21].strip() or 'A'
                    res_num = int(line[22:26].strip())
                    atom_name = line[12:16].strip()
                    try:
                        x = float(line[30:38])
                        y = float(line[38:46])
                        z = float(line[46:54])
                    except ValueError:
                        x, y, z = 0.0, 0.0, 0.0
                    
                    key = (chain, res_num)
                    if key not in kcx_residues:
                        kcx_residues[key] = {}
                    kcx_residues[key][atom_name] = (x, y, z)
    
    # Report findings
    results = {
        'found': len(kcx_residues) > 0,
        'count': len(kcx_residues),
        'residues': [],
        'warnings': [],
        'errors': []
    }
    
    if not kcx_residues:
        print("    No KCX residues found in PDB")
        print("    NOTE: If your protein has carboxylated lysine, ensure it's labeled as 'KCX'")
        print("          AlphaFold structures may label KCX as 'LYS' - manual editing required")
        results['warnings'].append("No KCX residues found - check if lysines should be carboxylated")
        return results
    
    print(f"    Found {len(kcx_residues)} KCX residue(s)")
    
    for (chain, res_num), atoms in kcx_residues.items():
        res_id = f"KCX-{chain}{res_num}"
        atom_names = set(atoms.keys())
        
        # Check for required carbamate atoms
        missing_required = REQUIRED_CARBAMATE_ATOMS - atom_names
        missing_optional = OPTIONAL_ATOMS - atom_names
        
        res_info = {
            'chain': chain,
            'resnum': res_num,
            'atoms': list(atom_names),
            'valid': True
        }
        
        if missing_required:
            print(f"    ✗ {res_id}: MISSING carbamate atoms: {missing_required}")
            print(f"      This KCX may not be parameterized correctly!")
            results['errors'].append(f"{res_id} missing atoms: {missing_required}")
            res_info['valid'] = False
        else:
            print(f"    ✓ {res_id}: All carbamate atoms present")
            
            # Geometry check - CY-OQ bond length should be ~1.25 Å
            if 'CY' in atoms and 'OQ1' in atoms:
                cy = atoms['CY']
                oq1 = atoms['OQ1']
                dist = ((cy[0]-oq1[0])**2 + (cy[1]-oq1[1])**2 + (cy[2]-oq1[2])**2) ** 0.5
                if 1.15 <= dist <= 1.35:
                    print(f"      ✓ CY-OQ1 bond length: {dist:.2f} Å (expected ~1.25 Å)")
                else:
                    print(f"      ⚠ CY-OQ1 bond length: {dist:.2f} Å (expected ~1.25 Å) - unusual!")
                    results['warnings'].append(f"{res_id} unusual CY-OQ1 distance: {dist:.2f} Å")
        
        if missing_optional:
            print(f"      Note: Missing optional atoms {missing_optional} (will be added by reduce)")
        
        results['residues'].append(res_info)
    
    # pH-specific warnings for KCX
    if pH < 6.0:
        print(f"    ⚠ WARNING: pH {pH} is near/below KCX carbamate pKa (~5.5)")
        print(f"      The carbamate group may be partially protonated at this pH")
        print(f"      Current parameters assume deprotonated (charged) carbamate")
        results['warnings'].append(f"pH {pH} may affect KCX protonation state")
    
    # Summary
    valid_count = sum(1 for r in results['residues'] if r['valid'])
    print(f"    Summary: {valid_count}/{len(kcx_residues)} KCX residues validated successfully")
    
    if results['errors']:
        print(f"    ⚠ {len(results['errors'])} error(s) found - simulation may fail!")
    
    return results

def prepare_ligand(lig_file, lig_charge, output_dir):
    """Prepare ligand with Antechamber and GAFF2"""
    print(f"--> Preparing Ligand: {lig_file}")
    os.makedirs(output_dir, exist_ok=True)
    mol2_out = os.path.join(output_dir, 'ligand.mol2')
    frcmod_out = os.path.join(output_dir, 'ligand.frcmod')
    
    cmd = [
        'antechamber', '-i', lig_file, '-fi', 'sdf',
        '-o', mol2_out, '-fo', 'mol2',
        '-c', 'bcc', '-at', 'gaff2',
        '-nc', str(lig_charge), '-pf', 'y'
    ]
    run_command(cmd, "Antechamber")
    
    cmd = ['parmchk2', '-i', mol2_out, '-f', 'mol2', '-o', frcmod_out, '-s', 'gaff2']
    run_command(cmd, "Parmchk2")
    
    print(f"--> Ligand prepared: {mol2_out}")
    return mol2_out, frcmod_out

def prepare_protein(pdb_path, output_dir):
    """Prepare protein structure with pdb4amber and pH-dependent protonation"""
    print(f"--> Preparing Protein: {pdb_path}")
    print(f"    Target pH: {pH}")
    os.makedirs(output_dir, exist_ok=True)
    
    # Step 1: Clean PDB with pdb4amber (remove waters, alternate conformations)
    cleaned_pdb = os.path.join(output_dir, 'protein_cleaned.pdb')
    cmd = ['pdb4amber', '-i', pdb_path, '-o', cleaned_pdb, '--dry', '--nohyd']
    run_command(cmd, "PDB4AMBER (clean)")
    
    # Step 2: Add hydrogens with reduce at target pH
    reduced_pdb = os.path.join(output_dir, 'protein_reduced.pdb')
    try:
        # reduce uses -build to add hydrogens, -HIS to flip histidines optimally
        cmd = ['reduce', '-build', '-his', cleaned_pdb]
        print(f"--> Running: Add hydrogens with reduce")
        print(f"    Command: {' '.join(cmd)}")
        process = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        stdout, stderr = process.communicate()
        
        if process.returncode == 0:
            with open(reduced_pdb, 'w') as f:
                f.write(stdout)
            print(f"    Hydrogens added successfully")
        else:
            print(f"    WARNING: reduce failed, continuing with cleaned PDB")
            print(f"    STDERR: {stderr[:200]}")
            reduced_pdb = cleaned_pdb
    except FileNotFoundError:
        print(f"    WARNING: reduce not found, skipping hydrogen addition")
        reduced_pdb = cleaned_pdb
    
    # Step 3: Assign histidine protonation states based on pH
    protein_out = os.path.join(output_dir, 'protein_fixed.pdb')
    assign_histidine_protonation(reduced_pdb, protein_out, pH)
    
    # pH warnings
    if pH < 5.0:
        print(f"    WARNING: Low pH ({pH}) - Asp/Glu may need protonation (ASH/GLH)")
        print(f"             Consider manual review of carboxylate residues")
    elif pH > 9.0:
        print(f"    WARNING: High pH ({pH}) - Lys/Cys/Tyr may need deprotonation")
        print(f"             Consider manual review of titratable residues")
    
    # KCX-specific warning
    if pH < 6.0:
        print(f"    WARNING: KCX carbamate is pH-sensitive (pKa ~5.5)")
        print(f"             At pH {pH}, KCX may be partially protonated")
    
    print(f"--> Protein prepared: {protein_out}")
    return protein_out

def assign_histidine_protonation(input_pdb, output_pdb, target_pH):
    """
    Assign histidine protonation states based on target pH.
    
    Histidine pKa ~6.0:
    - pH < 6.0: HIP (doubly protonated, +1 charge)
    - pH 6.0-7.0: Mixed HID/HIE (singly protonated, neutral)
    - pH > 7.0: HIE preferred (epsilon-protonated, neutral)
    
    TLeap recognizes: HIS -> HIE, HID -> HID, HIP -> HIP
    """
    print(f"    Assigning histidine protonation for pH {target_pH}")
    
    # Determine histidine state based on pH
    if target_pH < 6.0:
        his_state = 'HIP'  # Doubly protonated (charged)
        print(f"      pH < 6.0: Using HIP (doubly protonated, +1)")
    elif target_pH < 7.0:
        his_state = 'HID'  # Delta-protonated (neutral) - good for H-bonding
        print(f"      pH 6.0-7.0: Using HID (delta-protonated, neutral)")
    else:
        his_state = 'HIE'  # Epsilon-protonated (neutral) - most common at physiological pH
        print(f"      pH >= 7.0: Using HIE (epsilon-protonated, neutral)")
    
    # Read and modify PDB
    his_count = 0
    with open(input_pdb, 'r') as f_in:
        lines = f_in.readlines()
    
    with open(output_pdb, 'w') as f_out:
        for line in lines:
            if line.startswith(('ATOM', 'HETATM')):
                res_name = line[17:20].strip()
                if res_name == 'HIS':
                    # Replace HIS with appropriate protonation state
                    line = line[:17] + f'{his_state:<3}' + line[20:]
                    his_count += 1
            f_out.write(line)
    
    if his_count > 0:
        # Count unique residues (divide by ~10 atoms per His)
        unique_his = his_count // 10 + (1 if his_count % 10 > 0 else 0)
        print(f"      Renamed ~{unique_his} histidine residue(s) to {his_state}")
    else:
        print(f"      No histidine residues found")
    
    return output_pdb

# =============================================================================
# TLEAP & SYSTEM BUILDING
# =============================================================================

def get_water_box_name(water_model):
    """Get the correct water box name for TLeap"""
    water_boxes = {
        'tip3p': 'TIP3PBOX',
        'tip4pew': 'TIP4PEWBOX',
        'opc': 'OPCBOX',
        'opc3': 'OPC3BOX',
        'spce': 'SPCBOX',
    }
    return water_boxes.get(water_model.lower(), 'TIP3PBOX')

def calculate_ion_pairs(ionic_strength_M, box_padding_A, estimated_protein_size_A=60):
    """
    Calculate number of Na+/Cl- ion pairs needed for target ionic strength.
    
    Args:
        ionic_strength_M: Target ionic strength in Molar (e.g., 0.15 for 150 mM)
        box_padding_A: Box padding in Angstroms
        estimated_protein_size_A: Estimated protein diameter in Angstroms
    
    Returns:
        Number of ion pairs to add (after neutralization)
    """
    if ionic_strength_M <= 0:
        return 0
    
    # Estimate box edge length (protein + padding on each side)
    estimated_box_edge_A = estimated_protein_size_A + 2 * box_padding_A
    
    # Volume in liters: (Angstrom * 1e-8 cm/A)^3 * 1e-3 L/cm^3
    volume_liters = (estimated_box_edge_A * 1e-8) ** 3 * 1e-3
    
    # Number of ion pairs: concentration * volume * Avogadro's number
    # For NaCl: ionic_strength = concentration (since both ions are monovalent)
    avogadro = 6.022e23
    num_ion_pairs = int(ionic_strength_M * volume_liters * avogadro)
    
    # Ensure at least 1 pair if ionic strength is requested
    num_ion_pairs = max(1, num_ion_pairs)
    
    print(f"    Ionic strength calculation:")
    print(f"      Target: {ionic_strength_M*1000:.0f} mM NaCl")
    print(f"      Est. box edge: ~{estimated_box_edge_A:.0f} A")
    print(f"      Est. volume: {volume_liters*1e24:.0f} A^3")
    print(f"      Ion pairs to add: {num_ion_pairs}")
    
    return num_ion_pairs

def build_system(protein_pdb, lig_mol2, lig_frcmod, kcx_frcmod, kcx_lib, output_dir):
    """Build solvated system with TLeap"""
    print("--> Building System with TLeap")
    os.makedirs(output_dir, exist_ok=True)
    
    tleap_in = os.path.join(output_dir, 'tleap.in')
    prmtop = os.path.join(output_dir, 'system.prmtop')
    inpcrd = os.path.join(output_dir, 'system.inpcrd')
    water_box = get_water_box_name(water_model)
    
    script = f"""# TLeap input script for {job_name}
source leaprc.protein.{force_field}
source leaprc.water.{water_model}
source leaprc.gaff2
"""
    
    if kcx_frcmod and kcx_lib:
        script += f"""
# Load KCX parameters
loadamberparams {kcx_frcmod}
loadoff {kcx_lib}
"""
    
    script += f"""
# Load protein
mol = loadpdb {protein_pdb}
"""
    
    if lig_mol2 and os.path.exists(lig_mol2):
        script += f"""
# Load ligand
loadamberparams {lig_frcmod}
LIG = loadmol2 {lig_mol2}
system = combine {{mol LIG}}
"""
    else:
        script += "system = mol\n"
    
    # Calculate ion pairs for ionic strength
    num_ion_pairs = calculate_ion_pairs(ionic_strength, box_size)
    
    script += f"""
# Solvate and add ions
solvatebox system {water_box} {box_size}
# First neutralize the system
addions system Na+ 0
addions system Cl- 0
"""
    
    # Add additional ions for ionic strength (if requested)
    if num_ion_pairs > 0:
        script += f"""
# Add ions for {ionic_strength*1000:.0f} mM ionic strength
addions system Na+ {num_ion_pairs}
addions system Cl- {num_ion_pairs}
"""
    
    script += f"""check system
saveamberparm system {prmtop} {inpcrd}
savepdb system {os.path.join(output_dir, 'system_solvated.pdb')}
quit
"""
    
    with open(tleap_in, 'w') as f:
        f.write(script)
    
    print(f"--> TLeap script written to: {tleap_in}")
    run_command(['tleap', '-f', tleap_in], "TLeap")
    
    if not os.path.exists(prmtop) or not os.path.exists(inpcrd):
        print("ERROR: TLeap failed to generate topology files")
        sys.exit(1)
    
    prmtop_size = os.path.getsize(prmtop) / 1024
    inpcrd_size = os.path.getsize(inpcrd) / 1024
    print(f"--> System built: {prmtop} ({prmtop_size:.1f} KB), {inpcrd} ({inpcrd_size:.1f} KB)")
    
    return prmtop, inpcrd

# =============================================================================
# OPENMM SIMULATION
# =============================================================================

def get_constraint_type(constraint_str):
    """Convert constraint string to OpenMM constraint type"""
    import openmm.app as app
    constraint_map = {
        'hbonds': app.HBonds,
        'allbonds': app.AllBonds,
        'hangles': app.HAngles,
        'none': None,
    }
    return constraint_map.get(constraint_str.lower(), app.HBonds)

def run_simulation(prmtop_path, inpcrd_path):
    """Run OpenMM simulation optimized for NVIDIA A100"""
    print(f"\n{'='*60}\nSTARTING OPENMM SIMULATION\n{'='*60}")
    
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    print("--> Available OpenMM Platforms:")
    for i in range(mm.Platform.getNumPlatforms()):
        platform = mm.Platform.getPlatform(i)
        print(f"    {i}: {platform.getName()}")

    print("--> Loading topology and coordinates")
    prmtop = app.AmberPrmtopFile(prmtop_path)
    inpcrd = app.AmberInpcrdFile(inpcrd_path)
    print(f"    Atoms: {prmtop.topology.getNumAtoms()}")
    print(f"    Residues: {prmtop.topology.getNumResidues()}")

    constraint_type = get_constraint_type(constraints)
    
    print("--> Creating system")
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=constraint_type,
        rigidWater=True,
        hydrogenMass=None
    )
    
    system.addForce(mm.MonteCarloBarostat(
        pressure*unit.bar,
        temperature*unit.kelvin,
        25
    ))

    print("--> Setting up integrator")
    integrator = mm.LangevinMiddleIntegrator(
        temperature*unit.kelvin,
        1.0/unit.picosecond,
        timestep*unit.femtosecond
    )

    print("--> Configuring platform")
    platform = None
    properties = {}
    
    for platform_name in [default_platform, 'CUDA', 'OpenCL', 'CPU']: 
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            if platform_name == 'CUDA':
                properties = {'Precision': cuda_precision}
            elif platform_name == 'OpenCL':
                properties = {'Precision': cuda_precision}
            print(f"--> [SUCCESS] Using {platform_name} Platform")
            break
        except Exception as e:
            print(f"--> [INFO] {platform_name} not available: {e}")
            continue
    
    if platform is None:
        print("--> [ERROR] No suitable platform found!")
        sys.exit(1)

    if properties: 
        sim = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    else:
        sim = app.Simulation(prmtop.topology, system, integrator, platform)
    
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors:
        sim.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    print(f"--> [PLATFORM] {platform.getName()}")
    if platform.getName() == 'CUDA':
        try:
            print(f"    Device: {platform.getPropertyValue(sim.context, 'DeviceName')}")
            print(f"    Precision: {platform.getPropertyValue(sim.context, 'Precision')}")
        except Exception as e:
            print(f"    (Could not retrieve device info: {e})")

    print(f"--> Minimizing energy ({minimization_steps} steps)...")
    sim.minimizeEnergy(maxIterations=minimization_steps)
    state = sim.context.getState(getEnergy=True)
    print(f"--> Minimization complete. Potential Energy: {state.getPotentialEnergy()}")

    equil_steps = int((equilibration_time * 1e6) / timestep)
    print(f"--> Running equilibration ({equil_steps} steps, {equilibration_time} ns)...")
    sim.step(equil_steps)
    print("--> Equilibration complete")

    total_steps = int((production_time * 1e6) / timestep)
    print(f"--> Running Production MD:")
    print(f"    Duration: {production_time} ns, Steps: {total_steps}")
    
    dcd_file = os.path.join(OUTPUT_DIR, 'trajectory.dcd')
    log_file = os.path.join(OUTPUT_DIR, 'simulation.log')
    checkpoint_file = os.path.join(OUTPUT_DIR, 'checkpoint.chk')
    
    sim.reporters.append(app.DCDReporter(dcd_file, prod_traj_freq))
    sim.reporters.append(app.StateDataReporter(
        log_file, prod_traj_freq,
        step=True, time=True, potentialEnergy=True, kineticEnergy=True,
        totalEnergy=True, temperature=True, volume=True, density=True, speed=True
    ))
    sim.reporters.append(app.CheckpointReporter(checkpoint_file, prod_traj_freq * 10))
    sim.reporters.append(app.StateDataReporter(
        sys.stdout, prod_traj_freq * 10,
        step=True, time=True, speed=True, remainingTime=True, totalSteps=total_steps
    ))
    
    print("--> Starting production run...")
    sim.step(total_steps)
    
    final_pdb = os.path.join(OUTPUT_DIR, 'final_structure.pdb')
    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    with open(final_pdb, 'w') as f:
        app.PDBFile.writeFile(sim.topology, state.getPositions(), f)
    
    sim.saveCheckpoint(os.path.join(OUTPUT_DIR, 'final_checkpoint.chk'))
    
    print(f"\n{'='*60}")
    print("--> Simulation Finished Successfully")
    print(f"    Trajectory: {dcd_file}")
    print(f"    Final structure: {final_pdb}")
    print(f"{'='*60}\n")
    
    return dcd_file, final_pdb

# =============================================================================
# ANALYSIS
# =============================================================================

def run_analysis(dcd_path, pdb_path):
    """Run trajectory analysis (RMSD)"""
    print("--> Running Trajectory Analysis")
    
    try:
        import mdtraj as md
        import matplotlib.pyplot as plt
        import numpy as np
        
        print(f"    Loading trajectory: {dcd_path}")
        traj = md.load(dcd_path, top=pdb_path)
        print(f"    Loaded: {traj.n_frames} frames, {traj.n_atoms} atoms")
        
        print("    Calculating backbone RMSD...")
        backbone = traj.topology.select('backbone')
        if len(backbone) == 0:
            print("    WARNING: No backbone atoms, using heavy atoms")
            backbone = traj.topology.select('not element H')
        
        rmsd = md.rmsd(traj, traj[0], atom_indices=backbone)
        rmsd_angstrom = rmsd * 10
        time_ns = traj.time / 1000
        
        plt.figure(figsize=(10, 6))
        plt.plot(time_ns, rmsd_angstrom, linewidth=1.5, color='#2E86AB')
        plt.xlabel('Time (ns)', fontsize=12)
        plt.ylabel('Backbone RMSD (Angstrom)', fontsize=12)
        plt.title(f'Backbone RMSD vs Time\nMean: {rmsd_angstrom.mean():.2f} A, Max: {rmsd_angstrom.max():.2f} A')
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        rmsd_plot = os.path.join(OUTPUT_DIR, 'rmsd_analysis.png')
        plt.savefig(rmsd_plot, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"--> RMSD Analysis: Mean={rmsd_angstrom.mean():.2f}A, Max={rmsd_angstrom.max():.2f}A")
        print(f"    Plot saved: {rmsd_plot}")
        
        rmsd_data = os.path.join(OUTPUT_DIR, 'rmsd_data.txt')
        np.savetxt(rmsd_data, np.column_stack([time_ns, rmsd_angstrom]),
                   header='Time(ns) RMSD(Angstrom)', fmt='%.4f')
        
    except ImportError as e:
        print(f"--> Analysis skipped: Missing library ({e})")
    except Exception as e: 
        print(f"--> Analysis failed: {e}")
        traceback.print_exc()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    try:
        print(f"\n{'='*60}")
        print("OpenMM MD Simulation with KCX Support")
        print(f"Job Name: {job_name}")
        print(f"{'='*60}\n")
        
        print_gpu_info()
        
        print("\n--> Simulation Parameters:")
        print(f"    Force Field: {force_field}, Water: {water_model}")
        print(f"    Box: {box_size}A, Temp: {temperature}K, Press: {pressure}bar")
        print(f"    Timestep: {timestep}fs, Equil: {equilibration_time}ns, Prod: {production_time}ns")
        
        print("\n--> Setting up directories")
        for d in [OUTPUT_DIR, PREP_DIR, PARAMS_DIR]: 
            os.makedirs(d, exist_ok=True)
            print(f"    Created: {d}")
        
        print("\n--> Validating inputs")
        if not os.path.exists(PDB_FILE_INPUT):
            print(f"FATAL ERROR: Missing PDB file at {PDB_FILE_INPUT}")
            if os.path.exists('inputs'):
                print(f"    Contents of inputs/: {os.listdir('inputs')}")
            sys.exit(1)
        print(f"    Found PDB: {PDB_FILE_INPUT}")
        
        has_ligand = os.path.exists(LIGAND_FILE_INPUT)
        if has_ligand: 
            print(f"    Found Ligand: {LIGAND_FILE_INPUT}")
        else:
            print("    No ligand file (protein-only simulation)")
        
        print("\n--> Loading KCX parameters")
        kcx_frcmod, kcx_lib = get_kcx_parameters()
        
        # Validate KCX residues in input PDB
        print("\n--> Checking for KCX residues")
        kcx_validation = validate_kcx_residues(PDB_FILE_INPUT)
        if kcx_validation['errors']:
            print("\n    WARNING: KCX validation errors detected!")
            print("    The simulation will continue, but results may be incorrect.")
            print("    Consider fixing KCX atom names to match kcx.lib before proceeding.")
        
        lig_mol2, lig_frcmod = None, None
        if has_ligand: 
            print("\n--> Preparing ligand")
            lig_mol2, lig_frcmod = prepare_ligand(LIGAND_FILE_INPUT, ligand_charge, PREP_DIR)
        
        print("\n--> Preparing protein")
        fixed_protein = prepare_protein(PDB_FILE_INPUT, PREP_DIR)
        
        print("\n--> Building system topology")
        prmtop, inpcrd = build_system(
            fixed_protein, lig_mol2, lig_frcmod,
            kcx_frcmod, kcx_lib, OUTPUT_DIR
        )
        
        print("\n--> Running OpenMM simulation")
        traj_path, final_pdb = run_simulation(prmtop, inpcrd)
        
        print("\n--> Analyzing trajectory")
        run_analysis(traj_path, final_pdb)
        
        print(f"\n{'='*60}")
        print("SIMULATION COMPLETE")
        print(f"{'='*60}")
        print(f"Job: {job_name}")
        print(f"Results: {OUTPUT_DIR}/")
        print(f"{'='*60}\n")
        
    except Exception as e: 
        print(f"\nCRITICAL FAILURE: {e}")
        traceback.print_exc()
        sys.exit(1)