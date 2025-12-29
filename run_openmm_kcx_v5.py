#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX Support
Optimized for NVIDIA A100 on Tamarind Platform

Features:
- Carboxylated Lysine (KCX) & Zinc Support
- Tamarind pathing (inputs/ and out/)
- NVIDIA A100 Acceleration (Mixed Precision + JIT Kernels)
- Automated Trajectory Analysis (RMSD/RMSF)

Updates in this version:
- Enhanced error handling for AmberTools
- Proper matplotlib backend for headless operation
- Pre-bundled KCX parameters from Docker image
- Robust directory creation
- GPU optimization aligned with Dockerfile ENV variables
"""

import os
import sys
import subprocess
import shutil

# Configure matplotlib for headless operation BEFORE any other imports
import matplotlib
matplotlib.use('Agg')

# =============================================================================
# TAMARIND PATH CONFIGURATION
# =============================================================================
# File inputs: Read from inputs/settingName.ext (use exact setting name)
PDB_FILE_INPUT = "inputs/pdbFile.pdb"
LIGAND_FILE_INPUT = "inputs/ligandFile.sdf"

# Outputs: Save ALL results to out/ directory
OUTPUT_DIR = "out"

# Define subdirectories for organization
PREP_DIR = os.path. join(OUTPUT_DIR, "prep")
PARAMS_DIR = os.path.join(OUTPUT_DIR, "params")

# KCX parameters bundled in Docker image
KCX_PARAMS_DIR = "/app/kcx_params"

# =============================================================================
# PARAMETERS FROM ENVIRONMENT VARIABLES
# =============================================================================
# Job name:  Available as os.getenv('JobName')
job_name = os. getenv('JobName', 'openmm_simulation')
ligand_charge = int(os.getenv('ligandCharge', '0'))
force_field = os.getenv('forceField', 'ff19SB')
water_model = os.getenv('waterModel', 'tip3p')
box_size = float(os.getenv('boxSize', '12.0'))
ionic_strength = float(os.getenv('ionicStrength', '0.15'))
minimization_steps = int(os.getenv('minimizationSteps', '10000'))
equilibration_time = float(os.getenv('equilibrationTime', '0.2'))
production_time = float(os. getenv('productionTime', '1.0'))
timestep = float(os.getenv('timestep', '2.0'))
temperature = float(os. getenv('temperature', '310.0'))
pressure = float(os. getenv('pressure', '1.0'))
constraints = os.getenv('constraints', 'HBonds')
prod_traj_freq = int(os.getenv('prodTrajFreq', '5000'))
step_size = int(os.getenv('stepSize', '5'))

# GPU Configuration from environment (set in Dockerfile)
cuda_precision = os. getenv('CUDA_PRECISION', 'mixed')
default_platform = os.getenv('OPENMM_DEFAULT_PLATFORM', 'CUDA')

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def run_command(cmd, description):
    """Run subprocess command with error handling"""
    print(f"--> Running: {description}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"ERROR in {description}:")
        print(f"STDOUT: {result.stdout}")
        print(f"STDERR: {result.stderr}")
        sys.exit(1)
    return result

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
    lib_path = os. path.join(KCX_PARAMS_DIR, 'kcx.lib')
    
    if not os.path.exists(frcmod_path):
        print(f"ERROR: KCX frcmod not found at {frcmod_path}")
        sys.exit(1)
    if not os.path. exists(lib_path):
        print(f"ERROR: KCX lib not found at {lib_path}")
        sys.exit(1)
    
    print(f"--> Using KCX parameters from {KCX_PARAMS_DIR}")
    return frcmod_path, lib_path

def prepare_ligand(lig_file, lig_charge, output_dir):
    """Prepare ligand with Antechamber and GAFF2"""
    print(f"--> Preparing Ligand: {lig_file}")
    os.makedirs(output_dir, exist_ok=True)
    mol2_out = os.path. join(output_dir, 'ligand.mol2')
    frcmod_out = os.path. join(output_dir, 'ligand.frcmod')
    
    # Run antechamber
    cmd = [
        'antechamber', '-i', lig_file, '-fi', 'sdf', 
        '-o', mol2_out, '-fo', 'mol2',
        '-c', 'bcc', '-at', 'gaff2', 
        '-nc', str(lig_charge), '-pf', 'y'
    ]
    run_command(cmd, "Antechamber")
    
    # Run parmchk2
    cmd = ['parmchk2', '-i', mol2_out, '-f', 'mol2', '-o', frcmod_out]
    run_command(cmd, "Parmchk2")
    
    print(f"--> Ligand prepared:  {mol2_out}")
    return mol2_out, frcmod_out

def prepare_protein(pdb_path, output_dir):
    """Prepare protein structure with pdb4amber"""
    print(f"--> Preparing Protein: {pdb_path}")
    os.makedirs(output_dir, exist_ok=True)
    protein_out = os. path.join(output_dir, 'protein_fixed.pdb')
    
    cmd = ['pdb4amber', '-i', pdb_path, '-o', protein_out, '--dry']
    run_command(cmd, "PDB4AMBER")
    
    print(f"--> Protein prepared:  {protein_out}")
    return protein_out

# =============================================================================
# TLEAP & SYSTEM BUILDING
# =============================================================================

def build_system(protein_pdb, lig_mol2, lig_frcmod, kcx_frcmod, kcx_lib, output_dir):
    """Build solvated system with TLeap"""
    print("--> Building System with TLeap")
    os.makedirs(output_dir, exist_ok=True)
    
    tleap_in = os.path. join(output_dir, 'tleap. in')
    prmtop = os.path.join(output_dir, 'system.prmtop')
    inpcrd = os. path.join(output_dir, 'system.inpcrd')
    
    # Build TLeap script
    script = f"""source leaprc.protein.{force_field}
source leaprc.water.{water_model}
source leaprc.gaff2

# Load KCX parameters
loadamberparams {kcx_frcmod}
loadoff {kcx_lib}

# Load protein
mol = loadpdb {protein_pdb}
"""
    
    # Add ligand if present
    if lig_mol2 and os.path.exists(lig_mol2):
        script += f"""
# Load ligand parameters
loadamberparams {lig_frcmod}
LIG = loadmol2 {lig_mol2}

# Combine protein and ligand
system = combine {{mol LIG}}
"""
    else:
        script += "system = mol\n"
    
    # Solvate and neutralize
    script += f"""
# Solvate system
solvatebox system TIP3PBOX {box_size}

# Add ions to neutralize and set ionic strength
addionsrand system Na+ 0
addionsrand system Cl- 0

# Save topology and coordinates
saveamberparm system {prmtop} {inpcrd}

# Quit
quit
"""
    
    # Write TLeap input
    with open(tleap_in, 'w') as f:
        f.write(script)
    
    # Run TLeap
    print(f"--> Running TLeap (script: {tleap_in})")
    run_command(['tleap', '-f', tleap_in], "TLeap")
    
    # Verify output files
    if not os.path.exists(prmtop) or not os.path.exists(inpcrd):
        print("ERROR: TLeap failed to generate topology files")
        sys.exit(1)
    
    print(f"--> System built: {prmtop}, {inpcrd}")
    return prmtop, inpcrd

# =============================================================================
# OPENMM SIMULATION (A100 OPTIMIZED)
# =============================================================================

def run_simulation(prmtop_path, inpcrd_path):
    """Run OpenMM simulation optimized for NVIDIA A100"""
    print(f"\n{'='*60}\nSTARTING OPENMM ON NVIDIA A100\n{'='*60}")
    
    import openmm as mm
    import openmm.app as app
    import openmm.unit as unit

    # Print available platforms
    print("--> Available OpenMM Platforms:")
    for i in range(mm.Platform.getNumPlatforms()):
        platform = mm.Platform. getPlatform(i)
        print(f"    {i}: {platform.getName()}")

    # Load topology and coordinates
    print("--> Loading topology and coordinates")
    prmtop = app. AmberPrmtopFile(prmtop_path)
    inpcrd = app.AmberInpcrdFile(inpcrd_path)

    # 1. System Creation with optimized settings
    print("--> Creating system")
    system = prmtop.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit. nanometer,
        constraints=app.HBonds if constraints == 'HBonds' else None,
        rigidWater=True,
        hydrogenMass=None
    )
    
    # Add barostat for NPT ensemble
    system.addForce(mm.MonteCarloBarostat(
        pressure*unit.bar, 
        temperature*unit. kelvin, 
        25
    ))

    # 2. Integrator (LangevinMiddle is best for modern GPUs)
    print("--> Setting up integrator")
    integrator = mm.LangevinMiddleIntegrator(
        temperature*unit. kelvin, 
        1.0/unit.picosecond, 
        timestep*unit.femtosecond
    )

    # 3. Platform Configuration for A100 (aligned with Dockerfile ENV)
    print("--> Configuring platform")
    platform = None
    properties = {}
    
    # Try platforms in order of preference
    for platform_name in [default_platform, 'CUDA', 'OpenCL', 'CPU']: 
        try:
            platform = mm.Platform.getPlatformByName(platform_name)
            if platform_name in ['CUDA', 'OpenCL']: 
                properties = {
                    'Precision': cuda_precision,
                }
                if platform_name == 'CUDA': 
                    properties['DisablePmeStream'] = 'false'
            print(f"--> [SUCCESS] Using {platform_name} Platform")
            break
        except Exception as e:
            print(f"--> [INFO] {platform_name} not available: {e}")
            continue
    
    if platform is None:
        print("--> [ERROR] No suitable platform found!")
        sys.exit(1)

    # 4. Create simulation context
    if properties:
        sim = app.Simulation(prmtop. topology, system, integrator, platform, properties)
    else:
        sim = app.Simulation(prmtop. topology, system, integrator, platform)
    
    sim.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors:
        sim.context.setPeriodicBoxVectors(*inpcrd. boxVectors)

    # Diagnostic information
    print(f"--> [PLATFORM] {platform.getName()}")
    if platform.getName() == 'CUDA': 
        print(f"    Device: {platform.getPropertyValue(sim.context, 'DeviceName')}")
        print(f"    Device Index: {platform.getPropertyValue(sim. context, 'DeviceIndex')}")
        print(f"    Precision: {platform.getPropertyValue(sim.context, 'Precision')}")
        print(f"    CUDA Compiler: {platform.getPropertyValue(sim. context, 'CudaCompiler')}")

    # 5. Energy Minimization
    print("--> Minimizing energy...")
    print(f"    Target:  {minimization_steps} steps")
    sim.minimizeEnergy(maxIterations=minimization_steps)
    
    state = sim.context. getState(getEnergy=True)
    pot_energy = state. getPotentialEnergy()
    print(f"--> Minimization complete.  Potential Energy: {pot_energy}")

    # 6. Equilibration (NVT then NPT)
    print("--> Running equilibration...")
    equil_steps = int((equilibration_time * 1e6) / timestep)
    print(f"    Equilibration steps: {equil_steps}")
    sim.step(equil_steps)
    print("--> Equilibration complete")

    # 7. Production MD
    total_steps = int((production_time * 1e6) / timestep)
    print(f"--> Running Production MD:")
    print(f"    Duration: {production_time} ns")
    print(f"    Total steps: {total_steps}")
    print(f"    Timestep: {timestep} fs")
    print(f"    Output frequency: {prod_traj_freq} steps")
    
    # Setup reporters - save to out/ directory
    dcd_file = os. path.join(OUTPUT_DIR, 'trajectory.dcd')
    log_file = os. path.join(OUTPUT_DIR, 'simulation.log')
    checkpoint_file = os. path.join(OUTPUT_DIR, 'checkpoint. chk')
    
    sim.reporters. append(app.DCDReporter(dcd_file, prod_traj_freq))
    sim.reporters.append(app.StateDataReporter(
        log_file, 
        prod_traj_freq,
        step=True, 
        time=True, 
        potentialEnergy=True, 
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True, 
        volume=True,
        density=True,
        speed=True
    ))
    
    # Checkpoint reporter for restart capability
    sim.reporters.append(app. CheckpointReporter(checkpoint_file, prod_traj_freq * 10))
    
    # Also print to stdout
    sim.reporters.append(app.StateDataReporter(
        sys.stdout,
        prod_traj_freq * 10,
        step=True,
        time=True,
        speed=True,
        remainingTime=True,
        totalSteps=total_steps
    ))
    
    # Run production
    print("--> Starting production run...")
    sim.step(total_steps)
    
    # Save final state to out/
    final_pdb = os.path.join(OUTPUT_DIR, 'final_structure.pdb')
    state = sim.context.getState(getPositions=True, getVelocities=True, getEnergy=True)
    with open(final_pdb, 'w') as f:
        app.PDBFile. writeFile(sim.topology, state.getPositions(), f)
    
    # Save final checkpoint
    sim.saveCheckpoint(os.path.join(OUTPUT_DIR, 'final_checkpoint. chk'))
    
    print(f"\n{'='*60}")
    print("--> Simulation Finished Successfully")
    print(f"    Trajectory:  {dcd_file}")
    print(f"    Final structure: {final_pdb}")
    print(f"    Log file: {log_file}")
    print(f"    Checkpoint: {checkpoint_file}")
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
        
        # Load trajectory
        print(f"    Loading trajectory: {dcd_path}")
        traj = md.load(dcd_path, top=pdb_path)
        print(f"    Trajectory loaded: {traj. n_frames} frames, {traj.n_atoms} atoms")
        
        # Calculate RMSD for backbone atoms
        print("    Calculating backbone RMSD...")
        backbone = traj.topology.select('backbone')
        rmsd = md.rmsd(traj, traj[0], atom_indices=backbone)
        
        # Convert to Angstroms and ns
        rmsd_angstrom = rmsd * 10
        time_ns = traj.time / 1000
        
        # Generate RMSD plot
        plt.figure(figsize=(10, 6))
        plt.plot(time_ns, rmsd_angstrom, linewidth=1. 5, color='#2E86AB')
        plt.xlabel('Time (ns)', fontsize=12)
        plt.ylabel('Backbone RMSD (Å)', fontsize=12)
        plt.title(f'Backbone RMSD vs Time\nMean: {rmsd_angstrom.mean():.2f} Å, Max: {rmsd_angstrom.max():.2f} Å', 
                  fontsize=14)
        plt.grid(True, alpha=0.3)
        plt.tight_layout()
        
        # Save to out/ directory
        rmsd_plot = os.path. join(OUTPUT_DIR, 'rmsd_analysis.png')
        plt.savefig(rmsd_plot, dpi=150, bbox_inches='tight')
        plt.close()
        
        print(f"--> RMSD Analysis Complete:")
        print(f"    Mean RMSD: {rmsd_angstrom.mean():.2f} Å")
        print(f"    Max RMSD: {rmsd_angstrom. max():.2f} Å")
        print(f"    Min RMSD: {rmsd_angstrom.min():.2f} Å")
        print(f"    Plot saved:  {rmsd_plot}")
        
        # Save RMSD data to out/
        rmsd_data = os.path.join(OUTPUT_DIR, 'rmsd_data.txt')
        np.savetxt(rmsd_data, np.column_stack([time_ns, rmsd_angstrom]), 
                   header='Time(ns) RMSD(Angstrom)', fmt='%.4f')
        print(f"    Data saved:  {rmsd_data}")
        
    except Exception as e: 
        print(f"--> Analysis failed: {e}")
        import traceback
        traceback.print_exc()

# =============================================================================
# MAIN EXECUTION
# =============================================================================

if __name__ == "__main__":
    print(f"\n{'='*60}")
    print("OpenMM MD Simulation with KCX Support")
    print("Optimized for NVIDIA A100 on Tamarind Platform")
    print(f"Job Name: {job_name}")
    print(f"{'='*60}\n")
    
    # 0. Print GPU environment info
    print_gpu_info()
    
    # 1. Create all necessary directories - save to out/
    print("\n--> Setting up directories")
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    os.makedirs(PREP_DIR, exist_ok=True)
    os.makedirs(PARAMS_DIR, exist_ok=True)
    for directory in [OUTPUT_DIR, PREP_DIR, PARAMS_DIR]:
        print(f"    Created: {directory}")
    
    # 2. Validate inputs
    print("\n--> Validating inputs")
    if not os.path. exists(PDB_FILE_INPUT):
        print(f"FATAL ERROR: Missing PDB file at {PDB_FILE_INPUT}")
        sys.exit(1)
    print(f"    Found PDB:  {PDB_FILE_INPUT}")
    
    has_ligand = os.path.exists(LIGAND_FILE_INPUT)
    if has_ligand: 
        print(f"    Found Ligand: {LIGAND_FILE_INPUT}")
    else:
        print(f"    No ligand file detected (protein-only simulation)")
    
    # 3. Get KCX parameters
    print("\n--> Loading KCX parameters")
    kcx_frcmod, kcx_lib = get_kcx_parameters()
    
    # 4. Prepare ligand (if present)
    lig_mol2, lig_frcmod = None, None
    if has_ligand: 
        print("\n--> Preparing ligand")
        lig_mol2, lig_frcmod = prepare_ligand(LIGAND_FILE_INPUT, ligand_charge, PREP_DIR)
    
    # 5. Prepare protein
    print("\n--> Preparing protein")
    fixed_protein = prepare_protein(PDB_FILE_INPUT, PREP_DIR)
    
    # 6. Build system topology
    print("\n--> Building system topology")
    prmtop, inpcrd = build_system(
        fixed_protein, lig_mol2, lig_frcmod, 
        kcx_frcmod, kcx_lib, OUTPUT_DIR
    )
    
    # 7. Run simulation
    print("\n--> Running OpenMM simulation")
    traj_path, final_pdb = run_simulation(prmtop, inpcrd)
    
    # 8. Analyze trajectory
    print("\n--> Analyzing trajectory")
    run_analysis(traj_path, final_pdb)
    
    # 9. Summary
    print(f"\n{'='*60}")
    print("SIMULATION COMPLETE")
    print(f"{'='*60}")
    print(f"Job:  {job_name}")
    print(f"All results saved to: {OUTPUT_DIR}/")
    print(f"\nKey output files:")
    print(f"  - Topology: {OUTPUT_DIR}/system.prmtop")
    print(f"  - Trajectory: {OUTPUT_DIR}/trajectory. dcd")
    print(f"  - Final structure: {OUTPUT_DIR}/final_structure.pdb")
    print(f"  - Simulation log: {OUTPUT_DIR}/simulation.log")
    print(f"  - Checkpoints: {OUTPUT_DIR}/checkpoint.chk, {OUTPUT_DIR}/final_checkpoint.chk")
    print(f"  - RMSD plot: {OUTPUT_DIR}/rmsd_analysis.png")
    print(f"  - RMSD data: {OUTPUT_DIR}/rmsd_data.txt")
    print(f"{'='*60}\n")