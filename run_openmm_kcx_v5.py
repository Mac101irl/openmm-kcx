#!/usr/bin/env python
"""
OpenMM MD Simulation Script with KCX (N6-carboxylysine) Support
For Tamarind Bio Platform

This script follows Tamarind conventions:
- Input files: inputs/<settingName> (with or without extension)
- Settings: os.getenv('<settingName>') or CLI args
- Output files: out/<filename>
"""

import os
import sys
import argparse
import shutil

def find_input_file(setting_name, extensions=None):
    """
    Find an input file in Tamarind's inputs/ directory.
    Tamarind may upload files with or without extensions.
    
    Args:
        setting_name: The base name (e.g., 'pdbFile')
        extensions: List of extensions to try (e.g., ['.pdb', '.cif'])
    
    Returns:
        Path to the found file, or None if not found
    """
    if extensions is None:
        extensions = []
    
    base_path = f'inputs/{setting_name}'
    
    # First check if file exists without extension (Tamarind default)
    if os.path.exists(base_path) and os.path.isfile(base_path):
        print(f"Found input file: {base_path}")
        return base_path
    
    # Then check with each extension
    for ext in extensions:
        path_with_ext = f'{base_path}{ext}'
        if os.path.exists(path_with_ext) and os.path.isfile(path_with_ext):
            print(f"Found input file: {path_with_ext}")
            return path_with_ext
    
    # List what's in the inputs directory for debugging
    if os.path.exists('inputs'):
        print(f"Contents of inputs/ directory: {os.listdir('inputs')}")
    else:
        print("Warning: inputs/ directory does not exist")
    
    return None

def get_settings():
    """Get settings from environment variables (Tamarind convention) or CLI args."""
    parser = argparse.ArgumentParser(description='OpenMM MD with KCX support')
    parser.add_argument('--pdbFile', help='Input PDB file path')
    parser.add_argument('--ligandFile', help='Ligand file path (optional)')
    parser.add_argument('--productionTime', type=float, help='Production time in ns')
    parser.add_argument('--temperature', type=float, default=300, help='Temperature in K')
    parser.add_argument('--timestep', type=float, default=2.0, help='Timestep in fs')
    parser.add_argument('--friction', type=float, default=1.0, help='Friction coefficient')
    parser.add_argument('--reportInterval', type=int, default=10000, help='Report interval')
    
    args = parser.parse_args()
    
    # Priority: CLI args > environment variables > defaults
    settings = {
        'pdbFile': args.pdbFile or os.getenv('pdbFile'),
        'ligandFile': args.ligandFile or os.getenv('ligandFile'),
        'productionTime': args.productionTime or float(os.getenv('productionTime', '1.0')),
        'temperature': args.temperature or float(os.getenv('temperature', '300')),
        'timestep': args.timestep or float(os.getenv('timestep', '2.0')),
        'friction': args.friction or float(os.getenv('friction', '1.0')),
        'reportInterval': args.reportInterval or int(os.getenv('reportInterval', '10000')),
    }
    
    return settings

def prepare_input_files(settings):
    """
    Locate and prepare input files from Tamarind's inputs/ directory.
    Returns paths to the actual files.
    """
    # Find PDB file
    pdb_path = None
    if settings['pdbFile']:
        # If a path was provided, check if it's the Tamarind pattern
        if settings['pdbFile'].startswith('inputs/'):
            # Extract setting name from path
            setting_name = os.path.basename(settings['pdbFile'])
            pdb_path = find_input_file(setting_name, ['.pdb', '.cif', '.ent'])
        else:
            pdb_path = settings['pdbFile']
    else:
        # Try to find pdbFile in inputs/
        pdb_path = find_input_file('pdbFile', ['.pdb', '.cif', '.ent'])
    
    if not pdb_path:
        raise FileNotFoundError(
            f"Could not find PDB file. Checked: inputs/pdbFile and variants. "
            f"Settings pdbFile value: {settings['pdbFile']}"
        )
    
    # Find ligand file (optional)
    ligand_path = None
    if settings['ligandFile']:
        if settings['ligandFile'].startswith('inputs/'):
            setting_name = os.path.basename(settings['ligandFile'])
            ligand_path = find_input_file(setting_name, ['.sdf', '.mol2', '.pdb'])
        else:
            ligand_path = settings['ligandFile']
    else:
        # Try to find ligandFile in inputs/
        ligand_path = find_input_file('ligandFile', ['.sdf', '.mol2', '.pdb'])
    
    return pdb_path, ligand_path

def run_simulation(pdb_path, ligand_path, settings):
    """Run OpenMM MD simulation with KCX support."""
    
    # Import OpenMM and related packages
    try:
        from openmm import app, unit, LangevinMiddleIntegrator, Platform
        from openmm.app import PDBFile, ForceField, Modeller, Simulation, PDBReporter, StateDataReporter, DCDReporter
        import openmm
    except ImportError:
        import simtk.openmm as openmm
        from simtk.openmm import app, unit, LangevinMiddleIntegrator, Platform
        from simtk.openmm.app import PDBFile, ForceField, Modeller, Simulation, PDBReporter, StateDataReporter, DCDReporter
    
    import parmed
    import mdtraj
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    
    # Create output directory
    os.makedirs('out', exist_ok=True)
    
    print(f"Loading PDB file: {pdb_path}")
    
    # Load structure with parmed for KCX support
    structure = parmed.load_file(pdb_path)
    
    # Check for KCX residues
    kcx_residues = [res for res in structure.residues if res.name == 'KCX']
    print(f"Found {len(kcx_residues)} KCX residues")
    
    # Check for zinc ions
    zn_atoms = [atom for atom in structure.atoms if atom.name == 'ZN' or atom.element_name == 'Zn']
    print(f"Found {len(zn_atoms)} zinc atoms")
    
    # Load force field with KCX parameters
    print("Setting up force field with KCX parameters...")
    
    # Use ff14SB for protein
    ff_files = ['amber14/protein.ff14SB.xml', 'amber14/tip3p.xml']
    
    # Check for KCX force field files
    kcx_frcmod = '/app/kcx.frcmod'
    kcx_lib = '/app/kcx.lib'
    
    if os.path.exists(kcx_frcmod) and os.path.exists(kcx_lib):
        print("Loading KCX parameters from frcmod/lib files...")
        # Load KCX parameters into parmed structure
        try:
            from parmed.amber import AmberParameterSet
            params = AmberParameterSet(kcx_frcmod)
            # Apply parameters
            structure.load_parameters(params)
        except Exception as e:
            print(f"Warning: Could not load KCX parameters: {e}")
    
    # Convert to OpenMM system
    print("Creating OpenMM system...")
    
    # Use parmed to create the system with proper parameters
    system = structure.createSystem(
        nonbondedMethod=app.PME,
        nonbondedCutoff=1.0*unit.nanometer,
        constraints=app.HBonds,
        rigidWater=True
    )
    
    # Set up integrator
    temperature = settings['temperature'] * unit.kelvin
    friction = settings['friction'] / unit.picosecond
    timestep = settings['timestep'] * unit.femtosecond
    
    integrator = LangevinMiddleIntegrator(temperature, friction, timestep)
    
    # Try to use CUDA, fall back to CPU
    try:
        platform = Platform.getPlatformByName('CUDA')
        properties = {'DeviceIndex': '0', 'Precision': 'mixed'}
        print("Using CUDA platform")
    except Exception:
        try:
            platform = Platform.getPlatformByName('OpenCL')
            properties = {}
            print("Using OpenCL platform")
        except Exception:
            platform = Platform.getPlatformByName('CPU')
            properties = {}
            print("Using CPU platform")
    
    # Create simulation
    simulation = Simulation(structure.topology, system, integrator, platform, properties if properties else {})
    simulation.context.setPositions(structure.positions)
    
    # Minimize energy
    print("Minimizing energy...")
    simulation.minimizeEnergy(maxIterations=1000)
    
    # Save minimized structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open('out/minimized.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    print("Saved minimized structure to out/minimized.pdb")
    
    # Set up reporters
    report_interval = settings['reportInterval']
    production_steps = int(settings['productionTime'] * 1e6 / settings['timestep'])  # ns to steps
    
    simulation.reporters.append(DCDReporter('out/trajectory.dcd', report_interval))
    simulation.reporters.append(StateDataReporter(
        'out/simulation_log.csv',
        report_interval,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True
    ))
    simulation.reporters.append(StateDataReporter(
        sys.stdout,
        report_interval * 10,
        step=True,
        time=True,
        potentialEnergy=True,
        temperature=True,
        speed=True
    ))
    
    # Run equilibration
    print(f"Running equilibration (100 ps)...")
    equilibration_steps = int(100000 / settings['timestep'])
    simulation.step(equilibration_steps)
    
    # Run production
    print(f"Running production ({settings['productionTime']} ns, {production_steps} steps)...")
    simulation.step(production_steps)
    
    # Save final structure
    positions = simulation.context.getState(getPositions=True).getPositions()
    with open('out/final.pdb', 'w') as f:
        PDBFile.writeFile(simulation.topology, positions, f)
    print("Saved final structure to out/final.pdb")
    
    # Analyze trajectory
    print("Analyzing trajectory...")
    analyze_trajectory(pdb_path)
    
    print("Simulation complete!")
    return True

def analyze_trajectory(pdb_path):
    """Analyze the MD trajectory and generate plots."""
    import mdtraj
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import numpy as np
    import pandas as pd
    
    # Load trajectory
    try:
        traj = mdtraj.load('out/trajectory.dcd', top=pdb_path)
        
        # Calculate RMSD
        rmsd = mdtraj.rmsd(traj, traj, 0) * 10  # Convert to Angstroms
        
        # Calculate RMSF
        rmsf = mdtraj.rmsf(traj, traj, 0) * 10  # Convert to Angstroms
        
        # Save RMSD data
        np.savetxt('out/rmsd.csv', rmsd, delimiter=',', header='RMSD (Angstrom)')
        
        # Save RMSF data
        np.savetxt('out/rmsf.csv', rmsf, delimiter=',', header='RMSF (Angstrom)')
        
        # Plot RMSD
        plt.figure(figsize=(10, 6))
        time_ns = np.arange(len(rmsd)) * traj.timestep / 1000  # Convert to ns
        plt.plot(time_ns, rmsd)
        plt.xlabel('Time (ns)')
        plt.ylabel('RMSD (Å)')
        plt.title('Backbone RMSD')
        plt.savefig('out/rmsd_plot.png', dpi=150)
        plt.close()
        
        # Plot RMSF
        plt.figure(figsize=(12, 6))
        plt.plot(rmsf)
        plt.xlabel('Residue Index')
        plt.ylabel('RMSF (Å)')
        plt.title('Per-Residue RMSF')
        plt.savefig('out/rmsf_plot.png', dpi=150)
        plt.close()
        
        print("Analysis complete. Generated RMSD and RMSF plots.")
        
    except Exception as e:
        print(f"Warning: Could not analyze trajectory: {e}")
    
    # Parse simulation log
    try:
        df = pd.read_csv('out/simulation_log.csv')
        
        # Plot energy
        plt.figure(figsize=(10, 6))
        plt.plot(df['#"Step"'], df['Potential Energy (kJ/mole)'])
        plt.xlabel('Step')
        plt.ylabel('Potential Energy (kJ/mol)')
        plt.title('Potential Energy vs Time')
        plt.savefig('out/energy_plot.png', dpi=150)
        plt.close()
        
        # Plot temperature
        plt.figure(figsize=(10, 6))
        plt.plot(df['#"Step"'], df['Temperature (K)'])
        plt.xlabel('Step')
        plt.ylabel('Temperature (K)')
        plt.title('Temperature vs Time')
        plt.savefig('out/temperature_plot.png', dpi=150)
        plt.close()
        
    except Exception as e:
        print(f"Warning: Could not parse simulation log: {e}")

def main():
    """Main entry point."""
    print("=" * 60)
    print("OpenMM MD Simulation with KCX Support")
    print("=" * 60)
    
    # Get settings
    settings = get_settings()
    print(f"Settings: {settings}")
    
    # List inputs directory contents
    if os.path.exists('inputs'):
        print(f"Files in inputs/: {os.listdir('inputs')}")
    
    # Prepare input files
    pdb_path, ligand_path = prepare_input_files(settings)
    
    print(f"PDB file: {pdb_path}")
    print(f"Ligand file: {ligand_path}")
    
    # Run simulation
    run_simulation(pdb_path, ligand_path, settings)

if __name__ == '__main__':
    main()
