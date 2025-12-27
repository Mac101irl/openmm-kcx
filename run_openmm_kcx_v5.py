#!/usr/bin/env python
"""
OpenMM MD Simulation with KCX (Carboxylated Lysine) Support - v7
Tamarind Platform Compatible Version

Conventions:
- File inputs: Filenames in env vars or inputs.json, fetch from Tamarind API
- Job Settings: Access via os.getenv('settingName')
- Outputs: Save ALL results to out/ directory
- Job name: Available as os.getenv('JobName')

Author: Generated for BMS Hydantoinase Project
Version: 7.0 (Tamarind Compatible with API file fetching)
"""

import os
import sys
import json
import subprocess
import shutil
import glob
import base64
import urllib.request
import urllib.error
import urllib.parse
from pathlib import Path

# =============================================================================
# TAMARIND FILE HANDLING
# =============================================================================

def load_inputs_json():
    """Load inputs.json which contains file references from Tamarind"""
    inputs_json_path = 'inputs/inputs.json'
    if os.path.exists(inputs_json_path):
        with open(inputs_json_path, 'r') as f:
            data = json.load(f)
            print(f"Loaded inputs.json: {json.dumps(data, indent=2)}")
            return data
    return {}

def get_api_key():
    """
    Get Tamarind API key from various sources.
    Priority: env var > inputs.json setting > None
    """
    # Check environment variables
    for key_name in ['TAMARIND_API_KEY', 'API_KEY', 'apiKey', 'x-api-key']:
        key = os.environ.get(key_name)
        if key:
            print(f"Found API key in env var: {key_name}")
            return key
    
    # Check if passed as a job setting in inputs.json
    inputs_data = load_inputs_json()
    api_key = inputs_data.get('apiKey') or inputs_data.get('API_KEY')
    if api_key:
        print("Found API key in inputs.json")
        return api_key
    
    # Check regular env var from job settings
    api_key = os.getenv('apiKey')
    if api_key:
        print("Found API key in apiKey env var")
        return api_key
    
    return None

def fetch_file_from_tamarind(filename, output_path, api_key=None):
    """
    Fetch a file from Tamarind's file storage.
    
    Uses the Tamarind API: https://app.tamarind.bio/api/files/{filename}
    """
    if api_key is None:
        api_key = get_api_key()
    
    if not api_key:
        print(f"  No API key available to fetch {filename}")
        print("  Tip: Pass 'apiKey' in job settings to enable file fetching")
        return False
    
    base_url = 'https://app.tamarind.bio/api'
    
    # Try to download the file
    try:
        # URL encode the filename for the request
        encoded_filename = urllib.parse.quote(filename, safe='')
        url = f"{base_url}/files/{encoded_filename}"
        print(f"  Attempting to fetch: {url}")
        
        req = urllib.request.Request(url)
        req.add_header('x-api-key', api_key)
        
        with urllib.request.urlopen(req, timeout=60) as response:
            content = response.read()
            
            os.makedirs(os.path.dirname(output_path) or '.', exist_ok=True)
            with open(output_path, 'wb') as f:
                f.write(content)
            
            print(f"  Successfully downloaded {filename} to {output_path}")
            return True
            
    except urllib.error.HTTPError as e:
        print(f"  HTTP Error fetching {filename}: {e.code} {e.reason}")
    except urllib.error.URLError as e:
        print(f"  URL Error fetching {filename}: {e.reason}")
    except Exception as e:
        print(f"  Error fetching {filename}: {e}")
    
    return False

def search_for_file_globally(filename, extensions=None):
    """
    Search entire filesystem for a file by name.
    This is a last resort to find where Tamarind might place files.
    """
    if extensions is None:
        extensions = []
    
    search_dirs = [
        '/', '/app', '/home', '/data', '/input', '/inputs', '/files',
        '/var', '/tmp', '/mnt', '/opt', '/root', '/workspace',
        os.path.expanduser('~'), os.getcwd()
    ]
    
    # Get just the basename
    basename = os.path.basename(filename)
    name_without_ext = os.path.splitext(basename)[0]
    
    patterns_to_search = [basename, name_without_ext]
    for ext in extensions:
        patterns_to_search.append(f"{name_without_ext}{ext}")
    
    found_files = []
    
    for search_dir in search_dirs:
        if not os.path.exists(search_dir):
            continue
        try:
            for root, dirs, files in os.walk(search_dir):
                # Don't recurse too deep
                depth = root.replace(search_dir, '').count(os.sep)
                if depth > 4:
                    dirs.clear()
                    continue
                    
                for pattern in patterns_to_search:
                    if pattern in files:
                        found_path = os.path.join(root, pattern)
                        found_files.append(found_path)
                        print(f"  FOUND: {found_path}")
        except PermissionError:
            continue
        except Exception as e:
            continue
    
    return found_files


def find_input_file(setting_name, extensions=None):
    """
    Find an input file for Tamarind jobs.
    
    Strategy:
    1. Check if file exists directly in inputs/ directory
    2. Get filename from environment variable
    3. Get filename from inputs.json
    4. Search common directories
    5. Try to fetch from Tamarind API (requires API key in env)
    6. Global filesystem search
    
    Args:
        setting_name: The setting name (e.g., 'pdbFile')
        extensions: List of extensions to try (e.g., ['.pdb', '.cif'])
    
    Returns:
        Path to the found/downloaded file, or None if not found
    """
    if extensions is None:
        extensions = []
    
    print(f"\nSearching for {setting_name}...")
    
    # Strategy 1: Check if file exists directly in inputs/
    base_path = f'inputs/{setting_name}'
    if os.path.exists(base_path) and os.path.isfile(base_path):
        print(f"  Found at: {base_path}")
        return base_path
    
    for ext in extensions:
        path_with_ext = f'{base_path}{ext}'
        if os.path.exists(path_with_ext) and os.path.isfile(path_with_ext):
            print(f"  Found at: {path_with_ext}")
            return path_with_ext
    
    # Strategy 2: Get filename from environment variable
    filename = os.getenv(setting_name)
    print(f"  Env var {setting_name} = {filename}")
    
    # Strategy 3: Get filename from inputs.json
    if not filename:
        inputs_data = load_inputs_json()
        filename = inputs_data.get(setting_name)
        print(f"  inputs.json {setting_name} = {filename}")
    
    if not filename:
        print(f"  No filename found for {setting_name}")
        return None
    
    # Check if the filename itself is a path that exists
    if os.path.exists(filename):
        print(f"  File exists at: {filename}")
        return filename
    
    # Check in inputs/ directory
    inputs_path = f'inputs/{filename}'
    if os.path.exists(inputs_path):
        print(f"  File exists at: {inputs_path}")
        return inputs_path
    
    # Check in current directory
    if os.path.exists(filename):
        print(f"  File exists at: {filename}")
        return filename
    
    # Strategy 4: Check common directories
    common_dirs = ['/data', '/files', '/input', '/app/inputs', '/app/data', 
                   '/home', '/tmp', '/workspace', '/mnt']
    for dir_path in common_dirs:
        test_path = os.path.join(dir_path, filename)
        if os.path.exists(test_path):
            print(f"  File exists at: {test_path}")
            return test_path
        for ext in extensions:
            test_path_ext = os.path.join(dir_path, os.path.splitext(filename)[0] + ext)
            if os.path.exists(test_path_ext):
                print(f"  File exists at: {test_path_ext}")
                return test_path_ext
    
    # Strategy 5: Try to fetch from Tamarind API (only if API key available in env)
    output_path = f'inputs/{filename}'
    if fetch_file_from_tamarind(filename, output_path):
        return output_path
    
    # Strategy 6: Global filesystem search
    print(f"  Performing global filesystem search for {filename}...")
    found = search_for_file_globally(filename, extensions)
    if found:
        return found[0]  # Return first match
    
    print(f"  Could not find or fetch {setting_name}")
    return None

# =============================================================================
# SETTINGS
# =============================================================================

def get_env(name, default=None, type_fn=str):
    """Get environment variable with type conversion"""
    val = os.getenv(name)
    if val is None or val == '' or val == 'undefined':
        return default
    if type_fn == bool:
        return str(val).lower() in ('true', 'yes', '1')
    try:
        return type_fn(val)
    except (ValueError, TypeError):
        return default

def get_settings():
    """Get all settings from environment variables and inputs.json"""
    
    # Load inputs.json for any settings not in env vars
    inputs_data = load_inputs_json()
    
    def get_setting(name, default=None, type_fn=str):
        """Get setting from env var or inputs.json"""
        val = os.getenv(name)
        if val is None or val == '' or val == 'undefined':
            val = inputs_data.get(name)
        if val is None:
            return default
        if type_fn == bool:
            return str(val).lower() in ('true', 'yes', '1')
        try:
            return type_fn(val)
        except (ValueError, TypeError):
            return default
    
    # Find input files
    pdb_file = find_input_file('pdbFile', ['.pdb', '.cif', '.ent'])
    ligand_file = find_input_file('ligandFile', ['.sdf', '.mol2', '.pdb'])
    
    settings = {
        # Job info
        'jobName': get_setting('JobName', 'openmm_simulation'),
        
        # Input files
        'pdbFile': pdb_file,
        'ligandFile': ligand_file,
        'ligandCharge': get_setting('ligandCharge', 0, int),
        
        # Force field
        'forceField': get_setting('forceField', 'ff19SB'),
        'waterModel': get_setting('waterModel', 'tip3p'),
        
        # Solvation
        'boxSize': get_setting('boxSize', 12.0, float),
        'boxShape': get_setting('boxShape', 'cube'),
        'ionicStrength': get_setting('ionicStrength', 0.15, float),
        'positiveIon': get_setting('positiveIon', 'Na+'),
        'negativeIon': get_setting('negativeIon', 'Cl-'),
        'pH': get_setting('pH', 7.0, float),
        'removeWaters': get_setting('removeWaters', 'no'),
        
        # Nonbonded
        'nonbondedMethod': get_setting('nonbondedMethod', 'PME'),
        'nonbondedCutoff': get_setting('nonbondedCutoff', 1.0, float),
        'ewaldErrorTolerance': get_setting('ewaldErrorTolerance', 0.0005, float),
        'switchDistance': get_setting('switchDistance', None, float),
        
        # Constraints
        'constraints': get_setting('constraints', 'HBonds'),
        'rigidWater': get_setting('rigidWater', 'yes'),
        'hydrogenMass': get_setting('hydrogenMass', None, float),
        
        # Simulation
        'minimizationSteps': get_setting('minimizationSteps', 10000, int),
        'minimizationTolerance': get_setting('minimizationTolerance', 10.0, float),
        'equilibrationTime': get_setting('equilibrationTime', 0.2, float),
        'productionTime': get_setting('productionTime', 1.0, float),
        'timestep': get_setting('timestep', 2.0, float),
        
        # Temperature & Pressure
        'temperature': get_setting('temperature', 310.0, float),
        'pressure': get_setting('pressure', 1.0, float),
        'frictionCoeff': get_setting('frictionCoeff', 1.0, float),
        'barostatInterval': get_setting('barostatInterval', 25, int),
        
        # Integrator
        'integrator': get_setting('integrator', 'LangevinMiddle'),
        
        # Output
        'equilTrajFreq': get_setting('equilTrajFreq', 1000, int),
        'prodTrajFreq': get_setting('prodTrajFreq', 1000, int),
        'checkpointFreq': get_setting('checkpointFreq', 10000, int),
        'trajectoryFormat': get_setting('trajectoryFormat', 'DCD'),
        
        # Analysis
        'stepSize': get_setting('stepSize', 5, int),
        'rmsdMask': get_setting('rmsdMask', 'backbone'),
    }
    
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

def prepare_ligand(ligand_file, ligand_charge, output_dir='prep'):
    """Prepare ligand with antechamber using GAFF2 and AM1-BCC charges"""
    print(f"\n{'='*60}")
    print(f"LIGAND PREPARATION")
    print(f"{'='*60}")
    print(f"  Ligand file: {ligand_file}")
    print(f"  Net charge: {ligand_charge}")
    
    os.makedirs(output_dir, exist_ok=True)
    
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
    
    print("  Running antechamber...")
    cmd = [
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

def prepare_protein(pdb_file, remove_waters, output_dir='prep'):
    """Prepare protein with pdb4amber"""
    print(f"\n{'='*60}")
    print(f"PROTEIN PREPARATION")
    print(f"{'='*60}")
    print(f"  PDB file: {pdb_file}")
    
    os.makedirs(output_dir, exist_ok=True)
    
    with open(pdb_file, 'r') as f:
        pdb_content = f.read()
    
    has_kcx = 'KCX' in pdb_content
    has_zn = ' ZN ' in pdb_content or 'ZN2' in pdb_content
    
    print(f"  Contains KCX: {has_kcx}")
    print(f"  Contains Zn: {has_zn}")
    
    protein_out = os.path.join(output_dir, 'protein.pdb')
    
    cmd = ['pdb4amber', '-i', pdb_file, '-o', protein_out]
    if remove_waters == 'yes':
        cmd.append('--dry')
    
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"  WARNING: pdb4amber warning: {result.stderr}")
        shutil.copy(pdb_file, protein_out)
    
    print(f"  Protein prepared: {protein_out}")
    return protein_out, has_kcx, has_zn


# =============================================================================
# TLEAP SYSTEM BUILDING
# =============================================================================

def create_tleap_input(settings, protein_pdb, ligand_mol2, ligand_frcmod, 
                       kcx_frcmod, kcx_lib, has_kcx, has_zn):
    """Create tleap input script"""
    
    ff_protein = 'leaprc.protein.ff19SB' if settings.forceField == 'ff19SB' else 'leaprc.protein.ff14SB'
    
    water_models = {
        'tip3p': 'leaprc.water.tip3p',
        'opc': 'leaprc.water.opc',
        'spce': 'leaprc.water.spce',
    }
    ff_water = water_models.get(settings.waterModel, 'leaprc.water.tip3p')
    
    n_ions = max(1, int(settings.ionicStrength * 60))
    
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
    
    if settings.boxShape == 'octahedron':
        script += f"solvateOct complex TIP3PBOX {settings.boxSize}\n"
    else:
        script += f"solvatebox complex TIP3PBOX {settings.boxSize}\n"
    
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

def run_simulation(settings):
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
    constraints = {'HBonds': app.HBonds, 'AllBonds': app.AllBonds, 'None': None}
    
    system = prmtop.createSystem(
        nonbondedMethod=nb_methods.get(settings.nonbondedMethod, app.PME),
        nonbondedCutoff=settings.nonbondedCutoff * unit.nanometer,
        constraints=constraints.get(settings.constraints, app.HBonds),
        rigidWater=settings.rigidWater == 'yes'
    )
    
    system.addForce(mm.MonteCarloBarostat(
        settings.pressure * unit.bar,
        settings.temperature * unit.kelvin,
        settings.barostatInterval
    ))
    
    print(f"Creating {settings.integrator} integrator...")
    timestep = settings.timestep * unit.femtosecond
    temperature = settings.temperature * unit.kelvin
    friction = settings.frictionCoeff / unit.picosecond
    
    integrator = mm.LangevinMiddleIntegrator(temperature, friction, timestep)
    
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
    
    simulation = app.Simulation(prmtop.topology, system, integrator, platform, properties)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)
    
    # Minimization
    print(f"\nMinimizing ({settings.minimizationSteps} steps)...")
    simulation.minimizeEnergy(maxIterations=settings.minimizationSteps)
    
    state = simulation.context.getState(getPositions=True)
    with open('out/minimized.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    # Equilibration
    equil_steps = int(settings.equilibrationTime * 1e6 / settings.timestep)
    print(f"\nEquilibrating ({settings.equilibrationTime} ns)...")
    
    simulation.reporters.append(app.StateDataReporter(
        'out/equilibration.log', settings.equilTrajFreq,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True
    ))
    simulation.step(equil_steps)
    simulation.reporters.clear()
    
    # Production
    prod_steps = int(settings.productionTime * 1e6 / settings.timestep)
    print(f"\nProduction ({settings.productionTime} ns)...")
    
    simulation.reporters.append(app.StateDataReporter(
        'out/production.log', settings.prodTrajFreq,
        step=True, time=True, potentialEnergy=True, temperature=True, speed=True
    ))
    simulation.reporters.append(app.DCDReporter('out/production.dcd', settings.prodTrajFreq))
    simulation.reporters.append(app.CheckpointReporter('out/checkpoint.chk', settings.checkpointFreq))
    
    simulation.step(prod_steps)
    
    # Save final
    state = simulation.context.getState(getPositions=True)
    with open('out/final.pdb', 'w') as f:
        app.PDBFile.writeFile(simulation.topology, state.getPositions(), f)
    
    print("\nSimulation complete!")


# =============================================================================
# ANALYSIS
# =============================================================================

def run_analysis(settings):
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
    traj = md.load('out/production.dcd', top='system.pdb', stride=settings.stepSize)
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
        plt.savefig('out/rmsf.png', dpi=150)
        plt.close()
    
    print("  Analysis complete!")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*70)
    print("  OpenMM MD Simulation with KCX Support - v7")
    print("="*70)
    
    # Debug output
    print("\nEnvironment variables (file-related):")
    for key in ['pdbFile', 'ligandFile', 'ligandCharge', 'productionTime', 
                'JobName', 'TAMARIND_API_KEY', 'API_KEY']:
        print(f"  {key} = {os.getenv(key)}")
    
    print("\nFiles in inputs/:")
    if os.path.exists('inputs'):
        for f in os.listdir('inputs'):
            print(f"  {f}")
    
    if os.path.exists('inputs/inputs.json'):
        print("\ninputs.json content:")
        with open('inputs/inputs.json') as f:
            print(f.read())
    
    # Get settings
    settings = Settings(get_settings())
    
    print("\n" + "="*70)
    print("SIMULATION PARAMETERS:")
    print("-"*70)
    print(f"  PDB file:       {settings.pdbFile}")
    print(f"  Ligand file:    {settings.ligandFile}")
    print(f"  Production:     {settings.productionTime} ns")
    print(f"  Temperature:    {settings.temperature} K")
    print("="*70)
    
    # Check PDB file
    if settings.pdbFile is None or not os.path.exists(settings.pdbFile):
        print("\nERROR: PDB file not found!")
        print("The file needs to be fetched from Tamarind's file storage.")
        print("Make sure TAMARIND_API_KEY is set in the container environment.")
        sys.exit(1)
    
    # Setup and run
    kcx_frcmod, kcx_lib = write_kcx_parameters('params')
    
    ligand_mol2, ligand_frcmod = None, None
    if settings.ligandFile and os.path.exists(settings.ligandFile):
        ligand_mol2, ligand_frcmod = prepare_ligand(
            settings.ligandFile, settings.ligandCharge, 'prep'
        )
    
    protein_pdb, has_kcx, has_zn = prepare_protein(
        settings.pdbFile, settings.removeWaters, 'prep'
    )
    
    tleap_file = create_tleap_input(
        settings, protein_pdb, ligand_mol2, ligand_frcmod,
        kcx_frcmod, kcx_lib, has_kcx, has_zn
    )
    
    run_tleap(tleap_file)
    run_simulation(settings)
    run_analysis(settings)
    
    print("\n" + "="*70)
    print("  ALL DONE!")
    print("="*70)


if __name__ == '__main__':
    main()
    
