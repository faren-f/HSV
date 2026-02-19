#!/home/bbc8731/miniconda3/envs/modulediscovery/bin/python

import os
import subprocess
import time
from pathlib import Path
import sys

#./run_sectioned_seeds.py /home/bbc8731/HSV1-2_pipeline/2_seed_gene_refinement/mtopf_intersection


def setup_environment():
    """Set up the required environment variables"""
    env_vars = {
        'TMPDIR': '/home/bbc8731/tmp',
        'NXF_TEMP': '/home/bbc8731/tmp',
        'NXF_WORK': '/home/bbc8731/nextflow_work',
        'SINGULARITY_TMPDIR': '/home/bbc8731/singularity/tmp',
        'SINGULARITY_CACHEDIR': '/home/bbc8731/singularity/cache',
        'SINGULARITY_BUILD_TMPDIR': '/home/bbc8731/singularity/tmp',
        'SINGULARITY_SESSIONDIR': '/home/bbc8731/singularity/tmp',
        'TMP': '/home/bbc8731/singularity/tmp',
        'TEMP': '/home/bbc8731/singularity/tmp'
    }
    
    conda_env = "/home/bbc8731/miniconda3/envs/modulediscovery"
    os.environ['PATH'] = f"{conda_env}/bin:{os.environ['PATH']}"
    
    for key, value in env_vars.items():
        os.environ[key] = value

def create_screen_command(seed_file, output_dir):
    """Create the screen and nextflow command"""
    nextflow_script = "/home/bbc8731/diseasemodulediscovery/main.nf"
    
    nextflow_cmd = (
        f"cd /home/bbc8731/diseasemodulediscovery && "
        f"nextflow run {nextflow_script} "
        f"-profile singularity "
        f"--seeds {seed_file} "
        f"--network /home/bbc8731/diseasemodulediscovery/tests/uniprot_ppi.csv "
        f"--outdir {output_dir} "
        f"--skip_annotation "
        f"--skip_evaluation "
        f"-executor.cpus=20 -executor.memory=100GB "
        f"-c fernando_nextflow.config -resume"
    )
    
    screen_name = f"nextflow_{Path(seed_file).stem}"
    # screen_cmd = f"screen -dmS {screen_name} bash -c 'source /home/bbc8731/miniconda3/etc/profile.d/conda.sh && conda activate modulediscovery && {nextflow_cmd}'"


    screen_cmd = (
        f"screen -dmS {screen_name} bash -lc "
        f"'source /home/bbc8731/miniconda3/etc/profile.d/conda.sh && "
        f"conda activate modulediscovery && "
        f"{nextflow_cmd} "
        f"> log.log 2>&1'"
    )




    return screen_cmd, screen_name

def is_screen_running(screen_name):
    """Check if a screen session is still running"""
    result = subprocess.run(['screen', '-ls'], capture_output=True, text=True)
    return screen_name in result.stdout

def wait_for_screen(screen_name):
    """Wait for a screen session to complete"""
    while is_screen_running(screen_name):
        time.sleep(30)  # Check every 30 seconds
        sys.stdout.write(f"\rWaiting for {screen_name} to complete...")
        sys.stdout.flush()
    print(f"\n{screen_name} has completed!")

def process_directory(base_dir):
    """Process all CSV files in the BP and CC directories"""
    bp_dir = os.path.join(base_dir, 'BP')
    cc_dir = os.path.join(base_dir, 'CC')
    output_base = "/home/bbc8731/HSV1-2_pipeline/3_module_expansion/results_full_PPI"
    
    if not os.path.exists(bp_dir):
        print(f"Warning: BP directory does not exist: {bp_dir}")
    if not os.path.exists(cc_dir):
        print(f"Warning: CC directory does not exist: {cc_dir}")
    
    os.makedirs(output_base, exist_ok=True)
    
    all_files = []
    
    # Process BP directory
    if os.path.exists(bp_dir):
        for csv_file in sorted(Path(bp_dir).glob('*.csv')):
            output_dir = os.path.join(output_base, f"BP_{csv_file.stem}")
            all_files.append((str(csv_file), output_dir))
    
    # Process CC directory
    if os.path.exists(cc_dir):
        for csv_file in sorted(Path(cc_dir).glob('*.csv')):
            output_dir = os.path.join(output_base, f"CC_{csv_file.stem}")
            all_files.append((str(csv_file), output_dir))
    
    return all_files

def main():
    if len(sys.argv) != 2:
        print("Usage: ./run_sectioned_seeds.py /path/to/mtopf_intersection")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    setup_environment()
    files_to_process = process_directory(base_dir)
    
    total_files = len(files_to_process)
    if total_files == 0:
        print("No files found to process!")
        sys.exit(1)
        
    print(f"\nFound {total_files} files to process")
    print(f"All outputs will be saved under: /home/bbc8731/HSV1-2_pipeline/3_module_expansion\n")
    
    for i, (seed_file, output_dir) in enumerate(files_to_process, 1):
        file_name = Path(seed_file).stem
        print(f"\nProcessing file {i}/{total_files}: {file_name}")
        print(f"Output directory: {output_dir}")
        
        screen_cmd, screen_name = create_screen_command(seed_file, output_dir)
        
        print(f"Starting screen session: {screen_name}")
        subprocess.run(screen_cmd, shell=True)
        print(f"To monitor: screen -r {screen_name}")
        
        # Wait for this job to complete before starting the next one
        wait_for_screen(screen_name)
        print(f"Completed processing {file_name}")
        print("-" * 80)

    print("\nAll pipeline runs have completed!")

if __name__ == "__main__":
    main()