#!/home/bbc8731/miniconda3/envs/modulediscovery/bin/python
# -*- coding: utf-8 -*-

import subprocess
import argparse
import os
import sys
import glob

def find_nextflow():
    nextflow_path = "/home/bbc8731/miniconda3/envs/modulediscovery/bin/nextflow"
    if os.path.exists(nextflow_path):
        return nextflow_path
    print("Error: Nextflow not found. Please ensure it's installed and in your PATH.")
    sys.exit(1)


def run_nextflow_pipeline(input_files, output_dir, network_file=None):
    main_nf_path = "/home/bbc8731/diseasemodulediscovery/main.nf"
    nextflow_exec = find_nextflow()

    # Set JAVA_HOME environment variable
    env = os.environ.copy()

    for input_file in input_files:
        # Get the base name of the input file without extension
        seed_name = os.path.splitext(os.path.basename(input_file))[0]
        # Create a subdirectory in the output directory for this seed
        seed_output_dir = os.path.join(output_dir, seed_name)
        os.makedirs(seed_output_dir, exist_ok=True)

        command = [
            nextflow_exec, 'run', main_nf_path,
            '-profile', 'singularity',
            '--input', input_file,
            '--outdir', seed_output_dir,
            '--id_space', 'uniprot'
            
        ]

        if network_file:
            command.extend(['--network', network_file])

        try:
            print(f"Running pipeline for seed file: {input_file}")
            subprocess.run(command, check=True, env=env)
            print(f"Pipeline completed successfully for {input_file}!")
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while running the pipeline for {input_file}: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run the modulediscovery Nextflow pipeline for multiple seed files")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('--input_dir', help="Path to the directory containing input seed files")
    group.add_argument('--input_files', nargs='+', help="List of input seed files")
    parser.add_argument('--pattern', default='*', help="Pattern to match seed files in the input directory (default: '*')")
    parser.add_argument('--outdir', required=True, help="Path to the output directory")
    parser.add_argument('--network', help="Path to the network file (optional)")

    args = parser.parse_args()

    if args.input_dir:
        # Collect all seed files from the input directory matching the pattern
        input_files = glob.glob(os.path.join(args.input_dir, args.pattern))
        if not input_files:
            print(f"No input seed files found in {args.input_dir} matching pattern '{args.pattern}'")
            sys.exit(1)
    else:
        # Use the list of input files provided
        input_files = args.input_files

    run_nextflow_pipeline(input_files, args.outdir, args.network)
