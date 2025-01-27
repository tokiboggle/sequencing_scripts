#!/usr/bin/env python3
import os
import sys
import argparse
import subprocess
import venv
import pkg_resources
import shutil

def create_virtual_environment(project_dir):
    """Create virtual environment"""
    venv_path = os.path.join(project_dir, 'venv')
    if not os.path.exists(venv_path):
        print(f"Creating virtual environment in {venv_path}")
        venv.create(venv_path, with_pip=True)
    return os.path.join(venv_path, 'bin', 'python')

def install_dependencies(venv_python):
    """Install required dependencies"""
    subprocess.run([venv_python, '-m', 'pip', 'install', '-r', 'requirements.txt'], check=True)

def setup_project_files(project_dir, input_dir):
    """Setup project files including codon library"""
    # Create analysis results directory
    results_dir = os.path.join(input_dir, 'analysis_results')
    os.makedirs(results_dir, exist_ok=True)

    # Ensure codon library is in the correct location
    codon_lib_source = os.path.join(project_dir, 'codon_lib.csv')
    codon_lib_dest = os.path.join(os.path.dirname(input_dir), 'codon_lib.csv')
    
    if not os.path.exists(codon_lib_dest):
        if os.path.exists(codon_lib_source):
            shutil.copy2(codon_lib_source, codon_lib_dest)
        else:
            raise FileNotFoundError("Could not find codon_lib.csv in project directory")
    
    return results_dir

def run_analysis(input_dir, reference_dna, wild_type_protein):
    """Run analysis on sequencing files"""
    from src_analysis import process_files
    
    # Get project directory
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    # Setup results directory and ensure codon library is available
    results_dir = setup_project_files(project_dir, input_dir)
    
    # Run the analysis
    process_files(input_dir, results_dir, reference_dna, wild_type_protein)

def main():
    """Main entry point for the analysis pipeline"""
    parser = argparse.ArgumentParser(description='P450 Mutagenesis Analysis')
    parser.add_argument('input_dir', help='Directory with sequencing files')
    parser.add_argument('reference_dna', help='Reference DNA sequence')
    parser.add_argument('wild_type_protein', help='Wild-type protein sequence')
    
    args = parser.parse_args()
    
    project_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    venv_python = create_virtual_environment(project_dir)
    
    install_dependencies(venv_python)
    
    run_analysis(
        args.input_dir, 
        args.reference_dna, 
        args.wild_type_protein
    )

if __name__ == '__main__':
    main()
