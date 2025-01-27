import os
import json
import glob
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqUtils import seq1
import numpy as np
import csv

def load_codon_table():
    """Load codon table from CSV file"""
    codon_table = {}
    with open('codon_lib.csv', 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            codon_table[row['Triplet']] = row['Amino acid']
    return codon_table

def translate_dna(sequence, codon_table):
    """Translate DNA sequence using custom codon table"""
    protein = []
    for i in range(0, len(sequence), 3):
        codon = str(sequence[i:i+3])
        if len(codon) == 3:  # Only process complete codons
            amino_acid = codon_table.get(codon, 'X')  # 'X' for unknown/invalid codons
            protein.append(amino_acid)
    return ''.join(protein)

def parse_sequences(input_file):
    """Parse sequences from FASTQ or PHD file"""
    file_ext = input_file.lower().split('.')[-1]
    if file_ext in ['fastq', 'fq']:
        return list(SeqIO.parse(input_file, "fastq"))
    elif file_ext == 'phd':
        sequences = []
        with open(input_file, 'r') as f:
            content = f.read()
            sequences = [SeqIO.SeqRecord(Seq(seq)) for seq in content.split()]
        return sequences
    else:
        raise ValueError(f"Unsupported file type: {file_ext}")

def detect_mutations(reference_dna, sample_dna, codon_table):
    """Detect DNA and protein mutations using custom codon table"""
    dna_mutations = [
        {'position': pos+1, 'reference': ref, 'sample': sample} 
        for pos, (ref, sample) in enumerate(zip(reference_dna, sample_dna)) 
        if ref != sample
    ]
    
    # Use custom codon table for translation
    ref_protein = translate_dna(reference_dna, codon_table)
    sample_protein = translate_dna(sample_dna, codon_table)
    
    protein_mutations = [
        {'position': pos+1, 'wild_type': wt, 'mutant': sample}  
        for pos, (wt, sample) in enumerate(zip(ref_protein, sample_protein)) 
        if wt != sample
    ]
    
    return {
        'dna_mutations': dna_mutations,
        'protein_mutations': protein_mutations
    }

def calculate_quality_metrics(sequences):
    """Calculate sequencing quality metrics"""
    lengths = [len(seq) for seq in sequences]
    base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
    total_bases = 0
    
    for seq in sequences:
        for base in seq:
            if base.upper() in base_counts:
                base_counts[base.upper()] += 1
                total_bases += 1
    
    return {
        'total_sequences': len(sequences),
        'sequence_lengths': {
            'mean': np.mean(lengths),
            'min': min(lengths),
            'max': max(lengths)
        },
        'base_composition': {
            base: (count / total_bases * 100) 
            for base, count in base_counts.items()
        }
    }

def process_files(input_dir, results_dir, reference_dna, wild_type_protein):
    """Process all sequencing files in directory"""
    reference_dna = Seq(reference_dna)
    codon_table = load_codon_table()
    
    # Find input files
    input_files = (
        glob.glob(os.path.join(input_dir, '*.fastq')) + 
        glob.glob(os.path.join(input_dir, '*.fq')) + 
        glob.glob(os.path.join(input_dir, '*.phd'))
    )
    
    for input_file in input_files:
        # Parse sequences
        sequences = parse_sequences(input_file)
        
        # Prepare results
        results = {
            'input_file': input_file,
            'quality_metrics': calculate_quality_metrics(sequences),
            'mutations': []
        }
        
        # Analyze each sequence
        for seq in sequences:
            mutation_details = detect_mutations(reference_dna, seq.seq, codon_table)
            results['mutations'].append({
                'sequence_id': seq.id,
                'mutations': mutation_details
            })
        
        # Write results
        output_filename = os.path.join(
            results_dir, 
            f"{os.path.splitext(os.path.basename(input_file))[0]}_results.json"
        )
        with open(output_filename, 'w') as f:
            json.dump(results, f, indent=2)
        
        print(f"Analysis complete for {input_file}. Results in {output_filename}")
