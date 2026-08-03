#!/usr/bin/env python3
"""
Description: Python script to get the longest protein sequence for a list of 
NCBI GeneIDs, excluding pseudogenes, using NCBI Entrez.

- Input: TSV file with 2 columns: GeneIDs / Organism
- Output: FASTA file

Usage:	python3 -m pip install -r requirements.txt
		python3 get_biggest_geneids.py i- input_name.tsv -o output_name.fasta

Don't forget to add a valid NCBI Entrez email (Entrez.mail)

Author: Lucas Rocha
"""

import os
import time
import concurrent.futures
from Bio import Entrez, SeqIO
import pandas as pd
import threading
import argparse

##### ------------ Email config  ------------ #####

Entrez.email = "your@email.com"  # add your own

##### --------------------------------------- #####



# argparse config
parser = argparse.ArgumentParser(description="Fetches the longest non-pseudogene protein by GeneID from NCBI")
parser.add_argument("-i", "--input", help="Input TSV file", required=True)
parser.add_argument("-o", "--output", help="Output FASTA file", required=True)
args = parser.parse_args()

print_lock = threading.Lock()

def read_geneids_file(file_path):
    """Reads the TSV file containing GeneIDs and Org_names."""
    try:
        df = pd.read_csv(file_path, sep='\t', header=None)
        # Assuming the file has at least two columns: GeneID and Org_name
        if len(df.columns) >= 2:
            df.columns = ['GeneID', 'Org_name'] + [f'column_{i}' for i in range(2, len(df.columns))]
        else:
            df.columns = ['GeneID']
            df['Org_name'] = None
        return df
    except Exception as e:
        with print_lock:
            print(f"Error reading GeneID file: {e}")
        return None

def is_pseudogene(record):
    """
    Checks if a protein record is associated with a pseudogene.
    
    Args:
        record: A Biopython SeqRecord object with protein information
    
    Returns:
        bool: True if it is a pseudogene, False otherwise
    """
    # Keywords that may indicate a pseudogene
    keywords = ['pseudogene', 'pseudo', 'nonfunctional', 'non-functional']
    
    # Check in title/description
    if any(keyword in record.description.lower() for keyword in keywords):
        return True
        
    # Check in features
    for feature in record.features:
        # Check qualifiers like note, product, gene_desc
        for qualifier in ['note', 'product', 'gene_desc']:
            if qualifier in feature.qualifiers:
                for value in feature.qualifiers[qualifier]:
                    if any(keyword in value.lower() for keyword in keywords):
                        return True
                        
        # Check in specific pseudogene qualifier
        if 'pseudogene' in feature.qualifiers:
            return True
            
    # Check in other annotations
    for annotation_key, annotation_value in record.annotations.items():
        if isinstance(annotation_value, str):
            if any(keyword in annotation_value.lower() for keyword in keywords):
                return True
        elif isinstance(annotation_value, list):
            for item in annotation_value:
                if isinstance(item, str) and any(keyword in item.lower() for keyword in keywords):
                    return True
    
    return False

def process_gene(gene_info):
    """
    Processes a single gene to find the longest non-pseudogene protein.
    
    Args:
        gene_info: Tuple containing (GeneID, Org_name)
    
    Returns:
        tuple: (GeneID, protein, protein_id) or (GeneID, None, None) if not found
    """
    GeneID, Org_name = gene_info
    
    with print_lock:
        print(f"Processing gene {GeneID}...", end=" ", flush=True)
    
    try:
        # Build the query
        query = f"{GeneID}[Gene ID]"
        if Org_name:
            query += f" AND {Org_name}[Organism]"
        query += " NOT pseudogene[All Fields]"  # Attempt to exclude pseudogenes in query
        
        # First fetch protein IDs associated with the gene
        search_handle = Entrez.esearch(db="protein", term=query, retmax=100)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        protein_ids = search_results.get("IdList", [])
        
        if not protein_ids:
            with print_lock:
                print("No proteins found")
            return GeneID, None, None
        
        # Fetch information for all found proteins
        fetch_handle = Entrez.efetch(db="protein", id=",".join(protein_ids), rettype="gb", retmode="text")
        records = list(SeqIO.parse(fetch_handle, "genbank"))
        fetch_handle.close()
        
        # Filter out pseudogenes
        filtered_records = []
        for record in records:
            if not is_pseudogene(record):
                filtered_records.append(record)
            else:
                with print_lock:
                    print(f"  Ignored pseudogene: {record.id}")
        
        # Find the longest sequence among the filtered records
        if filtered_records:
            longest_protein = max(filtered_records, key=lambda rec: len(rec.seq))
            # Add GeneID to the sequence description
            longest_protein.description = f"GeneID:{GeneID} | {longest_protein.description}"
            with print_lock:
                print(f"OK - Protein: {longest_protein.id} ({len(longest_protein.seq)} aa)")
            return GeneID, longest_protein, longest_protein.id
        else:
            with print_lock:
                print("No valid proteins (all pseudogenes)")
            return GeneID, None, None
            
    except Exception as e:
        with print_lock:
            print(f"Error: {e}")
        return GeneID, None, None

def main():
    input_file = args.input
    output_file = args.output
    
    # Number of workers for parallel processing
    # Adjust this value based on your internet connection and CPU
    # A value between 4-10 is usually good to balance performance
    # and respect NCBI request limits
    num_workers = 5
    
    start_time = time.time()
    
    # Read input file
    df = read_geneids_file(input_file)
    if df is None:
        return
    
    print(f"Processing {len(df)} genes in parallel with {num_workers} workers...")
    
    # Create list of tuples (GeneID, Org_name) for processing
    genes_to_process = [(str(row['GeneID']), row.get('Org_name')) for _, row in df.iterrows()]
    
    # Dictionary to store results
    results = {}
    
    # Execute parallel processing
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        # Submit all tasks
        future_to_gene = {
            executor.submit(process_gene, gene_info): gene_info[0] 
            for gene_info in genes_to_process
        }
        
        # Process results as they are completed
        completed = 0
        for future in concurrent.futures.as_completed(future_to_gene):
            GeneID = future_to_gene[future]
            try:
                GeneID, protein, protein_id = future.result()
                results[GeneID] = protein
                
                # Update progress counter
                completed += 1
                if completed % 5 == 0 or completed == len(genes_to_process):
                    print(f"Progress: {completed}/{len(genes_to_process)} genes processed")
                
            except Exception as e:
                print(f"Error processing gene {GeneID}: {e}")
    
    # Filter only found proteins
    found_proteins = [protein for protein in results.values() if protein is not None]
    
    # Save results to FASTA file
    if found_proteins:
        with open(output_file, "w") as output_handle:
            SeqIO.write(found_proteins, output_handle, "fasta")
        
        end_time = time.time()
        duration = end_time - start_time
        
        print(f"\nProcessing completed in {duration:.2f} seconds.")
        print(f"{len(found_proteins)} proteins were saved to '{output_file}'.")
        print(f"Genes not found or containing only pseudogenes: {len(df) - len(found_proteins)}")
    else:
        print("\nNo proteins were found.")

if __name__ == "__main__":
    main()