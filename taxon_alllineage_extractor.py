#!/usr/bin/env python

from ete3 import NCBITaxa

# Initialize NCBITaxa instance
ncbi = NCBITaxa()

# Paths to the input and output files
input_file = 'TGTin_uniprot_taxid_uniq.txt'
output_file = 'TGTin_uniprot_alllineage.txt'

# Open the input file for reading and the output file for writing
with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    for line in infile:
        taxon_id = line.strip()
        if taxon_id.isdigit():
            taxon_id = int(taxon_id)
            try:
                # Get the full lineage (as taxon IDs)
                lineage = ncbi.get_lineage(taxon_id)
                # Get the scientific names for the lineage
                names = ncbi.get_taxid_translator(lineage)
                # Create a string of lineage names
                lineage_names = "; ".join([names[taxid] for taxid in lineage])
                # Write the taxon ID and its lineage to the output file
                outfile.write(f"{taxon_id}: {lineage_names}\n")
            except Exception as e:
                # If there's an error (e.g., taxon ID not found), write an error message
                outfile.write(f"{taxon_id}: Error retrieving lineage\n")
        else:
            outfile.write(f"{line}: Invalid taxon ID\n")

print("Lineage extraction completed. Results written to", output_file)
