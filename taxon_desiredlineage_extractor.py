#!/usr/bin/env python

from ete3 import NCBITaxa

# Initialize the NCBITaxa instance
ncbi = NCBITaxa()

# Define the ranks you're interested in
desired_ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']

# Read taxon IDs from input file
input_file = "TGTin_uniprot_taxid_uniq.txt"
output_file = "TGTin_uniprot_taxlineage.txt"

with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
    # Write header to the output file
    outfile.write("\t".join(['taxon_id'] + desired_ranks) + "\n")

    # Process each taxon ID
    for line in infile:
        taxon_id = line.strip()
        if taxon_id:
            try:
                # Get the full lineage for the taxon ID
                lineage = ncbi.get_lineage(taxon_id)

                # Get the names and ranks of the lineage
                names = ncbi.get_taxid_translator(lineage)
                ranks = ncbi.get_rank(lineage)

                # Prepare a dictionary to store the ranks
                lineage_dict = {ranks[taxid]: names[taxid] for taxid in lineage}

                # Extract the desired ranks
                output_line = [taxon_id]
                for rank in desired_ranks:
                    output_line.append(lineage_dict.get(rank, 'NA'))  # 'NA' if rank is missing

                # Write the result to the output file
                outfile.write("\t".join(output_line) + "\n")

            except Exception as e:
                print(f"Error processing taxon ID {taxon_id}: {e}")

print(f"Results written to {output_file}")
