#!/usr/bin/env python
import sys
from Bio.PDB import PDBParser

# Load the PDB file.
parser = PDBParser()

def print_pLDDT(id) :
    structure = parser.get_structure("protein", id)
    plddt_scores = []
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if atom.name == "CA":  # AlphaFold pLDDT scores are typically for C-alpha atoms.
                        plddt_scores.append((residue.id[1], atom.bfactor))
    # Print residue ID and pLDDT score
    for residue_id, score in plddt_scores:
        print(f"Res{residue_id}: {score}")


#uid = sys.argv[1]
#pdbfile = "AF-" + uid + "-F1-model_v4.pdb"
pdbfile = sys.argv[1]
print_pLDDT(pdbfile)
