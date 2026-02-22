## Description
This project contains commands and scripts that haven been generated and used by the de Crecy lab at University of Florida for codon usage analysis in the manuscript:

“tRNA Modification Landscapes in Streptococci: Shared Losses and Clade-Specific Adaptations”

Authors: Ho-Ching Tiffany Tsui, Chi-Kong Chan, Yifeng Yuan, Roba Elias, Jingjing Sun, Virginie Marchand, Marshall Jaroch, Guangxin Sun, Irfan Manzoor, Ana Kutchuashvili, Yuri Motorin, Grazyna Leszczynska, Kinda Seaton, Kelly C. Rice, Manal A. Swairjo, Malcolm E. Winkler, Peter C. Dedon and Valérie de Crécy-Lagard

## Dependencies
MASH v2.3 https://github.com/marbl/mash

seqkit v2.8.0 https://bioinf.shenwei.me/seqkit/

python v3.12 https://www.python.org/

## Usage
Preparation
1. Retrive CDSs sequences of S. pnuemoniae strains and other Streptococcus strains in fasta format from NCBI (GCF_xxxx.ffn and GCA_xxx.ffn files) and put them in the working directory, for example $dir=work_dir, where GCF_xxx and GCA_xxx are NCBI assembly IDs.

2. Prepare two txt file of NCBI assembly IDs for the S.pnuemoniae and other Streptococcus, for example, list_Spne.txt and list_other.txt.

3. calculate the distance between S. pnuemoniae strains and other Streptococcus genomes using MASH.

```bash
# -p threads
# -l Lines in each <input> specify paths to sequence files, one per line.
mash sketch -p 4 -o Spne_sketch -l list_Spne.txt

mash sketch -p 4 -o other_sketch -l list_other.txt

mash dist -p 4 Spne_sketch.msh other_sketch.msh > mash_result.txt)

# check mash result
# remove duplicate

sed 's#working_dir/##g' mash_result.txt | sort -k2,2n > mash_result_sort.txt
```

4. Clean the MASH result. Then, calculate the distance between genomes so the sum of distance is minimal.

```
grep -v -E '^\s*$' remove_Qpattern_gid.txt | grep -F -v -f - mash_result_sort.txt > mash_result_clean.txt

python select_minimal_value_per_gid1.py mash_result_clean.txt mash_res_minimum_by_gidSpne.txt

python select_pair_with_global_minimal_distance.py mash_result_clean.txt mash_res_global_minimum.txt
```

## Help and Issues
Please contact Yifeng Yuan at yuanyifeng@ufl.edu

## Authors
Yifeng Yuan, Ph.D. Valerie de Crecy-Lagard (Principal Investigator)

Version History
v1.0 --2024/02/26 --creation of the project.
