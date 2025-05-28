# This file contains scripts and commands used in the manuscript below.
# The oncogene SLC35F2 is a high-specificity transporter for the micronutrients queuine and queuosine
# Authors: Lyubomyr Burtnyaka, Yifeng Yuanb, Xiaobei Panc, Lankani Gunaratned, Gabriel Silveira D’Almeidad, \
# Maria Martinellib&, Colbie Reedc$, Jessie Fernandez Garciab, Bhargesh Indravadan Patele, Isaac Marquezf, \
# Ann E. Ehrenhofer-Murraye, Manal A. Swairjof, Juan D. Alfonzod, Brian D. Greenc, Vincent P. Kellya* and Valérie de Crécy-Lagardb,g*

# Author information
# aSchool of Biochemistry & Immunology, Trinity Biomedical Sciences Institute, Trinity College Dublin, Dublin 2, Ireland
# bDepartment of Microbiology and Cell Science, University of Florida, Gainesville, Florida 32611, USA
# cSchool of Biological Sciences, Institute for Global Food Security, Queen's University Belfast, Belfast, UK
# dDepartment of Molecular Biology, Cell Biology and Biochemistry and the Brown RNA Center, Brown University, Providence, Rhode Island 02903, USA
# eInstitut für Biologie, Lebenswissenschaftliche Fakultät, Humboldt-Universität zu Berlin, 10115 Berlin, Germany
# fDepartment of Chemistry and Biochemistry, San Diego State University, San Diego, California 92182, USA
# gGenetic Institute, University of Florida, Gainesville, Florida 32611, USA

# &Current address: eSTEAMed Learning Inc., Maitland, FL 32751, USA
# $Current address: Data Science Institute, Medical College of Wisconsin, Milwaukee, Wisconsin 53226, USA

# *Corresponding authors: Valérie de Crécy-Lagard and Vincent P. Kelly
# Email: vcrecy@ufl.edu; kellyvp@tcd.ie

# scripts below performed on hpc.
# copy right and contact yuanyifeng@ufl.edu.

# I. select Q+ and Q- organisms.
# 1. download protein sequences.
cd /blue/lagard/yuanyifeng/Qtaxon

wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
wget https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz

gunzip uniprot_sprot.dat.gz
gunzip uniprot_trembl.dat.gz

# 2. search QTRT1 using diamond
# submit the script to SLURM.

#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=3-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=8
#SBATCH --mem=10gb
###############################################################
module load diamond/2.1.8

db=/blue/lagard/yuanyifeng/Qtaxon/uniprot_sprot.fasta.gz

diamond makedb --in ${db} -d uniprotsp --threads 12

diamond blastp --db uniprotsp -q query.fasta --threads 12 --outfmt 6 --very-sensitive --matrix BLOSUM45 --out TGTin_sp.out --evalue 0.001 -k0

db=/blue/lagard/yuanyifeng/Qtaxon/uniprot_trembl.fasta.gz

diamond makedb --in ${db} -d uniprottrembl --threads 12

diamond blastp --db uniprottrembl -q query.fasta --threads 12 --outfmt 6 --very-sensitive --matrix BLOSUM45 --out TGTin_trembl.out --evalue 0.001 -k0

#### END ####

# 3. count TGT per taxon id.
> TGTin_sp_taxid.txt

grep '^TGT' TGTin_sp.out | cut -d '|' -f2 | while read name ;  do
  echo -e "$name\t\c" >> TGTin_sp_taxid.txt
  LC_ALL=C grep -m1 -A 20 -e "DE .*$name" -e "AC .*$name" uniprot_sprot.dat | grep -m1 '^OX ' >> TGTin_sp_taxid.txt
done

> TGTin_trembl_taxid.txt
cut -d $'\t' -f2 TGTin_trembl.out | while read title ;  do
  LC_ALL=C grep -m1 "^>$title" uniprot_seq_title.txt |cut -d '|' -f2,3 |sed 's/ .*OX=/#/' | sed 's/ .*$//' | sed 's/|.*#/\tOX\tNCBI_TaxID=/'>> TGTin_trembl_taxid.txt
done

cat TGTin_sp_taxid.txt TGTin_trembl_taxid.txt >> TGTin_uniprot_taxid.txt

cut -d '=' -f2 TGTin_uniprot_taxid.txt | sort | uniq -c > TGTin_uniprot_taxid_count.txt

# 4. filter taxon id with >=5000 sequences.
sed 's/^.* OX=//' uniprot_seq_title.txt | sed 's/ .*$//' > uniprot_seq_OX.txt

sort uniprot_seq_OX.txt | uniq -c > uniprot_seq_OX_count.txt

sort -k1,1gr uniprot_seq_OX_count.txt > uniprot_seq_OX_count_sorted.txt

awk '$1>=5000 {print $0}' uniprot_seq_OX_count_sorted.txt > uniprot_seq_OX_count5000.txt

# get taxon lineage.
python3 taxon_desiredlineage_extractor.py

# group results by Eukaryota kingdoms/clades.
grep ' root; Viruses; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Viruses
grep ' root; unclassified entries; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.unclassified
grep ' root; cellular organisms; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.cellular
grep ' root; cellular organisms; Bacteria; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Bacteria
grep ' root; cellular organisms; Archaea; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Archaea
grep ' root; cellular organisms; Eukaryota; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Eukaryota
grep ' root; cellular organisms; Eukaryota; Opisthokonta; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Opisthokonta
grep ' root; cellular organisms; Eukaryota; Viridiplantae; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Viridiplantae
grep ' root; cellular organisms; Eukaryota; Sar; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Sar
grep ' root; cellular organisms; Eukaryota; Discoba; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Discoba
grep ' root; cellular organisms; Eukaryota; Amoebozoa; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Amoebozoa
grep ' root; cellular organisms; Eukaryota; Metamonada; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Metamonada
grep ' root; cellular organisms; Eukaryota; Rhodophyta; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Rhodophyta
grep ' root; cellular organisms; Eukaryota; Haptista; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Haptista
grep ' root; cellular organisms; Eukaryota; Cryptophyceae; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Cryptophyceae
grep ' root; cellular organisms; Eukaryota; Apusozoa; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Apusozoa
grep ' root; cellular organisms; Eukaryota; Eukaryota incertae sedis; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.sedis
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Fungi; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Fungi
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Metazoa; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Metazoa
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Choanoflagellata; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Choanoflagellata
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Filasterea; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Filasterea
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Ichthyosporea; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Ichthyosporea
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Rotosphaerida; ' uniprot_seq_OX_count5000_alllineage.txt | cut -d ':' -f1 > uniprot_seq_OX_count5000_alllineage.txt.Rotosphaerida

grep ' root; Viruses; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Viruses
grep ' root; unclassified entries; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.unclassified
grep ' root; cellular organisms; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.cellular
grep ' root; cellular organisms; Bacteria; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Bacteria
grep ' root; cellular organisms; Archaea; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Archaea
grep ' root; cellular organisms; Eukaryota; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Eukaryota
grep ' root; cellular organisms; Eukaryota; Opisthokonta; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Opisthokonta
grep ' root; cellular organisms; Eukaryota; Viridiplantae; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Viridiplantae
grep ' root; cellular organisms; Eukaryota; Sar; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Sar
grep ' root; cellular organisms; Eukaryota; Discoba; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Discoba
grep ' root; cellular organisms; Eukaryota; Amoebozoa; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Amoebozoa
grep ' root; cellular organisms; Eukaryota; Metamonada; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Metamonada
grep ' root; cellular organisms; Eukaryota; Rhodophyta; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Rhodophyta
grep ' root; cellular organisms; Eukaryota; Haptista; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Haptista
grep ' root; cellular organisms; Eukaryota; Cryptophyceae; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Cryptophyceae
grep ' root; cellular organisms; Eukaryota; Apusozoa; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Apusozoa
grep ' root; cellular organisms; Eukaryota; Eukaryota incertae sedis; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.sedis
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Fungi; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Fungi
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Metazoa; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Metazoa
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Choanoflagellata; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Choanoflagellata
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Filasterea; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Filasterea
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Ichthyosporea; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Ichthyosporea
grep ' root; cellular organisms; Eukaryota; Opisthokonta; Rotosphaerida; ' TGTin_uniprot_alllineage.txt | cut -d ':' -f1 > TGTin_uniprot_alllineage.txt.Rotosphaerida

taxon=(Viruses unclassified cellular Bacteria Archaea Eukaryota Opisthokonta \
Viridiplantae Sar Discoba Amoebozoa Metamonada Rhodophyta Haptista Cryptophyceae Apusozoa sedis \
Fungi Metazoa Choanoflagellata Filasterea Ichthyosporea Rotosphaerida)

for tax in ${taxon[@]} ; do
  awk 'NR==FNR{a[$0];next} $0 in a' uniprot_seq_OX_count5000_alllineage.txt.${tax} TGTin_uniprot_alllineage.txt.${tax} > TGTin_uniprot_OX5000.txt.${tax}
done

for tax in ${taxon[@]} ; do
  wc -l TGTin_uniprot_OX5000.txt.${tax}
done

# taxon id only in OX5000 but not in TGT+ list.
taxon=(Viridiplantae Sar Fungi Metazoa)
for tax in ${taxon[@]} ; do
  grep -Fxv -f TGTin_uniprot_OX5000.txt.${tax} uniprot_seq_OX_count5000_alllineage.txt.${tax} > TGTabs_OX5000_taxid.txt.${tax}
done

# taxon id only in OX5000 and in TGT+ list.
for tax in ${taxon[@]} ; do
  awk 'NR==FNR{a[$0];next} $0 in a' TGTin_uniprot_OX5000.txt.${tax} uniprot_seq_OX_count5000_alllineage.txt.${tax} > TGTpre_OX5000_taxid.txt.${tax}
done

for tax in ${taxon[@]} ; do
  awk 'NR==FNR{a[$0];next} $2 in a' TGTpre_OX5000_taxid.txt.${tax} uniprot_seq_OX_count5000.txt > TGTpre_OX5000_taxid_count.txt.${tax}
done

sort -k1,1nr TGTpre_OX5000_taxid_count.txt.Fungi > TGTpre_OX5000_taxid.txt.Fungi.select
sort -k1,1nr TGTpre_OX5000_taxid_count.txt.Metazoa > TGTpre_OX5000_taxid.txt.Metazoa.select
sort -k1,1nr TGTpre_OX5000_taxid_count.txt.Sar > TGTpre_OX5000_taxid.txt.Sar.select
sort -k1,1nr TGTpre_OX5000_taxid_count.txt.Viridiplantae > TGTpre_OX5000_taxid.txt.Viridiplantae.select

# results:
#1130 TGTpre_OX5000_taxid_count.txt.Metazoa
#208 TGTpre_OX5000_taxid_count.txt.Sar
#303 TGTpre_OX5000_taxid_count.txt.Viridiplantae

# $1 is the number of entries.

#clean the list of taxon ids by tBLASTn, leaving true TGT- organisms by TGT length 300-500 aa.

# II. clustering of protein families.
cd /blue/lagard/yuanyifeng/Qtaxon/sequences

# download sequences for each of taxon ids from Uniprot.
taxons=(Sar Viridiplantae Metazoa Fungi)

for rank in ${taxons[@]} ; do
  sbatch -J wget -e errwget  wget_$rank.sh
done

# example of wget_Fungi.sh

#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=3-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --mem=1gb
###############################################################

for taxid in $(cat TGTpre_OX5000_taxid.txt.Fungi.select); do
  wget -O Fungi/taxid_${taxid}.tsv "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Corganism_id%2Cft_transmem%2Csequence&format=tsv&query=%28taxonomy_id%3A${taxid}%29"
done

for taxid in $(cat TGTabs_OX5000_taxid.txt.Fungi.select); do
  wget -O Fungi/noq_${taxid}.tsv "https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cid%2Cprotein_name%2Cgene_names%2Corganism_name%2Corganism_id%2Cft_transmem%2Csequence&format=tsv&query=%28taxonomy_id%3A${taxid}%29"
done

## END ##

# 2. calculate transmem domain.
for rank in ${taxons[@]} ; do
  # calculate the number and proportion of Transmembrane domains in Q+ taxon ids.
  sbatch -J transmem -e errtransmem transmem_q.sh $rank
  # calculate the number and proportion of Transmembrane domains in Q- taxon ids.
  sbatch -J transmem -e errtransmem transmem_n.sh $rank
done

# transmem_q.sh

#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=3-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
###############################################################

rank=$1

TGTprelist=TGTpre_OX5000_taxid.txt.${rank}.select

for taxid in $(cat ${TGTprelist}) ; do
  awk -F '\t' -v str="TRANSMEM" 'BEGIN{OFS="\t"}FNR==1{print $0,"#ofTRANSMEM4","length"} FNR>1{count=gsub(str, str, $7); print $0"\t"count"\t"length($8)}' ${rank}/taxid_${taxid}.tsv > ${rank}/taxid_${taxid}_transmemcount.tsv.tmp

  awk -F '\t' '{print $7}' ${rank}/taxid_${taxid}_transmemcount.tsv.tmp | while IFS= read -r line; do
    len=0
    if [ -n "$line" ] ; then
      len=$(echo "$line" | sed 's/^TRANSMEM //' | sed 's/; TRANSMEM /\n/g' | sed 's/; .*$//' | sed 's/\.\./#/' | awk -F '#' '{print ($2-$1+1)}' | paste -sd+ | bc)
    fi
    echo "$len" >> ${rank}/taxid_${taxid}_transmemcount.tsv.tmp2
  done

  paste -d '\t' ${rank}/taxid_${taxid}_transmemcount.tsv.tmp ${rank}/taxid_${taxid}_transmemcount.tsv.tmp2 | awk -F '\t' 'BEGIN{OFS="\t"} FNR==1{print "Entry","Entry Name","Protein names","Gene Names","Organism","Organism (ID)","Transmembrane","Sequence","#ofTRANSMEM4","length","length of TRANSMEM","TRANSMEM%"} FNR>1{print $0,($11*100/$10)}' > ${rank}/taxid_${taxid}_transmemcount.tsv

  rm ${rank}/taxid_${taxid}_transmemcount.tsv.tmp ${rank}/taxid_${taxid}_transmemcount.tsv.tmp2
done


## END

# transmem_n.sh

#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=3-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=1
#SBATCH --mem=4gb
###############################################################

rank=$1

for taxid in $(cat TGTabs_OX5000_taxid.txt.${rank}.select) ; do

 awk -F '\t' -v str="TRANSMEM" 'BEGIN{OFS="\t"}FNR==1{print $0,"#ofTRANSMEM4","length"} FNR>1{count=gsub(str, str, $7); print $0"\t"count"\t"length($8)}' ${rank}/noq_${taxid}.tsv > ${rank}/noq_${taxid}_transmemcount.tsv.tmp

 awk -F '\t' '{print $7}' ${rank}/noq_${taxid}_transmemcount.tsv.tmp | while IFS= read -r line; do
   len=0
   if [ -n "$line" ] ; then
     len=$(echo "$line" | sed 's/^TRANSMEM //' | sed 's/; TRANSMEM /\n/g' | sed 's/; .*$//' | sed 's/\.\./#/' | awk -F '#' '{print ($2-$1+1)}' | paste -sd+ | bc)
   fi
   echo "$len" >> ${rank}/noq_${taxid}_transmemcount.tsv.tmp2
 done

 paste -d '\t' ${rank}/noq_${taxid}_transmemcount.tsv.tmp ${rank}/noq_${taxid}_transmemcount.tsv.tmp2 | awk -F '\t' 'BEGIN{OFS="\t"} FNR==1{print "Entry","Entry Name","Protein names","Gene Names","Organism","Organism (ID)","Transmembrane","Sequence","#ofTRANSMEM4","length","length of TRANSMEM","TRANSMEM%"} FNR>1{print $0, ($11*100/$10)}' > ${rank}/noq_${taxid}_transmemcount.tsv

 rm ${rank}/noq_${taxid}_transmemcount.tsv.tmp ${rank}/noq_${taxid}_transmemcount.tsv.tmp2
done

## END

# filter Q+ organisms by >= 4 Transmembrane domains covering >= 25% protein length.
taxons=(Sar Viridiplantae Metazoa Fungi)

for rank in ${taxons[@]} ; do
  for taxid in $(cat TGTpre_OX5000_taxid.txt.${rank}.select) ; do
    awk -F '\t' '$9>=4&&$12>=25 {print $9}' ${rank}/taxid_${taxid}_transmemcount.tsv | wc -l
  done
  for taxid in $(cat TGTabs_OX5000_taxid.txt.${rank}.select) ; do
    awk -F '\t' '$9>=4&&$12>=25 {print $9}' ${rank}/noq_${taxid}_transmemcount.tsv | wc -l
  done
done

# set filter of number of sequences depending on the taxon rank.
# copy and paste taxon ids to file : TGTpre_OX5000_taxid.txt.${rank}.select and TGTabs_OX5000_taxid.txt.${rank}.select.

# extract sequences as fasta format.
for rank in ${taxons[@]} ; do
  for taxid in $(cat TGTpre_OX5000_taxid.txt.${rank}.select) ; do
    awk -F '\t' -v t=$taxid '$9>=4&&$12>=25 {print ">q_"t"_"$1"\n"$8}' ${rank}/taxid_${taxid}_transmemcount.tsv > ${rank}/q_${taxid}_transmem4p25.fasta
  done
  for taxid in $(cat TGTabs_OX5000_taxid.txt.${rank}.select) ; do
    awk -F '\t' -v t=$taxid '$9>=4&&$12>=25 {print ">n_"t"_"$1"\n"$8}' ${rank}/noq_${taxid}_transmemcount.tsv > ${rank}/n_${taxid}_transmem4p25.fasta
  done
done

# compile sequences.
for rank in ${taxons[@]} ; do
  > ${rank}/${rank}_merge_transmem4p25.fasta
  for taxid in $(cat TGTpre_OX5000_taxid.txt.${rank}.select); do
    cat ${rank}/q_${taxid}_transmem4p25.fasta >> ${rank}/${rank}_merge_transmem4p25.fasta
  done
  for taxid in $(cat TGTabs_OX5000_taxid.txt.${rank}.select); do
    cat ${rank}/n_${taxid}_transmem4p25.fasta >> ${rank}/${rank}_merge_transmem4p25.fasta
  done
done

# clustering using mmseqs2 with identity and coverage from 0.05 to 0.5 to determine the cutoff of identity and coverage.

for rank in ${taxons[@]} ; do
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.05
  sleep 0.5
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.1
  sleep 0.5
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.2
  sleep 0.5
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.3
  sleep 0.5
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.4
  sleep 0.5
  sbatch -J ${rank} -e errmmseq${rank} mmseqs_transmem4p25.sh ${rank} 0.5
  sleep 0.5
done

# mmseqs_transmem4p25.sh

#!/bin/bash
#SBATCH --mail-type=END,FAIL          # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=yuanyifeng@ufl.edu     # Where to send mail
#SBATCH --nodes=1
#SBATCH --ntasks=1                    # Run on a single CPU
#SBATCH --time=3-00:00:00               # Time limit hrs:min:sec
#SBATCH --cpus-per-task=4
#SBATCH --mem=10gb
###############################################################

module load mmseqs2/14

rank=$1
job=transmem4p25
ct=$2

fas=${rank}/${rank}_merge_${job}.fasta

mmseqs easy-cluster ${fas} mmseq_${rank}_${job}_${ct} tmp_${rank}_${job}_${ct} \
        --min-seq-id ${ct} -c ${ct} --cov-mode 0 \
         --threads 4 \
         #--cluster-reassign 1

#### END

# III.
# count the occurence of Q+ and Q- taxon ids in each cluster.
for rank in ${taxons[@]} ; do
  for ct in 0.05 0.1 0.2 0.3 0.4 0.5; do
    sbatch -J ${rank} -e errmmseq${rank} --wrap "python3 score_taxid.py mmseq_${rank}_transmem4p25_${ct}_cluster.tsv ${rank}_transmem4p25_${ct}_rep.txt"
    sleep 0.3
  done
done

# calculate scores in Excel.
