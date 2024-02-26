## An example that only works for the hipergator server.
# module load blast/2.6.0

## 1. build blastDB.
db=/blue/lagard/yuanyifeng/HOM/allfaa/ALL_genomes.faa
ls ${db} | xargs -i makeblastdb -in {} -dbtype prot -parse_seqids -out {}_blastdb
>&2 echo "db done" 

 
## 2. blastp.
qfaa=$1

out=Qgene_in_hom.out
blastp -query "$qfaa" -db "$db"_blastdb -out "$out" \
         -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" \
         -num_threads 10
