## An example that only works for the luria server.
# module load blast/2.6.0 
# export PERL5LIB=/home/yfyuan/perl5/lib/perl5 
# module load perl/5.24.1 

## build blastDB
# path_HGM_genomes=~/data/genomes/genome13663fasta
#ls ${path_HGM_genomes}/*.fasta | xargs -i makeblastdb -in {} -dbtype nucl -parse_seqids -out {}_blastdb 
#>&2 echo "db done" 

 

## tblastn
qfaa=$1
mkdir outfiles || true
for db in $(ls ${path_HGM_genomes}/*fasta) ; do 
  bname=$(basename "$db" | sed 's/\..*//' ) 
  out=outfiles/${bname}.out 
  tblastn -query "$qfaa" -db "$db"_blastdb -out "$out" \
         -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" \
         -num_threads 10
done
