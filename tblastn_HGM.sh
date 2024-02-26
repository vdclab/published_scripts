## An example that only works for the luria server.
# module load blast/2.6.0 
# export PERL5LIB=/home/yfyuan/perl5/lib/perl5 
# module load perl/5.24.1 

## 1. build blastDB.
path_HGM_genomes=~/data/genomes/genome13663fasta
ls ${path_HGM_genomes}/*.fasta | xargs -i makeblastdb -in {} -dbtype nucl -parse_seqids -out {}_blastdb 
>&2 echo "db done" 

 
## 2. tblastn.
qfaa=$1
mkdir outfiles || true
for db in $(ls ${path_HGM_genomes}/*fasta) ; do 
  bname=$(basename "$db" | sed 's/\..*//' ) 
  out=outfiles/${bname}.out 
  tblastn -query "$qfaa" -db "$db"_blastdb -out "$out" \
         -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids salltitles" \
         -num_threads 10
done

## 3. re-format and summarize output.
fsum=summary_tblastn.txt

echo -e "genomeID\c" > $fsum
for Qgene in $(cat Q_genes_list.txt) ; do
  echo -e "\t$Qgene\c" >> $fsum
done

echo "" >> $fsum

for outf in $(ls outfiles/*out); do
  gname=$(basename "$outf" | sed 's/\..*//' )
  echo -e "$gname\c" >> $fsum
  for Qgene in $(cat Q_genes_list.txt) ; do
    n=$(grep -v '^#'  $outf  | awk -F '\t' -v g=$Qgene '$1==g&&$3>=20&&$11<=0.0000000001{print $1}' | wc -l)
    echo -e "\t$n\c" >> $fsum
  done
  echo "" >> $fsum
done
