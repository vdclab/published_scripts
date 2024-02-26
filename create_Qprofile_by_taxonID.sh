## for each Q protein, download InterPro entries.
IPRID=IPR004803
protein=TGT
python3 download_${IPRID}.py > ${protein}_entry.txt

# re-format output, keep taxonID only.
awk -F '\t' '{print $1}' ${protein}_entry.txt > ${protein}_uniprotid.txt

# merge by ids for each uniport id (column 4).
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' FolE1_entry.txt > FolE1_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' FolE2_entry.txt > FolE2_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueD_entry.txt > QueD_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueE_entry.txt > QueE_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueC_entry.txt > QueC_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' TGT_entry.txt > TGT_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueA_entry.txt > QueA_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueG_entry.txt > QueG_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueH_entry.txt > QueH_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueK_entry.txt > QueK_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' DUF2419_entry.txt > Qng1_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' YhhQ_entry.txt > YhhQ_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QueT_entry.txt > QueT_taxon.txt
awk -F '\t' 'BEGIN {OFS="\t";print "taxon","family"} {print $4,$1}' QrtT_entry.txt > QrtT_taxon.txt

# 1 remove duplicated taxon.
Qgenes=(FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH Qng1 QueK QrtT YhhQ QueT) 

for qgene in ${Qgenes[*]}; do
  awk '!a[$1]++' ${qgene}_taxon.txt > ${qgene}_taxon_uniq.txt
done

# 2 merge using join.
# join FolE1 and FolE2
out=FolE1_taxon_uniq.txt
outnew=FolE12_join.txt
echo 'taxon FolE1 FolE2'
true > ${outnew}
join -t $'\t' -o 0,1.2,2.2 -a1 -a2 -e na <(sort $out) <(sort FolE2_taxon_uniq.txt) >> ${outnew} 

# append_QueD.
appd=QueD_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueD_join.txt
echo 'taxon FolE1 FolE2 QueD'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,2.2 -a1 -a2 -e na $out <(sort ${appd}) >> ${outnew}

# append_QueE.
appd=QueE_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDE_join.txt
echo 'taxon FolE1 FolE2 QueD QueE'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueC.
appd=QueC_taxon_uniq.txt 
out=${outnew}
outnew=FolE12_QueDEC_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueF.
appd=QueF_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew} 

# append_TGT.
appd=TGT_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueA.
appd=QueA_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueA_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueG.
appd=QueG_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAG_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueH.
appd=QueH_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueK.
appd=QueK_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_QueK_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH QueK'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_Qng1.
appd=Qng1_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_QueKQng1_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH QueK Qng1'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QrtT.
appd=QrtT_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_QueKQng1_QrtT_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH QueK Qng1 QrtT'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_QueT.
appd=QueT_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_QueKQng1_QrtTQueT_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH QueK Qng1 QrtT QueT'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}

# append_YhhQ.
appd=YhhQ_taxon_uniq.txt
out=${outnew}
outnew=FolE12_QueDEC_QueF_TGT_QueAGH_QueKQng1_QrtTQueTYhhQ_join.txt
echo 'taxon FolE1 FolE2 QueD QueE QueC QueF TGT QueA QueG QueH QueK Qng1 QrtT QueT YhhQ'
true > ${outnew}
join -t $'\t' -o 0,1.2,1.3,1.4,1.5,1.6,1.7,1.8,1.9,1.10,1.11,1.12,1.13,1.14,1.15,2.2 -a1 -a2 -e na <(cat $out) <(sort ${appd}) >> ${outnew}
