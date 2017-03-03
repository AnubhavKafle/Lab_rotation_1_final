mkdir scripts
curl -o scripts/firehose_get_latest.zip http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip
unzip -d scripts scripts/firehose_get_latest.zip
./scripts/firehose_get -b -only Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data data latest
mkdir rna_seq
for file in stddata__*/*/*/*RSEM*.Level_3*.tar.gz; do tar -zxf $file -C rna_seq; done
ls -d rna_seq/gdac* | perl -ne 'chomp; ($t)=m/gdac.broadinstitute.org_(\w+)/; print "mv $_ rna_seq/$t\n"' | bash
rm -rf {rna_seq,rppa_expr}/{COAD,READ}

