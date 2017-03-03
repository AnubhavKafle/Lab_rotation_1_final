#To download raw counts from Firehose data portal using their python script

mkdir scripts
curl -o scripts/firehose_get_latest.zip http://gdac.broadinstitute.org/runs/code/firehose_get_latest.zip
unzip -d scripts scripts/firehose_get_latest.zip
./scripts/firehose_get -b -only mRNAseq_Preprocess.Level_3  data latest


#To extract the raw_count files using shell script from the tar

for i in $(ls); do  tar -xzf $i/20160128/*.tar.gz --wildcards --no-anchored '*raw_counts*'; done

#To move files in the respective name folder
ls -d ./gdac* | perl -ne 'chomp; ($t)=m/gdac.broadinstitute.org_(\w+)/; print "mv $_ ./$t\n"' | bash

#To move all raw_counts.txt in one folder 
for i in $(find -iname *uncv2*.txt) ; do mv $i ./Raw_counts/; done

#Remove Header ‘HYBRIDIZATION R’ and replace by ‘Gene_symbol’ so that there is no interference of the space in R

for i in $(ls *.txt) ; do sed -i 's/HYBRIDIZATION R/Gene_symbol/g' $i; done


