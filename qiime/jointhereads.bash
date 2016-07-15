#in the microbiomeQIIME all the sequence start with Argonne 
#autism_mic_backup
#/bin/bash: i: command not found
#indir will be /scratch/users/mmdavid/autism_mic_backup
#outdir will be /scratch/users/mmdavid/microbiome_Qiime
#OR /scratch/users/ctataru5/microbiome/qiime
while getopts "i:o:"; do
	case "$opt" in
		i) 
			indir=$OPTARG
			;;
		o) 
			outdir=$OPTARG
			;;
	esac
done
cd $indir
for folder in A*; do
	if  [! -f $outdir/${folder%.}.paired];
    	then 
	join_paired_ends.py -f $indir/${folder}/*R1* -r $indir/${folder}/*R2* -o $indir/${folder%.}.paired/ -b ${folder}/*I1*
	fi
done

echo "Reads join finished"
now=$(date +"%m_%d_%Y")
#this script assume that you are running ONE batch of new sequences at a time 
cat $indir/*.paired/*join.fastq > $outdir/cat.fastq.paired.$now  

