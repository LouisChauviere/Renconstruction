#!/bin/bash

##########################################
#     Louis Chauviere 07/03/2018         #
#       Download SRA data                #
# transform sra files in 2 fastq files   #
#    quality control using fastq         #
#     read mapping using HISAT2          #
#          sort mapped reads             #       
#           index bam files              #
#  gene reconstruction using stringtie   #
#  gene reconstruction using scallop     #
##########################################

#Reading necessary arguments from config file
if [ $# -lt 1 ]; then
  echo "Please give configuration file as a command line argument"
  echo "Usage: get_sra.sh <configuration file>"
  exit 3
else
  source $1
fi



#search samples ID : return a list
sra=($(awk '{ if (NR!=1) { print $10 } }' $metadata))

echo ${sra[@]}
echo "SRA download begining"

if [ ! -d qc_results ]; then #if directory doesn' exist
	mkdir qc_results
fi


echo $metadata
listMerge="${experiment}.bam"

#SRA download + fastq creation + mapping + reconstruction + comparison stringtie/scallop
for sample in "${sra[@]}" 
do
	#if [ ${sample} == "plouf" ]; then
	logfile=${sample}.log		
	outcompare=gffrefseq${sample}_scallop
	beginSample=($(echo $sample | awk -F '' '{print $1$2$3$4$5$6}')) 
	#if we need to download sra files
	if [ -f "${sample}.sra" ]
		then echo "${sample}.sra exists";
	else
		wget ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/ERR/${beginSample}/${sample}/${sample}.sra >> $logfile;
	fi
	echo "###############################"
	echo "#                             #"
	echo "# sra to fastq -> fastq-dump  #" 
	echo "# For paired-end experiments  #"
	echo "#                             #"
	echo "###############################"
	fastq-dump -I --split-files ${sample}.sra >> $logfile
	echo "SRA to fastq done"
	rm ${sample}.sra
	echo "fastqc for the two fastq files"
	fastqc -o qc_results ${sample}_1.fastq >> $logfile
	fastqc -o qc_results ${sample}_2.fastq >> $logfile
	echo "fastqc done"
	##echo "###########################"
	##echo "copy to cluster"
	##sshpass -p $password scp ${sample}_1.fastq moideen@goliath.sdv.univ-paris-diderot.fr:/scratch/user/moideen/k562/fastq/
	##sshpass -p $password scp ${sample}_2.fastq moideen@goliath.sdv.univ-paris-diderot.fr:/scratch/user/moideen/k562/fastq/
	##echo "copy done"
	##rm ${sample}_1.fastq
	##rm ${sample}_2.fastq
	echo "##########################"
	echo "#                        #"
	echo "# start HISAT2 alignment #"
	echo "#                        #"
	echo "##########################"
	hisat2 -p 8 --dta --novel-splicesite-outfile ${sample}_novel_splice_sites_HISAT2.txt -x ${hisat_index} -1 ${sample}_1.fastq -2 ${sample}_2.fastq -S ${sample}_hg19.sam --summary-file /summary_HISAT2/${sample}.summary.txt >> $logfile
	echo "########################################"
	echo "#                                      #"
	echo "# start samtools sort & samtools index #"
	echo "#                                      #"
	echo "########################################"
	samtools sort -@ 8 -o ${sample}_hg19.bam ${sample}_hg19.sam >> $logfile
	samtools index ${sample}_hg19.bam >> $logfile
	#rm ${sample}_hg19.sam
	echo "HISAT + samtools DONE"
	echo "#####################################"
	echo "#                                   #"
	echo "#  BEGIN stringtie reconstruction   #"
	echo "#                                   #"
	echo "#####################################"
	stringtie -p 8 -o stringtie_${sample}_hg19.gtf -l ${sample} ${sample}_hg19.bam >> $logfile
	echo "stringtie reconstruction DONE"
	
	echo "########################"
	echo "#                      #"
	echo "#    BEGIN scallop     #"
	echo "#                      #"
	echo "########################"
	~/ngs_bin/scallop-0.10.2_linux_x86_64/scallop -i ${sample}_hg19.bam -o scallop_${sample}_hg19.gtf --min_single_exon_coverage $min_single_exon_coverage >> $logfile
	gffcompare -o $outcompare -r $refgtf scallop_${sample}_hg19.gtf
	echo "scallop DONE"
	echo "all DONE"
	#fi
	listMerge="${listMerge} ${sample}_hg19.bam"
	#echo "?????????????????????,,;:zeofpzkeàoeafkezfk"	
	echo $sample
	echo $listMerge
done

echo "##############################"
echo "#                            #"
echo "# merging and reconstruction #"
echo "#                            #"
echo "##############################"

echo "merge all bam files"
samtools merge -f ${listMerge}
echo "sort merged bam file"
samtools sort -@ 8 -o ${experiment}.sorted.bam ${experiment}.bam
echo "index sorted bam file"
samtools index ${experiment}.sorted.bam
#rm ${experiment}.bam
echo "Merging DONE"

echo "reconstruction of merged files with stringtie and scallop"
stringtie -p 8 -o stringtie_${experiment}_merged_hg19.gtf ${experiment}.sorted.bam >> $logfile
~/ngs_bin/scallop-0.10.2_linux_x86_64/scallop -i ${experiment}.sorted.bam -o scallop_${experiment}_merged_hg19.gtf --min_single_exon_coverage $min_single_exon_coverage >> $logfile
gffcompare -o gffrefseq_${experiment}_scallop -r $refgtf scallop_${experiment}_merged_hg19.gtf

echo "all DONE"
