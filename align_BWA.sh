#!/bin/bash

#SBATCH -n 2                        # number of cores                                                         
#SBATCH -t 0-12:00                  # wall time (D-HH:MM)                                                     
##SBATCH -A ssnyde11                # Account hours will be pulled from (commented out with double # in front)                                                                                                            
#SBATCH -o slurm.%j.out             # STDOUT (%j = JobId)                                                     
#SBATCH -e slurm.%j.err             # STDERR (%j = JobId)                                                     
#SBATCH --mail-type=ALL             # Send a notification

module load star/2.5.2b


baseDir=/home/ssnyde11/scratch/Daphnia_RNAseq_090719
fileDir=/home/ssnyde11/scratch/Daphnia_RNAseq_090719/tagdusted
BWAdir=BWAdir
####### Before running the script, please enter path to desired output directory, below ####
WD=/home/ssnyde11/scratch/
BWA=/packages/spack/linux-centos7-x86_64/gcc-4.8.5/bwa-0.7.17-fusrahsomx5s4f64jk455gpi7i6x\
gaza/bin/bwa
SAMTOOLS=/packages/7x/samtools/1.9/bin/samtools
fqDir=fastqs
GENOME_DIR=/home/ssnyde11/scratch/rnaseq_genomes
GENOME_FILE=PA42.4.1.fasta
THREADS=16


echo "Starting job"

echo "Making symbolic links to fastq files."

cd $WD

if [ ! -d "$fqDir" ] ; then
    mkdir $fqDir
  fi

cd $fqDir

#update symlinks here
ln -s ${fileDir}/GSF2349-LPB-2014-32-5_S5_trno_tagdusted.fq LPB-2014-32-5.fq
ln -s ${fileDir}/GSF2349-LPB-2014-32-2-new_S6_trno_tagdusted.fq LPB-2014-32-2.fq
ln -s ${fileDir}/GSF2349-LPA-2014-32-5-new_S4_trno_tagdusted.fq LPA-2014-32-5.fq
ln -s ${fileDir}/GSF2349-LPA-2014-16-REP2_S2_trno_tagdusted.fq LPA-2014-16-REP2.fq
ln -s ${fileDir}/GSF2349-LPA-2014-16-REP1_S1_trno_tagdusted.fq LPA-2014-16-REP1.fq
ln -s ${fileDir}/GSF2349-LPA-2014-16-3_S3_trno_tagdusted.fq LPA-2014-16-3.fq
ln -s ${fileDir}/GSF2349-POV-2014-12-REP4_S12_trno_tagdusted.fq POV-2014-12-REP4.fq
ln -s ${fileDir}/GSF2349-POV-2014-12-REP3_S11_trno_tagdusted.fq POV-2014-12-REP3.fq

ln -s ${fileDir}/GSF2349-POV-2014-12-REP2_S10_trno_tagdusted.fq POV-2014-12-REP2.fq
ln -s ${fileDir}/GSF2349-POV-2014-12-REP1_S9_trno_tagdusted.fq POV-2014-12-REP1.fq
ln -s ${fileDir}/GSF2349-NFL3-REP4_S20_trno_tagdusted.fq NFL3-REP4.fq
ln -s ${fileDir}/GSF2349-NFL3-REP3_S19_trno_tagdusted.fq NFL3-REP3.fq
ln -s ${fileDir}/GSF2349-NFL3-REP2_S18_trno_tagdusted.fq NFL3-REP2.fq
ln -s ${fileDir}/GSF2349-NFL3-REP1_S17_trno_tagdusted.fq NFL3-REP1.fq
ln -s ${fileDir}/GSF2349-LPB-2014-32-4-old_S8_trno_tagdusted.fq LPB-2014-32-4-old.fq
ln -s ${fileDir}/GSF2349-LPB-2014-32-3_S7_trno_tagdusted.fq LPB-2014-32-3.fq
ln -s ${fileDir}/GSF2349-KAP-2013-114-REP1_S13_trno_tagdusted.fq KAP-2013-114-REP1.fq
ln -s ${fileDir}/GSF2349-KAP-2013-114-4_S16_trno_tagdusted.fq KAP-2013-114-4.fq
ln -s ${fileDir}/GSF2349-KAP-2013-114-3_S15_trno_tagdusted.fq KAP-2013-114-3.fq

cd $WD

#source 0README

cd $GENOME_DIR

echo "Indexing the Daphnia genome using bwa ..."
echo "bwa index ${GENOME_DIR}/${GENOME_FILE}"
#${BWA}    index ${GENOME_DIR}/${GENOME_FILE}

echo "Performing alignment on the STRIPE-seq files"

cd ${baseDir}

if [ ! -d "$BWAdir" ] ; then
    mkdir $BWAdir
  fi

cd $BWAdir

echo "Starting alignments ..."
for fq in ${WD}/${fqDir}/*.fq;
do

        echo "bwa aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${fq} .fq).sai ${fq}"
  	${BWA} aln -t ${THREADS} -n 3 ${GENOME_DIR}/${GENOME_FILE} -f $(basename ${fq} .fq).sai ${fq} 

	echo "$fq"

echo "${BWA} samse ${GENOME_DIR}/${GENOME_FILE} $(basename $fq .fq).sai $fq | ${SAMTOOLS} view -uS - | ${SAMTOOLS} sort -O BAM - > $(basename $fq .fq)_sorted.bam"
${BWA} samse ${GENOME_DIR}/${GENOME_FILE} $(basename $fq .fq).sai $fq | ${SAMTOOLS} view -uS - | ${SAMTOOLS} sort -O BAM - > $(basename $fq .fq)_sorted.bam

echo "samtools index -b $(basename $fq .fq)_sorted.bam "
${SAMTOOLS} index -b $(basename $fq .fq)_sorted.bam

#FILTERED_BAM=$(basename $fq .fq)_filtered.bam
#.. post-alignment filtering for proper alignments and MAPQ >= 10:
#
echo "${SAMTOOLS} view -f 2 -q 10 -u ${SORTED_BAM} | ${SAMTOOLS} sort -O BAM -@ 10 - > $(basename $fq .fq)_filtered.bam"
${SAMTOOLS} view -h -F 4 -q 10 -u $(basename $fq .fq)_sorted.bam | ${SAMTOOLS} sort -O BAM -@ 8 - > $(basename $fq .fq)_filtered.bam

echo "samtools index -b $(basename $fq .fq)_filtered.bam"
${SAMTOOLS} index -b $(basename $fq .fq)_filtered.bam

done

  echo ""
  echo " Done with step 3 (read mapping)."
  echo ""
  echo "================================================================================"


exit


