#! /bin/bash

#SBATCH --job-name=07
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00



for i in *_rmRep_sorted.bam;do
	n=${i/_rmRep_sorted.bam/}
        gatk --java-options "-XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=${n}_rmRep_sorted.bam --OUTPUT=${n}_rmRep_gatksorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true
done
