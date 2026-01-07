#! /bin/bash

#SBATCH --job-name=05
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --time=2-00:00:00



for i in *.bam;do
	n=${i/.bam/}
        gatk --java-options "-XX:+UseParallelOldGC -XX:ParallelGCThreads=4 -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" SortSam --INPUT=${n}.bam --OUTPUT=${n}_sorted.bam --VALIDATION_STRINGENCY=SILENT -SO=coordinate --CREATE_INDEX=true
done
