for bam in *_mapped.bam; do
    mapped=$(samtools view -c -F 4 "$bam")
    total=$(samtools view -c "$bam")
    unmapped=$((total - mapped))
    echo -e "$bam\t$total\t$mapped\t$unmapped"
done
