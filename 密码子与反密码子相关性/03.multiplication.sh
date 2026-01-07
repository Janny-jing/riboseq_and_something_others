awk '
NR==FNR {
    fileA_values[$1] = $2
    next
}
{
    if ($1 in fileA_values) {
        for (i = NR - 21; i <= NR; i++) {
            if (i in lines) {
                print lines[i]
                    for (j = 2; j <= 8; j += 2) {
                        printf "%.2f\t", fileA_values[$1] * $(j)
                    }
            }
        }
        print $1 "\t" fileA_values[$1]
        lines[NR] = $0
    } else {
        print $0
    }
}
' A549_RNASEQ_readcount_ctl.txt A549_RNASEQ_hg38_nameid_modify_1.txt >A549_RNASEQ_hg38_nameid_modify_ctl.txt

