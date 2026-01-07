awk '
{
    pos = index($0, "in")
    if(pos != 0) {
       print substr($0 , pos+1)
    } else{
       print $0
    }
}
' A549_RNASEQ_hg38_nameid_modify.txt >A549_RNASEQ_hg38_nameid_modify_1.txt
