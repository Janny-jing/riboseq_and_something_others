#去重
#! /bin/sh

awk '

BEGIN { FS="\n"; OFS="\n" }

{

    lines[NR] = $0

    if (NR % 22 == 0) {

        current_line = lines[NR]

        if (!seen[current_line]) {

            for (i = NR - 21; i <= NR; i++) {

                print lines[i]

            }

            seen[current_line] = 1

        }

    }

}

' GRCh38_cds_name_modify.txt > GRCh38_cds_name_modify_modified.txt


##取值
grep -B21 --binary-files=text -wFf gene_name.txt GRCh38_cds_name_modify_modified.txt >A549_RNASEQ_hg38_nameid_modify.txt

awk '!/--/'  A549_RNASEQ_hg38_nameid_modify_1.txt >A549_RNASEQ_hg38_nameid_modify_2.txt
awk 'NF'  A549_RNASEQ_hg38_nameid_modify_pic_calculation.txt > temp && mv temp A549_RNASEQ_hg38_nameid_modify_pic_calculation.txt
awk '!(NR % 17 == 0) {print $0}' A549_RNASEQ_hg38_nameid_modify_ctl.txt > output_file





grep -B21 -iwFf all_de.txt ../mm39_cds_modify_2.txt >1.txt
awk '!/--/' 1.txt >2.txt && awk 'NF'  2.txt >3.txt && awk '!(NR % 17 == 0) {print $0}' 3.txt >4.txt
cut -f1,2 4.txt >5.txt && cut -f3,4 4.txt >> 5.txt && cut -f5,6 4.txt >> 5.txt && cut -f7,8 4.txt >> 5.txt
mkdir all01 && cp 6.txt all01
python ../../06.corr/genecodons_trnas_anticodons_mouse.py all01 all02
python ../../06.corr/03.merge_mrna_trna_anticodons.py all02 ../trna/BMDM_LPS.txt alllps
python ../../06.corr/03.merge_mrna_trna_anticodons.py all02 ../trna/BMDM_LPS_ctl.txt allctl

python ../../06.corr/eachgene_expression.py 2.txt eachgene01
cd eachgene01 &&  mkdir ../eachgene02 && mv *_percentages.txt ../eachgene02 && cd ..
python ../../06.corr/genecodons_trnas_anticodons_mouse.py eachgene02 eachgene03
python ../../06.corr/03.merge_mrna_trna_anticodons.py eachgene03 ../trna/BMDM_LPS.txt  eachgene04lps
python ../../06.corr/03.merge_mrna_trna_anticodons.py eachgene03 ../trna/BMDM_LPS_ctl.txt eachgene04ctl

