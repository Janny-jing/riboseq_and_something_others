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
