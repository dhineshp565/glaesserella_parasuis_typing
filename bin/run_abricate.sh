#!/bin/bash
# $1 = SampleName
# $2 = genome assembly consensus
# $3 = directory containing databases
# Usage: ./run_abricate.sh <SampleName> <ConsensusFASTA> <DBDir>

set -euo pipefail

# Define default values
DefaultLine="$1\t$1_contig_1\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone\tNone"
HeaderOnly=1  # Expected line count if only the header is present

abricate --datadir $3 --db Gparasuis_serodb_Howell -minid 80  -mincov 80 --quiet $2 1> $1_serotype.csv

if [ "$(wc -l < "$1_serotype.csv")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$1_serotype.csv"
fi
abricate --datadir $3 --db Gparasuis_serodb_Jia -minid 80  -mincov 80 --quiet $2 1> $1_sero_Jia.csv
if [ "$(wc -l < "$1_sero_Jia.csv")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$1_sero_Jia.csv"
fi

abricate -datadir $3 --db Gparasuis_vfdb $2 1> $1_vf.csv

if [ "$(wc -l < "$1_vf.csv")" -eq $HeaderOnly ]; then
    	echo -e "$DefaultLine" >> "$1_vf.csv"
fi

abricate --db card $2 1> $1_AMR.csv

if [ "$(wc -l < "$1_AMR.csv")" -eq $HeaderOnly ]; then
    echo -e "$DefaultLine" >> "$1_AMR.csv"
fi
