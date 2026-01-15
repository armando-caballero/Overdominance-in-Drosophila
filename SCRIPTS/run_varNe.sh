# Path to a new folder
PATHR="/Users/sramos/Desktop/Armando-Overdominance2025/run_script"

mkdir ${PATHR}/Variability

# ################################
# These 4 files must be in the folder
VCF_FILE="Pb-B-ALL_merged_snps_excised_autosomic.nr_annovar.dm6nr_multianno.vcf"
REF_FILE="dmel-all-chromosome-r6.18.fasta"
FAI_FILE="dmel-all-chromosome-r6.18.fasta.fai"
GFF_FILE="dmel-all-r6.18.gtf"
# Include "calculate_sfs_mstatspop.R" in the folder ${PATHR}/Variability:
# ################################

OUT_FILE="Drosophila_Arm"

# Download software
git clone https://github.com/sramosonsins/gVCF2tFasta.git
cd ${PATHR}/gVCF2tFasta/
sh ./build.sh
cd ..

git clone https://github.com/CRAGENOMICA/fastaconvtr.git
cd ${PATHR}/fastaconvtr/
sh ./build.sh
cd ..

git clone https://github.com/CRAGENOMICA/mstatspop.git
cd ${PATHR}/mstatspop
sh ./build.sh
cd ..

#copy executables into folder /Variability
cp ${PATHR}/gVCF2tFasta/build/gVCF2tFasta ${PATHR}/Variability/
cp ${PATHR}/fastaconvtr/build/fastaconvtr ${PATHR}/Variability/
cp ${PATHR}/mstatspop/build/mstatspop ${PATHR}/Variability/
cp ${PATHR}/mstatspop/bin/collect_data_columns.pl ${PATHR}/Variability/

# ################################
# Estimate Variability
# ################################

cd ${PATHR}/Variability/

# convert to TFA
./gVCF2tFasta -v ${PATHR}/${VCF_FILE} -r ${PATHR}/${REF_FILE} -c 2 -o ${PATHR}/${OUT_FILE} -n ${PATHR}/${REF_FILE}.fai -i 1

# Obtain Silent/Syn/NSyn positions
./fastaconvtr -i ${PATHR}/${OUT_FILE}.tfa.gz -F t -f t -o ${PATHR}/${OUT_FILE}B.tfa.gz.tfa.gz -n ${PATHR}/${REF_FILE}.fai -t ${PATHR}/${OUT_FILE}.tfa.gz.Silent.WEIGHT.txt.gz -g ${PATHR}/dmel-all-r6.18.gtf silent        Nuclear_Universal -T 0
./fastaconvtr -i ${PATHR}/${OUT_FILE}.tfa.gz -F t -f t -o ${PATHR}/${OUT_FILE}B.tfa.gz.tfa.gz -n ${PATHR}/${REF_FILE}.fai -t ${PATHR}/${OUT_FILE}.tfa.gz.Syn.WEIGHT.txt.gz    -g ${PATHR}/dmel-all-r6.18.gtf synonymous    Nuclear_Universal -T 0
./fastaconvtr -i ${PATHR}/${OUT_FILE}.tfa.gz -F t -f t -o ${PATHR}/${OUT_FILE}B.tfa.gz.tfa.gz -n ${PATHR}/${REF_FILE}.fai -t ${PATHR}/${OUT_FILE}.tfa.gz.Nsyn.WEIGHT.txt.gz   -g ${PATHR}/dmel-all-r6.18.gtf nonsynonymous Nuclear_Universal -T 0

# obtain fSFS and Variability
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 2 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_TOTAL_Variability_o2.txt  -N 1 102 -u 0 -G 0 -w 100000000 
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 2 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_Silent_Variability_o2.txt -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.Silent.WEIGHT.txt.gz
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 2 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_Syn_Variability_o2.txt    -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.Syn.WEIGHT.txt.gz
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 2 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_NSyn_Variability_o2.txt   -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.NSyn.WEIGHT.txt.gz

./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 0 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_TOTAL_Variability_o0.txt  -N 1 102 -u 0 -G 0 -w 100000000 
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 0 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_Silent_Variability_o0.txt -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.Silent.WEIGHT.txt.gz
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 0 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_Syn_Variability_o0.txt    -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.Syn.WEIGHT.txt.gz
./mstatspop -f tfa -i ${PATHR}/${OUT_FILE}.tfa.gz -o 0 -n ${PATHR}/${REF_FILE}.fai -T ${PATHR}/${OUT_FILE}.tfa_NSyn_Variability_o0.txt   -N 1 102 -u 0 -G 0 -w 100000000 -E ${PATHR}/${OUT_FILE}.tfa.gz.NSyn.WEIGHT.txt.gz

# create a file called "cols.txt" containing the names of statistics to keep (here all)

# sort into a tabulated table
perl collect_data_columns.pl -in ${OUT_FILE}.tfa_NSyn_Variability_o2.txt   -fc cols.txt > ${OUT_FILE}.tfa_NSyn_Variability_o2.COLS.txt
perl collect_data_columns.pl -in ${OUT_FILE}.tfa_Syn_Variability_o2.txt    -fc cols.txt > ${OUT_FILE}.tfa_Syn_Variability_o2.COLS.txt
perl collect_data_columns.pl -in ${OUT_FILE}.tfa_Silent_Variability_o2.txt -fc cols.txt > ${OUT_FILE}.tfa_Silent_Variability_o2.COLS.txt
perl collect_data_columns.pl -in ${OUT_FILE}.tfa_TOTAL_Variability_o2.txt  -fc cols.txt > ${OUT_FILE}.tfa_TOTAL_Variability_o2.COLS.txt

# estimate global variability and plot fSFS
R --vanilla < ./calculate_sfs_mstatspop.R --args ${OUT_FILE}
cd ..

# ################################
# Estimate Ne using fSFS
# ################################

#first create Drosophila_blueprint file. 
#Include fSFS and the rest of parameters
#files:
# ./Drosophila_Silent.blueprint for fSFS Silent
# ./Drosophila_Syn.blueprint for fSFS Synonymous

#download software
git clone https://github.com/xiaoming-liu/stairway-plot-v2.git
cd ${PATHR}/stairway-plot-v2/stairway_plot_v2.1.3
#include Blueprint files into this folder

java -cp stairway_plot_es Stairbuilder ./Drosophila_Silent.blueprint
sh Drosophila_Silent.blueprint.sh

java -cp stairway_plot_es Stairbuilder ./Drosophila_Syn.blueprint
sh Drosophila_Syn.blueprint.sh

R --vanilla < ./plot_stairwayplot2_mstatspop.R
cd ../..