#!/bin/bash
#$ -cwd

rm script_SLIM4_3L_total_bottle.sh.*

if [ $# -ne 2 ]
then
	echo "Usage: $0 <slim3INPUTfile> <REPS>"
	exit 1
fi

#Set arguments
INPUT=$1
REPS=$2

#Working directory
WDIR=$PWD

#Output directory
if [ -d "$WDIR/OUTPUT_$INPUT" ]
then
rm -r $WDIR/OUTPUT_$INPUT
fi

mkdir -p $WDIR/OUTPUT_$INPUT

################################ REPLICATES #############################

for r in $(seq 1 $REPS)
do

###################### run SLIM4.3 #########################

num=$(date +%s)$RANDOM

START=$(date +%s)
slim -s $num -t -m -d "Chr=\"3L\"" ./$INPUT
END=$(date +%s)
DIFF=$(( $END - $START ))
echo "Slim4 took 		$DIFF seconds" >> timefile

cp 3Lbott_list_allsnps.txt list_allsnps
sed -i -e 's/3L/1/g' list_allsnps
cp 3Lbott_dataBP.ped dataBP.ped

###################### SIMOVERDOM ########################

echo "$num" > seedfile

START=$(date +%s)
./SIMOVERDOM_ines17<<@
0
-99
51    NINDNP (max 10000)
51     NIND (max 10000)
100             classes
24500000          nwindows
100000          window_size
1              Replicates
@
DIFF=$(( $END - $START ))
echo "SIMOVERDOM took 		$DIFF seconds" >> timefile

cp genfile.dat genfile$r.dat
cp table.txt table$r.txt
cp list_allsnps list_allsnps$r
cp dataBP.ped dataBP$r.ped
cp data.map data$r.map
cp data.ped data$r.ped
cp dfilename.dat dfilename$r.dat
cp frequencies.dat frequencies$r.dat

sleep 5

######################## AVERAGE #########################

START=$(date +%s)
./AVERAGE<<@
$r    NREPS (max 100)
245   NWINDOWS
@
DIFF=$(( $END - $START ))
echo "AVERAGE took 		$DIFF seconds" >> timefile

cp timefile timefile$r

sleep 5

######################## TRANSFER OF FILES TO DIRECTORY #########################

cp -r dataBP$r.ped $WDIR/OUTPUT_$INPUT
cp -r data$r.map $WDIR/OUTPUT_$INPUT
cp -r data$r.ped $WDIR/OUTPUT_$INPUT
cp -r list_allsnps$r $WDIR/OUTPUT_$INPUT
cp -r genfile$r.dat $WDIR/OUTPUT_$INPUT
cp -r table$r.txt $WDIR/OUTPUT_$INPUT
cp -r tableavg.txt $WDIR/OUTPUT_$INPUT
cp -r se.txt $WDIR/OUTPUT_$INPUT
cp -r timefile$r $WDIR/OUTPUT_$INPUT
cp -r tablepi* $WDIR/OUTPUT_$INPUT
cp -r tableDT* $WDIR/OUTPUT_$INPUT
cp -r tablenq* $WDIR/OUTPUT_$INPUT
cp -r tabler2* $WDIR/OUTPUT_$INPUT
cp -r tablec* $WDIR/OUTPUT_$INPUT
cp -r dfilename$r.dat $WDIR/OUTPUT_$INPUT
cp -r frequencies$r.dat $WDIR/OUTPUT_$INPUT
cp -r seedfile$r.dat $WDIR/OUTPUT_$INPUT
cp -r table_C_PI* $WDIR/OUTPUT_$INPUT

done

rm 3Lbott*
rm seedfile
rm data*
rm list*
rm dfilename*
rm frequencies*
rm kk*
rm pi*
rm plink.ld
rm plink.log
rm plink.vcf
rm r2file
rm se.txt
rm table*
rm Tajima*
rm timefile*
rm repfile.dat
rm genfile*



