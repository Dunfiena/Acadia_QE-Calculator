#!/bin/bash
#  $dirName calculator
#  Quantum Espresso
#  May 14th 2023
#  Acadia Physics

### Section 1 - Setup and Configuration ###

###1.1  Create and populate Kpoint directory

cd ../
dirName="$prefix-BandsCalculation"
mkdir $dirName

cp ./$prefix.tmp/*.in ./$dirName

cd $dirName
grep -P -i '^(\d( \d)+)' $prefix.scf.in > findk.txt
cat > kpoint.in << EOF
2
kpoint.append
EOF
awk '{ print $1 " " $2 " " $3 }' findk.txt >> kpoint.in
awk '{ print $4 " " $5 " " $6 }' findk.txt >> kpoint.in
echo f >>kpoint.in

sleep 2

kpoints.x 0< kpoint.in

cp kpoint.append ../$prefix.tmp
