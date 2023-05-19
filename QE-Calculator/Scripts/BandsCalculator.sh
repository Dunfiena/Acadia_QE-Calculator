#!/bin/bash
#  Bands calculator
#  Quantum Espresso
#  May 11th 2023
#  Acadia Physics

				####################################################
				###						 ###
				###   Section #1 - Setup and Configuration       ###
				###						 ###
				####################################################

### Create and populate Bands directory
readout(){
xterm -e 'tail -f '$currentfile'' &
sleep 1
}

Bands(){
DIR='../Bands/'
if [ -d "$DIR" ]; then
	mkdir ../Bands
	
else
	echo "Bands directory found"	
fi

cd ../
dirName="$prefix-BandsCalculation"
mkdir $dirName
mkdir ./$dirName/Ephemera
mkdir ./$dirName/GraphFiles
mkdir ./$dirName/OutFiles
case $path in
	1)
	cp ./$prefix.tmp/*.in ./$dirName
	;;
	3)
	inputfile="$prefix.scf.in"
	cp ./QE-Calculator/InputFiles/$inputfile  ./$dirName/
	;;
esac

cp ./$dirName/$prefix.scf.in ./$dirName/$prefix.scfBands.in
cp ./QE-Calculator/Bands.config/* ./$dirName

cd ./$dirName
### configure input file
sed -i "s|^calculation = 'scf',|calculation = 'bands',|g" $prefix.scfBands.in
sed -i "/^ecutwfc =/a nbnd = 8," $prefix.scfBands.in
sed -i "s|^K_POINTS automatic|K_POINTS {crystal_b}|g" $prefix.scfBands.in
sed -r -i "s|[0-9]{1} [0-9]{1} [0-9]{1} 1 1 1|11|g" $prefix.scfBands.in
cat $prefix.Bands.append >> $prefix.scfBands.in

echo 'Running pw.x from optimized'
echo '' > $prefix.scf.out
currentfile=$prefix.scf.out
readout
process=$!
sleep 1
pw.x -in $prefix.scf.in > $currentfile
kill $process

echo 'Running pw.x for band'
echo '' > $prefix.scfBands.out
currentfile=$prefix.scfBands.out
readout
process=$!
sleep 1
pw.x -in $prefix.scfBands.in > $currentfile
kill $process

echo 'Running bands.x'
echo '' > $prefix.BandsConfig.out
currentfile=$prefix.BandsConfig.out
readout
process=$!
bands.x -in $prefix.BandsConfig.in > $currentfile
kill $process

echo 'plotting'
python3 run.py

mv *.py *.in *.txt *.xml ./Ephemera
mv *.rap *.append *.dat *.gnu ./GraphFiles
mv *.out ./OutFiles
cd ..
}

Bands

	

