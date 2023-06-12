#!/bin/bash
#  Bands calculator
#  Quantum Espresso
#  May 11th 2023
#  Acadia Physics

### Section #1 - Setup and Configuration ###
### 1.1 Create and populate Bands directory
cd ../
dirName="$prefix-BandsCalculation"
mkdir $dirName
mkdir ./$dirName/Ephemera
mkdir ./$dirName/GraphFiles
mkdir ./$dirName/OutFiles

### 1.2 check run state and get input file
case $path in
	1)
	cp ./$prefix.tmp/* ./$dirName
	;;
	3)
	inputfile="$prefix.scf.in"
	cp ./QE-Calculator/InputFiles/$inputfile  ./$dirName/
	;;
esac

cp ./$dirName/$prefix.scf.in ./$dirName/$prefix.scfBands.in
cp ./QE-Calculator/Bands.config/* ./$dirName

cd ./$dirName

### 1.3 configure input file

sed -i "s|^calculation = 'scf',|calculation = 'bands',|g" $prefix.scfBands.in

sed -i "s|^K_POINTS automatic|K_POINTS {crystal_b}|g" $prefix.scfBands.in
sed -i '$d' $prefix.scfBands.in

cat $prefix.Bands.append >> $prefix.scfBands.in
sed -i "s|prefix='Si',|prefix='$prefix',|g" BandsConfig.in
sed -i "s|filband='Si.Bands.dat',|filband='$prefix.Bands.dat',|g" BandsConfig.in


### Main body
### 2.1 Readout function
readout(){
xterm -e 'tail -f '$currentfile'' &
sleep 1
}

### 2.2 Bands calculation
Bands(){

echo 'Running pw.x from optimized'
echo '' > $prefix.scf.out
currentfile=$prefix.scf.out
readout
process=$!
sleep 1
pw.x -in $prefix.scf.in > $currentfile
grep -e 'number of electrons' $prefix.scf.out|head -1 > num.electron
electron=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' num.electron)
sed -i "/^ecutwfc =/a nbnd = 8," $prefix.scfBands.in
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
echo '' > BandsConfig.out
currentfile=BandsConfig.out
readout
process=$!
bands.x -in BandsConfig.in > $currentfile
kill $process

### 2.3 Graphing and cleanup
sed -i -r "s|datafile=|datafile='$prefix.Bands.dat.gnu'|g" run.py
sed -i -r "s|symmetryfile=|symmetryfile='BandsConfig.out'|g" run.py
grep -e highest occupied  $prefix.scf.out$i >> $prefix.fermi
grep -e Fermi  $prefix.scf.out$i >> $prefix.fermi
x=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.fermi | tail -1)
echo $x
sed -i -r "s|fermi =|fermi = $x|g" run.py
sed -i -r "s|fig.savefig|fig.savefig('$prefix.band_structure.png')|g" run.py

echo 'plotting'
python3 run.py

mv *.py *.in *.txt *.xml ./Ephemera
mv *.rap *.append *.dat *.gnu ./GraphFiles
mv *.out ./OutFiles
cd ..
}

Bands

	

