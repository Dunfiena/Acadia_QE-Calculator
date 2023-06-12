#!/bin/bash
#  $dirName calculator
#  Quantum Espresso
#  May 14th 2023
#  Acadia Physics

### Section 1 - Setup and Configuration ###

###1.1  Create and populate Phonon directory

cd ../
dirName="$prefix-PhononCalculation"
mkdir $dirName
mkdir ./$dirName/Ephemera
mkdir ./$dirName/GraphFiles
mkdir ./$dirName/OutFiles
mkdir ./$dirName/DynFiles

### 1.2 check run state and get input file
case $path in
	1)
	cp ./$prefix.tmp/* ./$dirName
	;;
	4)
	inputfile="$prefix.scf.in"
	cp ./QE-Calculator/InputFiles/$inputfile  ./$dirName/
	;;
esac

cp ./QE-Calculator/Phonon.config/* ./$dirName
cd ./$dirName

### 1.3 configure input file
sed -i "s/Si.freqNU.gp/$prefix.freqNU.gp/g" DispersionRelation.plot.gnu

sed -i "s/Density of state of Si Crystal/Density of state of $prefix Crystal/g" DOS.plot.gnu
sed -i "s/Si.phdos/$prefix.phdos/g" DOS.plot.gnu

sed -i "s/Si.fc/$prefix.fc/g" q2r.in
grep -e ^$prefix $prefix.scf.in > automicmass.txt
amass=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' automicmass.txt|head -1)
echo $amass
sed -i "s/  prefix = 'Si',/  prefix = '$prefix',/g" ph.in
sed -i "s/  amass(1) = 28.086,/  amass(1) = $amass,/g" ph.in

sed -i "s/  amass(1) = 28.086,/  amass(1) = $amass,/g" matdyn-Uniform.in
sed -i "s/  flfrc = 'Si.fc',/  flfrc = '$prefix.fc',/g" matdyn-Uniform.in
sed -i "s/  flfrq = 'Si.freq',/  flfrq = '$prefix.freq',/g" matdyn-Uniform.in
sed -i "s/  fldos = 'Si.phdos',/  fldos = '$prefix.phdos',/g" matdyn-Uniform.in

sed -i "s/  amass(1) = 28.086,/  amass(1) = $amass,/g" $prefix.matdyn-nonUniform.in
sed -i "s/  flfrc = 'Si.fc',/  flfrc = '$prefix.fc',/g" $prefix.matdyn-nonUniform.in
sed -i "s/  flfrq = 'Si.freqNU',/  flfrq = '$prefix.freqNU',/g" $prefix.matdyn-nonUniform.in
#awk '{ print $1 }' kpoint.append|head -1 >> matdyn-nonUniform.in
#sed -i 's/\.00$/ /g' kpoint.append
#awk '{ print $2 "  "  $3 "  " $4 "  "  $5 }' kpoint.append > kpoint
#cat kpoint >> matdyn-nonUniform.in

### Main body
### 2.1 Readout function
readout(){
xterm -e "tail -f '$currentfile'" &
}

Phonon(){
### 2.2 Phonon calculation
echo "Running pw.x"
pw.x -in $prefix.scf.in > $prefix.scf.out

echo '' > $prefix.ph.out
currentfile=$prefix.ph.out
readout
process=$!
echo "Running ph.x"
ph.x -in ph.in > $prefix.ph.out
kill $process

echo "Running q2r"
q2r.x -in q2r.in > q2r.out

echo '' > $prefix.phdos
currentfile=$prefix.phdos
readout
process=$!
echo "Running matdyn"
matdyn.x -in matdyn-Uniform.in > $prefix.matdyn-Uniform.out
echo "Running maytdyn nonUniform"
matdyn.x -in $prefix.matdyn-nonUniform.in > $prefix.matdyn-nonUniform.out
kill $process

### 2.3 Graph and cleanup
gnuplot DOS.plot.gnu
gnuplot DispersionRelation.plot.gnu

mv *.in *.xml *.wfc1  ./Ephemera
mv *.gnu *.phdos *.gp *.freq *.freqNU *.fc ./GraphFiles
mv *.out *.modes ./OutFiles
mv dyn* ./DynFiles

cd ..
}
Phonon
