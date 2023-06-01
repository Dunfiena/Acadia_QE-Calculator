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
	cp ./$prefix.tmp/*.in ./$dirName
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
ph.x -in $prefix.ph.in > $prefix.ph.out
kill $process

echo "Running q2r"
q2r.x -in q2r.in > q2r.out

echo '' > $prefix.phdos
currentfile=$prefix.phdos
readout
process=$!
echo "Running matdyn"
matdyn.x -in $prefix.matdyn-Uniform.in > $prefix.matdyn-Uniform.out
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
