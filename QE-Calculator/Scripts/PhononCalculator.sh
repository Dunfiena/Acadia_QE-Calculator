#!/bin/bash
#  $dirName calculator
#  Quantum Espresso
#  May 14th 2023
#  Acadia Physics

				####################################################
				###						 ###
				###   Section #1 - Setup and Configuration       ###
				###						 ###
				####################################################

### Create and populate Phonon directory
readout(){
xterm -e "tail -f '$currentfile'" &
}
Phonon(){
cd ..
dirName="$prefix-PhononCalculation"
mkdir $dirName
mkdir ./$dirName/Ephemera
mkdir ./$dirName/GraphFiles
mkdir ./$dirName/OutFiles
mkdir ./$dirName/DynFiles

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
echo "Running pw.x"
pw.x -in $prefix.scf.in > $prefix.scf.out

echo '' > $prefix.ph.out
currentfile=$prefix.ph.out
readout
process=$!
echo "Running ph.x"
ph.x -in $prefix.ph.in > $prefix.ph.out
kill $process

echo '' > $prefix.q2r.out
currentfile=$prefix.q2r.out
readout
process=$!
echo "Running q2r"
q2r.x -in $prefix.q2r.in > $prefix.q2r.out
kill $process

echo '' > $prefix.phdos
currentfile=$prefix.phdos
readout
process=$!
echo "Running matdyn"
matdyn.x -in $prefix.matdyn-Uniform.in
matdyn.x -in $prefix.matdyn-nonUniform.in
kill $process

gnuplot DOS.plot.gnu
gnuplot DispersionRelation.plot.gnu

mv *.in *.xml *.wfc1  ./Ephemera
mv *.gnu *.phdos *.gp *.freq *.freqNU *.fc ./GraphFiles
mv *.out *.modes ./OutFiles
mv dyn* ./DynFiles

cd ..
}
Phonon
