#!/bin/bash
#  Dielectric calculator
#  Quantum Espresso
#  May 11th 2023
#  Acadia Physics

				####################################################
				###						 ###
				###   Section #1 - Setup and Configuration       ###
				###						 ###
				####################################################

### Create and populate Dielectric directory
readout(){
xterm -e "tail -f '$currentfile'" &
}

DIR='../Dielectric/'
if [ -d "$DIR" ]; then
	mkdir ../Dielectric
	
else
	echo "Dielectric directory found"	
fi

cd ../
dirName="$prefix-DielectricCalculation"
mkdir $dirName
mkdir ./$dirName/Ephemera
mkdir ./$dirName/GraphFiles
mkdir ./$dirName/OutFiles

case $path in
	1)
	cp ./$prefix.tmp/* ./$dirName
	;;
	5)
	inputfile="$prefix.scf.in"
	cp ./QE-Calculator/InputFiles/$inputfile  ./$dirName/
	;;
esac

cp ./$dirName/$prefix.scf.in ./$dirName/$prefix.scfDielectric.in
cp ./QE-Calculator/Dielectric.config/* ./$dirName

cd ./$dirName

sed -i "s|^calculation = 'scf',|calculation = 'nscf',|g" $prefix.scfDielectric.in
sed -i "/^ecutwfc =/a nbnd = 8," $prefix.scfDielectric.in
sed -i "/^nbnd =/a nosym = .true.," $prefix.scfDielectric.in
sed -i "/^nbnd =/a nosym = .true.," $prefix.scfDielectric.in

sed -i "s/  prefix = 'Si',/  prefix = '$prefix',/g" eps.in

sed -i "s|plot \"eels_Si.dat\" using 1:2 w lines|plot \"eels_$prefix.dat\" using 1:2 w lines|g" EelsvsEV.gnu
sed -i "s|plot \"epsi_Si.dat\" using 1:2 w lines, \"epsr_Si.dat\" using 1:2 w lines|plot \"epsi_$prefix.dat\" using 1:2 w lines, \"epsr_$prefix.dat\" using 1:2 w lines|g" EpsvsEV.gnu

echo 'Running pw.x from optimized'
echo '' > $prefix.scf.out
currentfile=$prefix.scf.out
readout
process=$!
sleep 1
pw.x -in $prefix.scf.in > $currentfile
kill $process

echo 'Running pw.x for Dielectric'
echo '' > $prefix.scfDielectric.out
currentfile=$prefix.scfDielectric.out
readout
process=$!
sleep 1
pw.x -in $prefix.scfDielectric.in > $currentfile
kill $process

echo "Running epsilon.x"
epsilon.x -in eps.in

### Plotting graphs
gnuplot < EpsvsEV.gnu
gnuplot < EelsvsEV.gnu

mv *.dat *.gnu ./GraphFiles
mv *.in *.xml *.wfc1 ./Ephemera
mv *.out ./OutFiles



