#!/bin/bash
#  Run calculation suite
#  Quantum Espresso
#  May 17th 2023
#  Acadia Physics


### 1 Check programs and setup ###
prog=('gnuplot', 'gnuplot-x11', 'python3', 'Xterm')

for i in ${prog[@]};
do
if which $i >/dev/null; then
	echo "$i is missing, please install it before proceeding"
	exit
else
	echo "$i located"
fi
done

di='dialog'
if which $di > /usr/share/doc/; then
	echo "$di is missing, please install it before proceeding"
	exit
else
	echo "$di located"
fi
### End of check ###


### 2 Configure run parameters ###
#echo "What is the run state"
#echo "(1) - Run full"
#echo "(2) - Run Self-Consistent Funtional Calculation"
#echo "(3) - Run Electric Bands Calculation"
#echo "(4) - Run Phonon Calculation"
#echo "(5) - Run Dielectric Calculation"

#read path
#export path

#echo "Enter the number of unique atoms (i.e, for SiO2 enter 2)"
#read numAtom
#if (($numAtom > 1)); then
#	for (( i=0; i<$numAtom; i++ ))
#	do
#		echo "Enter the prefix being used (i.e., Si)"
#		read atom$i
#	done
#fi
#echo "$atom0"
#echo "$atom1"

#echo "Please enter prefix"
#read prefix

#if ls ../Pseudo-Potential/$atom0.* && ls ../Pseudo-Potential/$atom1.*; then
#export prefix

d=$(date +"%Y.%m.%d_%H.%M.%S")
dirName="$prefix-test$d"
mkdir ../$dirName

### 3 Run functions ###
RunFull(){
echo "Running SCF Optimization Calculations"
./Scripts/SCFCalculator.sh

echo "Generating Kpoints"
./Scripts/KPointsCalculator.sh

echo "Running Band Calculations"
./Scripts/BandsCalculator.sh

echo "Running Phonon Calculations"
./Scripts/PhononCalculator.sh

echo "Running Dielectric Calculations"
./Scripts/DielectricCalculator.sh


mv ../$prefix* ../$dirName
}


RunSCF(){
echo "Running SCF Optimization Calculations"
./Scripts/SCFCalculator.sh
mv ../$prefix* ../$dirName
}

RunBands(){
echo "Running Band Calculations"
./Scripts/BandsCalculator.sh
mv ../$prefix* ../$dirName
}

RunPhonon(){
echo "Running Phonon Calculations"
./Scripts/PhononCalculator.sh
mv ../$prefix* ../$dirName
}

RunDielectric(){
echo "Running Dielectric Calculations"
./Scripts/DielectricCalculator.sh
mv ../$prefix* ../$dirName
}

case $path in

	1)
	xterm -e top &
	process=$!
	RunFull
	kill $process
	;;

	2)
	xterm -e top &
	process=$!
	RunSCF
	kill $process
	;;
	
	3)
	xterm -e top &
	process=$!
	RunBands
	kill $process
	;;
	
	4)
	xterm -e top &
	process=$!
	RunPhonon
	kill $process
	;;
	5)
	xterm -e top &
	process=$!
	RunDielectric
	kill $process
	;;
esac

else
	echo "File not found"
fi
