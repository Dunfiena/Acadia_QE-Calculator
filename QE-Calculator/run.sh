#!/bin/bash
#  Run calculation suite
#  Quantum Espresso
#  May 17th 2023
#  Acadia Physics


### 1.1 Check programs and set up to ensure working system ###
GPlot=`which gnuplot 2> /dev/null`
if [ "$GPlot" = "" ]; then
	echo "Gnuplot is missing, if you want to plot the data output please install gnuplot using sudo apt install gnuplot"
else
	echo "Gnuplot located"
fi

### 1.2 Configure run parameters ###
echo "What is the run state"
echo "(1) - Run full"
echo "(2) - Run Self-Consistent Funtional Calculation"
echo "(3) - Run Electric Bands Calculation"
echo "(4) - Run Phonon Calculation"
echo "(5) - Run Dielectric Calculation"

read path
export path

echo "Enter the prefix being used (i.e., Si)"
read prefix
export prefix

d=$(date +"%Y.%m.%d_%H.%M.%S")
dirName="$prefix-test$d"
mkdir ../$dirName

RunFull(){
echo "Running SCF Optimization Calculations"
./Scripts/SCFCalculator.sh

echo "Running Band Calculations"
./Scripts/BandsCalculator.sh

echo "Running Phonon Calculations"
./Scripts/PhononCalculator.sh

echo "Running Dielectric Calculations"
./Scripts/DielectricCalculator.sh


mv ../$prefix* ../$dirName
}

### 1.3 Individual runs ###
RunSCF(){
echo "Running SCF Optimization Calculations"
./Scripts/SCFCalculator.sh
}

RunBands(){
echo "Running Band Calculations"
./Scripts/BandsCalculator.sh
}

RunPhonon(){
echo "Running Phonon Calculations"
./Scripts/PhononCalculator.sh
}

RunDielectric(){
echo "Running Dielectric Calculations"
./Scripts/DielectricCalculator.sh
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
