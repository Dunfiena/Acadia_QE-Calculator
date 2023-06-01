#!/bin/bash
#  Recursive running for ecutwfc iterative and K points
#  Quantum Espresso
#  May 5th 2023
#  Acadia Physics


### 1.1 Setup test directory and files for program ###

cd ../
mkdir "$prefix.tmp"

echo "Test running directory"
dirName="$prefix-scfCalculation"
inputfile=$prefix.scf.in
mkdir $dirName
mkdir $dirName/EnergyOut
mkdir $dirName/KpointOut
mkdir $dirName/CellOut
mkdir $dirName/Graphs
mkdir $dirName/Totals
mkdir $dirName/Graphs/Datafiles
pwd
cp ./QE-Calculator/InputFiles/$inputfile  ./$dirName/
cp ./QE-Calculator/SCF.config/* ./$dirName/Graphs/

cd ./$dirName/
echo -e "1  2  3\n2  3  4" >> ./Graphs/energy.dat
echo -e "1  2  3\n2  3  4" >> ./Graphs/Kpoint.dat
echo -e "1  2  3\n2  3  4"> ./Graphs/Datafiles/celldmtmp.dat

### 1.2 Clean input files before run ###
echo "Clean input file"
sed -i "s|^ecutwfc =.*,|ecutwfc = 40,|g" $inputfile

### 1.3 Recursive array creation ###

# Adjustable section
estop=0.003
kstop=0.003

energymin=10
energymax=100
energystep=10

Kpointmin=2
Kpointmax=32
Kpointstep=2
# end of adjustable

celldm=()
Ecut=("$energymin")
Kpoints=("$Kpointmin")

cellLine=$(grep -e 'celldm(1) = *' $inputfile)
cell=$(echo $cellLine | grep -Eo '[0-9]{1,9}.[0-9]{1,9}')
scaling=("-20" "-15" "-10" "-5" "0" "5" "10" "15" "20")

for i in ${scaling[@]}
do
	cell1=$(bc <<< "scale=6; ((($i/100)*$cell)+$cell)")
	celldm+=("$cell1")
done

while [ $energymin -lt $energymax ]
do
	energymin=$(($energymin + $energystep))
	Ecut+=("$energymin")
done

while [ $Kpointmin -lt $Kpointmax ]
do
	Kpointmin=$(($Kpointmin + $Kpointstep))
	Kpoints+=("$Kpointmin")
done

### 1.4 Declare global variables ###
b=0
e1=${Ecut[0]}
k1=${Kpoints[0]}
c1=${celldm[0]}


### Section 2 - For multiplot info please see liveplot.gnu ###

### Section 3 - Main Body

Run(){

#####
### 3.1 Recursive loop for energy cutoff ###
#####

### 3.1.1 Running calculations
for i in "${Ecut[@]}"
do
	echo "Running for $i RY"
	sed -i "s|^ecutwfc =.*,|ecutwfc = $i,|g" $inputfile
	#mpirun -np 12 pw.x -in $inputfile > $prefix.scf.out$i	
	pw.x -in $inputfile > $prefix.scf.out$i
	
	### 3.1.2 Write data to files###	
	grep -e ^! $prefix.scf.out$i > $prefix.scf.total
	grep -e estimate $prefix.scf.out$i |tail -1 > $prefix.scf.accuracy
	grep -e PWSCF $prefix.scf.out$i |tail -1 > $prefix.scf.runtime

	#Write to Graph files#
	x=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.total | tail -1) 
	y=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.runtime| tail -1)
	clean=$(grep -Eo '1  2  3' ./Graphs/energy.dat)
	if [ "$clean" = "1  2  3" ]; then
		echo -n > ./Graphs/energy.dat
	fi
	echo "$i            -$x         $y" >> ./Graphs/energy.dat 
	a=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.total | tail -1)

	### 3.1.3 Calculte the energy change by percentage
	# 'a' is the energy o'f the current calculation
	# and 'b' is the calculation of the previous iteration
	if [ $i -gt $e1 ]; then
		c=$(bc<<< "scale=9;(($a-$b)/(($a+$b)/2))*100")
		echo "Percentage change is $c"
		sleep 1
		###Check the percent change of the previous calculation against cutoff points
		###'c' is percent < 'd' is the cutoff.  If cut off reached, go to K_Point caluculation
		if (( $(echo "$c  < $estop"|bc -l) )); then
			break
		fi
	fi

	sleep 1
	###Assign 'b' to energy of current run so that it can be called next iteration
	b=$a
done

#####
### 3.2 Recursive loop for Kpoints calculation
#####

### 3.2.1 Running calculations
for i in "${Kpoints[@]}"
do
	echo "Running for $i KPoints"

	sed -r -i "s|[1-9]+\s[1-9]+\s[1-9]+\s|$i $i $i |g" $inputfile
	sed -r -i "s|[1-9]+\s[1-9]+\s[1-9]+\s|$i $i $i |g" $inputfile
	#mpirun -np 12 pw.x -in $inputfile > $prefix.scf.Kout$i
	pw.x -in $inputfile > $prefix.scf.Kout$i
	
	### 3.2.2 Write data to files###
	grep -e ^! $prefix.scf.Kout$i >> $prefix.scf.Ktotal
	grep -e estimate $prefix.scf.Kout$i |tail -1 >> $prefix.scf.Kaccuracy
	grep -e PWSCF $prefix.scf.Kout$i |tail -1 >> $prefix.scf.Kruntime

	#Write K_Point Data to Graph#
	x=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.Ktotal | tail -1) 
	y=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.Kruntime| tail -1)
	clean=$(grep -Eo '1  2  3' ./Graphs/Kpoint.dat)
	if [ "$clean" = "1  2  3" ]; then
		echo -n > ./Graphs/Kpoint.dat
	fi
	echo "$i            -$x         $y" >>./Graphs/Kpoint.dat 

	a=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.Ktotal | tail -1)

	### 3.2.3 Calculate the energy change by percentage
	# Where 'a' is the energy of the current calculation
	# and 'b' is the calculation of the previous iteration
	if [ $i -gt $k1 ]; then
		c=$(bc<<< "scale=9;(($a-$b)/(($a+$b)/2))*100")
		echo "Percentage change is $c"
		sleep 1
		###Check the percent change of the previous calculation against cutoff points
		###'c' is percent < 'd' is the cutoff.  If cut off reached, go to K_Point caluculation
		if (( $(echo "$c  < $kstop"|bc -l) )); then
		break
		fi
	fi
	sleep 1
	###Assign 'b' to energy of current run so that it can be called next iteration
	b=$a
done


#####
### 3.3 Ev.x - Calculating volumn in Angstrom
####

### 3.3.1 Running Calculations
for i in "${celldm[@]}"
do
	echo "Running for $i Cell"
	sed -i "s|^celldm(1) =.*|celldm(1) = $i,|g" $inputfile
	#mpirun -np 12 pw.x -in $inputfile > $prefix.scf.cellout$i
	pw.x -in $inputfile > $prefix.scf.cellout$i
	
	### 3.3.2 Write data to files
	grep -e ^! $prefix.scf.cellout$i >> $prefix.scf.celltotal
	grep -e PWSCF $prefix.scf.cellout$i |tail -1 >> $prefix.scf.cellruntime

	x=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.celltotal | tail -1)
	y=$(bc <<< "scale=6; (($i*0.529177))")
	z=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.cellruntime| tail -1)
	echo "$y            -$x" >>./Graphs/celldm.dat
	a=$(grep -Eo '[0-9]{1,9}.[0-9]{1,9}' $prefix.scf.celltotal | tail -1)
	
	clean=$(grep -Eo '1  2  3' ./Graphs/Datafiles/celldmtmp.dat)
	if [ "$clean" = "1  2  3" ]; then
		echo -n > ./Graphs/Datafiles/celldmtmp.dat
	fi

	echo "$y        -$x        $z" >>./Graphs/Datafiles/celldmtmp.dat
	
	### 3.3.3 Calculate the energy change by percentage
	# Where 'a' is the energy of the current calculation
	# and 'b' is the calculation of the previous iteration
	if (( $(echo "$i  > $c1"|bc -l) )); then
		if (( $(echo "$a > $b"|bc -l) )); then
		cp $prefix.scf.in ../$prefix.tmp
		fi
	fi
	sleep 1
	###Assign 'b' to energy of current run so that it can be called next iteration
	b=$a
done


### 3.3.4 Ev.x setup and run
cd ./Graphs
ev.x < evxcalc.tmp


### Section #4 - Final Graph creation ###

#####
###Gnuplot for dat files
#####
echo "Plotting Graphs..."
### Ry Cutoff vs Total energy
gnuplot < EnergyvsCut.gnu

### Ry Cutoff vs Runtime
gnuplot < RunvsCut.gnu

### Kpoints vs Total Energy
gnuplot < EnergyvsKpoints.gnu

### KPoints vs Runtime
gnuplot < RunvsKpoints.gnu

### Total energy vs Volume ev.x Graph
gnuplot < EnergyvsLattice.gnu

kill $processID
}

### Section #5 - Execute and clean ###
gnuplot ./Graphs/liveplot.gnu &
processID=$!
sleep 1
Run


### Clean up ###
mv *.dat *.tmp *.gnu ./Datafiles

cd ../

mv $prefix.scf.cellout* ./CellOut/
mv $prefix.scf.out* ./EnergyOut/
mv $prefix.scf.Kout* ./KpointOut/
mv $prefix.scf.* *.xml ./Totals/
cd ../

