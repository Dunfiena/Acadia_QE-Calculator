#!/bin/bash
# setup.sh is the gui componant of QE-Calculator
INPUT=/tmp/menu.sh.$$
OUTPUT=/tmp/output.sh.$$


# trap and delete temp files
trap "rm $OUTPUT; rm $INPUT; exit" SIGHUP SIGINT SIGTERM

### Variables for runstate
Fileopt="Unselected"
Runopt="All"
Graphopt="All"
Readyopt="Please provide missing information"

#
# Purpose - select input file
#
function Select_Input(){
	File=$(dialog --title "Please select you input file" --stdout --title "Please choose a file to delete" --fselect InputFiles/ 14 48)
	prefix=${File:11:2}
	export prefix
}


function Select_Run(){
	dialog --title "Select run type" --menu "Select the run you want" 30 40 5 \
	RunAll "Run all scripts" \
	SCF "Run the SCF" \
	Bands "Run the Bands" \
	Phonon "Run the Phonon" \
	Dielectric "Run the Dielectric" 2>"${INPUT}"
	
	item=$(<"${INPUT}")
	case $item in
		RunAll) path=1;;
		SCF) path=2;;
		Bands) path=3;;
		Phonon) path=4;;
		Dielectric) path=5;;
	esac
	export path
	}

function Select_Graph(){
	dialog --title "Select graphs" --menu "Select all the graphs that you want to run \n\
	If all is selected it will proceed to the next step." 30 40 2 \
	All "Generate all Graphs"\
	Select "Select Generated graphs" 2>"${INPUT}"
	
	item2=$(<"${INPUT}")	
	case $item2 in
		All) Graph=0;;
		Select) Graphs;;
	esac
}

function Graphs(){
dialog --title "Select graphs" --checklist "Select all the graphs that you want to run \n\
	You can select as many as you want" 30 40 4 \
	SCF "SCF Graphs" off\
	Bands "Bands Graphs" off\
	Phonon "Phonon Graphs" off\
	Dielectric "Dielectric Graphs" off 2>../$prefix.graphsSel.txt

	}

Select_Input
Select_Run
Select_Graph

./Scripts/run.sh
echo $File
echo $path
echo $prefix
echo $Graph
echo $graphsSel






