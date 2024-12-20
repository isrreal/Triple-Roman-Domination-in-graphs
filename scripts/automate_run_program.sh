#!/bin/bash

while getopts "f:" opt; do
	case $opt in
		f) input_folder_name="$OPTARG" ;;
	esac
done

output_folder_array=(
  "heuristic1 1" 
  "heuristic2 2"
  "heuristic3 3"
  "mixedHeuristics 4"
)

for line in "${output_folder_array[@]}"
do
    read -r folder heuristic <<< "$line"
    ./run_program.sh -i "$input_folder_name" -o "$folder" -h "$heuristic"
done












