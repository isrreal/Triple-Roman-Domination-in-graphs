#!/bin/bash

while getopts "f:i:o:" opt; do
    case $opt in
		f) input_folder_name="$OPTARG" ;;
		i) input_file_name="$OPTARG" ;;
		o) output_folder_name="$OPTARG" ;;
   esac
done

executable="../app"

#input_directory="../input_files/$input_folder_name"

output_directory="../output_files/genetic_algorithm_outputs/GA-Two-Point-Crossover/BAI/mixedHeuristics"

mkdir -p "$output_directory"

#for input_file in "$input_directory"/*.txt; do
	$input_file_name="../input_files/BAI/bwm200.txt"
   	if [ -f "$input_file_name" ]; then

        echo "Processing file: $(basename "$input_file_name")"
		
        output_file="$output_directory/$(basename "$input_file_name" .txt)_output.csv"
		
        $executable "$input_file_name" "$(basename "$input_file_name" .txt)" 4 0 0 > "$output_file"

        echo "Results written to: $output_file"
    else
        echo "No text files found in the directory."
    fi
#done

