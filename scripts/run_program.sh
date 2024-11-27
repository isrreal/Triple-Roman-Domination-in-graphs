#!/bin/bash

while getopts "f:h:" opt; do
    case $opt in
        f) folder_name="$OPTARG" ;;
        h) heuristic="$OPTARG" ;;
        *) echo "Usage: $0 -f <folder_name> -h <GA_heuristic>" ; exit 1 ;;
    esac
done

if [[ -z "$folder_name" || -z "$heuristic" ]]; then
    echo "Usage: $0 -f <folder_name> -h <first_number>"
    exit 1
fi

executable="../app"

input_directory="../input_files/representative_set"

output_directory="../output_files/genetic_algorithm_tests"

tests_directory="$output_directory/$folder_name"
mkdir -p "$tests_directory"

for input_file in "$input_directory"/*.txt; do
    if [ -f "$input_file" ]; then
    	# basename gets the name of the file, ignores relative path
        echo "Processing file: $(basename "$input_file")"
		
    	# basename, in this case, removes .txt extension in the file $input_file
        output_file="$tests_directory/$(basename "$input_file" .txt)_output.csv"
		
        $executable "$input_file" "$(basename "$input_file" .txt)" "$heuristic" > "$output_file"

        echo "Results written to: $output_file"
    else
        echo "No text files found in the directory."
    fi
done

