#!/bin/bash

while getopts "i:o:h:" opt; do
    case $opt in
		i) input_folder_name="$OPTARG" ;;
        h) heuristic="$OPTARG" ;;
        o) output_folder_name="$OPTARG" ;;
    	*) echo "Usage: $0 -i <input_folder_name> -o <output_folder_name> -h <first_number>" ; exit 1 ;;
    esac
done

if [[ -z "$output_folder_name" || -z "$heuristic" ]]; then
    echo "Usage: $0 -i <input_folder_name> -o <output_folder_name> -h <first_number>"
    exit 1
fi

executable="../app"

input_directory="../input_files/$input_folder_name"

output_directory="../output_files/$input_folder_name/$output_folder_name"

mkdir -p "$output_directory"

for input_file in "$input_directory"/*.txt; do
    if [ -f "$input_file" ]; then
        echo "Processing file: $(basename "$input_file")"
		
        output_file="$output_directory/$(basename "$input_file" .txt)_output.csv"
		
        $executable "$input_file" "$(basename "$input_file" .txt)" "$heuristic" > "$output_file"

        echo "Results written to: $output_file"
    else
        echo "No text files found in the directory."
    fi
done

