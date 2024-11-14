#!/bin/bash

while getopts "f:h:R:" opt; do
    case $opt in
        f) folder_name="$OPTARG" ;;
        h) first_number="$OPTARG" ;;
        R) second_number="$OPTARG" ;;
        *) echo "Usage: $0 -f <folder_name> -h <GA_heuristic> -R <GA_with_RNVS_or_no>"; exit 1 ;;
    esac
done

# -z verify if the string is empty
if [[ -z "$folder_name" || -z "$first_number" || -z "$second_number" ]]; then
    echo "Usage: $0 -f <folder_name> -h <first_number> -r <second_number>"
    exit 1
fi

executable="./app"

input_directory="./representative_set"

output_directory="./output_files"

tests_directory="$output_directory/$folder_name"
mkdir -p "$tests_directory"

for input_file in "$input_directory"/*.txt; do
    if [ -f "$input_file" ]; then
    	# basename gets the name of the file, ignores relative path
        echo "Processing file: $(basename "$input_file")"
		
    	# basename, in this case, removes .txt extension in the file $input_file
        output_file="$tests_directory/$(basename "$input_file" .txt)_output.csv"
		
        $executable "$input_file" "$(basename "$input_file" .txt)" "$first_number" "$second_number" > "$output_file"

        echo "Results written to: $output_file"
    else
        echo "No text files found in the directory."
    fi
done

