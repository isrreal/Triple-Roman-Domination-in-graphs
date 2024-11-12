#!/bin/bash

# C++ executable
executable="./app"

# Directory containing the text files
input_directory="./input_files"

# Directory to store the output files
output_directory="./output_files"

# Create the output directory if it doesn't exist
mkdir -p "$output_directory"

# Loop through each text file in the input directory
for input_file in "$input_directory"/*.txt; do
    # Check if the file exists
    if [ -f "$input_file" ]; then
        echo "Processing file: $input_file"

        # Define the output file name (e.g., output_files/file1_output.csv)
        output_file="$output_directory/$(basename "$input_file" .txt)_output.csv"

        # Run the C++ executable and redirect its output to the corresponding output file
        $executable "$input_file" "$(basename "$input_file" .txt)" > "$output_file"

        echo "Results written to: $output_file"
    else
        echo "No text files found in the directory."
    fi
done
