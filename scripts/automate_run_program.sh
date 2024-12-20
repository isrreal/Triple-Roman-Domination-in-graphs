#!/bin/bash

output_folder_array=(
  "heuristic1 1" 
  "heuristic2 2"
  "heuristic3 3"
  "mixedHeuristics 4"
)

for line in "${output_folder_array[@]}"
do
    read -r folder heuristic <<< "$line"
    ./run_program.sh -i random_graphs -o "$folder" -h "$heuristic"
done












