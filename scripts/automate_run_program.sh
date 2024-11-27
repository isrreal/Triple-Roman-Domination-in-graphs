#!/bin/bash

folder_array=(
  "heuristic1 1" 
  "heuristic2 2"
  "heuristic3 3"
  "mixedHeuristics 4"
)

for line in "${folder_array[@]}"
do
    read -r folder heuristic <<< "$line"
    ./run_program.sh -f "$folder" -h "$heuristic"
done












