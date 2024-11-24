#!/bin/bash

folder_array=(
  "heuristic1_without_RVNS 1 0"
  "heuristic1_with_RVNS 1 1"  
  "heuristic2_without_RVNS 2 0"
  "heuristic2_with_RVNS 2 1"
  "heuristic3_without_RVNS 3 0"
  "heuristic3_with_RVNS 3 1"
)

for line in "${folder_array[@]}"
do
    read -r folder heuristic rvns <<< "$line"
    ./run_program.sh -f "$folder" -h "$heuristic" -R "$rvns"
done












