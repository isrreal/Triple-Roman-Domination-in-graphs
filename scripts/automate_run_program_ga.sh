#!/bin/bash

# Array de pastas de entrada
input_folder_array=(
  "BAI" 
  "Cycles" 
  "Harwell-Boeing" 
  "random_graphs"
  "Stars"
  "Trees"
)

# Array de pastas de saída e parâmetros
output_folder_array=(
  "heuristic1 1 0" 
  "heuristic2 2 0"
  "heuristic3 3 0"
  "mixedHeuristics 4 0"
)

# Loop sobre as pastas de entrada
for input in "${input_folder_array[@]}"; do
  # Loop sobre as pastas de saída e parâmetros
  for line in "${output_folder_array[@]}"; do
    # Divide a linha em variáveis
    read -r folder heuristic flag_aco <<< "$line"
    
    # Executa o programa com os parâmetros
    ./run_program_ga.sh -i "$input" -o "$folder" -h "$heuristic"
  done
done
