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

# Loop sobre as pastas de entrada
for input in "${input_folder_array[@]}"; do

  # Executa o programa com os parâmetros
  ./run_program_aco.sh -i "$input" -o "$input" -h 0
done
