#!/bin/bash

# Define the Python script and output CSV file
PYTHON_SCRIPT="pyomo2.py"

# Run the Python script and save the output to the CSV file
python3 "$PYTHON_SCRIPT"

echo "Python script executed, and output saved to $OUTPUT_CSV"
