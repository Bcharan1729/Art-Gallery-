#!/bin/bash

# Check if input file is provided
if [ $# -eq 0 ]; then
    echo "Usage: ./test.sh input.txt"
    exit 1
fi

input_file="$1"

# Check if file exists
if [ ! -f "$input_file" ]; then
    echo "Error: File '$input_file' not found"
    exit 1
fi

if [ ! -f "gallery_triangulation.py" ]; then
    echo "test.sh and required code is not in the same folder / the code is missing "
    exit 1    
fi

# Run the Python script with input redirection
python3 gallery_triangulation.py < "$input_file"