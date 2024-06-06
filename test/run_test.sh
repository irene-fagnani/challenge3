#!/bin/bash

# Navigate to the src directory

# Compile the program
cd ..
mpic++ -c -fopenmp -O3 -pedantic -I/home/irene/PACS/pacs-examples/Examples/include  -std=c++20 -Iinclude  src/main.cpp
mpic++ -std=c++20 -L/home/irene/PACS/pacs-examples/Examples/lib main.o -lmuparser -lgomp -o main
# Navigate back to the test directory
# cd test

# Run the program with 1, 2, and 4 cores
for cores in 1 2 4; do
    for i in 4 5 6 7 8; do
     k=$((2**$i))
     echo "Running with $cores core(s) and $k dimension of the grid"
     mpirun -np $cores ./main $k test=1
    done
done