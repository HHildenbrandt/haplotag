#!/bin/bash
output=$1
g++ -O3 -o ${output/.cpp/.o} $output -lgzstream -I/scratch/hb-1000G_str/local/gzstream -L/scratch/hb-1000G_str/local/gzstream -I/usr/include -std=gnu++11 -lz -Wall
