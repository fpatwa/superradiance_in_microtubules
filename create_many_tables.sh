#!/bin/bash

for n in 40 50 60 70 80 90;
do
	python create_spiral.py $n
	fname="trp_positions_dipoles_${n}_spiral.txt"
	python create_table.py $fname
	./reset_files.sh
done	
