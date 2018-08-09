#!/bin/bash

for i in `seq 400 50 1200`;
do

	echo "############## Radius: " $i "px ################"
	python Stars-Exp-Feb2017.py $i
done 
