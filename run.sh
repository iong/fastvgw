#!/bin/bash

export OMP_NUM_THREADS=8
for LAMBDA in 0.09; do
#for LAMBDA in 0.03 0.06 0.09 0.12 0.15 0.19 0.21; do
	rm mdp_deBoer=${LAMBDA}.dat
	for f in `ls *.mdp?.xyz | sort -n`; do
		echo "deBoer = $LAMBDA $f" 1>&2
		./clustergs $f $LAMBDA >> mdp_deBoer=${LAMBDA}.dat
	done
	rm ticos_deBoer=${LAMBDA}.dat
	for f in `ls *.ticos.xyz| sort -n`; do
		echo "deBoer = $LAMBDA $f" 1>&2
		./clustergs $f $LAMBDA >> ticos_deBoer=${LAMBDA}.dat
	done
done
