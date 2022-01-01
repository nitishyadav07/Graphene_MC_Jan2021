#!/bin/bash
#! /usr/bin/env python
#
#
#echo Diamagnetic subtraction initiated.
#echo
#for index in {0..3}; do
#for index in $(seq 0 1 3); do
#	python3 Diamag_Component.py --i $index
#	echo
#done
#echo Diamagnetic subtraction done.
#echo ------------------------------
#echo
#
#
echo MC simulation begun.
#python3 Ising_sim_1D.py --J1 0.0 --J2 0.0 --N 80
#python3 Ising_sim_1D.py --J1 0.2 --J2 0.001 --N 400
#python3 Ising_sim_1D.py --J1 0.2 --J2 0.001 --N 1000
#python3 Ising_sim_1D.py --J1 0.2 --J2 -0.001 --N 200
#python3 Ising_sim_1D.py --J1 0.2 --J2 -0.001 --N 500
python3 Ising_sim_1D.py --J1 0.0 --J2 0.0 --N 20
echo MC simulation complete.
echo ------------------------------
echo
echo
#echo
#
#
echo Plotting process initiated.
echo
for index in {0..3}; do
	python3 MH_add_and_plot.py --i $index
	echo
done
echo
echo Plotting process complete.
echo ------------------------------
echo


#for i in $(seq 0.0 0.1 1); do
#for i in $(seq 0.5 0.1 1); do						# intraedge ferromagnetic strong coupling
#for i in $(seq 2 0.2 3); do						# intraedge superparamag case
#for i in $(seq 0.6 0.2 3); do						# intraedge superparamag case
#	for j in $(seq 0.0 0.001 0.005); do
#	for j in $(seq -0.0001 -0.0001 -0.0005); do		# interedge antiferromagnetic coupling 		(for ferrimag case)
#	for j in $(seq -0.001 -0.001 -0.005); do		# interedge antiferromagnetic coupling 		(for ferrimag case)
#	for j in $(seq -0.0005 0.0001 0.0005); do		# interedge ferromagnetic weak coupling 	(for ferromag case)
#	for j in $(seq 0.001 0.001 0.005); do			# interedge ferromagnetic weak coupling 	(for ferromag case)
#	for j in $(seq 0 1 0); do						# for superparamag case
#		python3 trans.py $i $j

#		python3 MH_add_and_plot.py --J1 $i --J2 $j --N 100
#		python3 MH_add_and_plot.py --J1 $i --J2 $j --N 50
#		python3 MH_add_and_plot.py --J1 $i --J2 $j --N 40

#		python3 Ising_sim_1D.py --J1 $i --J2 $j --N 1000
#		python3 Ising_sim_1D.py --J1 $i --J2 $j --N 400
#		python3 Ising_sim_1D.py --J1 0.2 --J2 0.001 --N 1000
#		python3 Ising_sim_1D.py --J1 0.2 --J2 0.001 --N 400

#		python3 Ising_sim_1D.py --J1 $i --J2 0.0 --N 100
#		python3 Ising_sim_1D.py --J1 $i --J2 0.0 --N 50
#		python3 Ising_sim_1D.py --J1 $i --J2 0.0 --N 40
#		python3 Ising_sim_1D.py --J1 $i --J2 0.0 --N 20


#		python3 Ising_sim_1D.py --J1 0.0 --J2 0.0 --N 8			# best SPM, with g = 16*2.002

#	done						
# done

#for i in $(seq 0.0 0.1 1); do
#	for j in $(seq -0.001 -0.001 -0.005); do
#		python3 Ising_sim_1D.py --J1 $i --J2 $j --N 200
#		python3 Ising_sim_1D.py --J1 $i --J2 $j --N 400
#		python3 Ising_sim_1D.py --J1 $i --J2 $j --N 500
#	done						
# done

