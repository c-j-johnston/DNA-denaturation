"""
input_file_create.py - Cameron Johnston
This program will create a LAMMPS input file, 'lammps_input', describing
a polymer with randomly generated initial positions for atoms. A single
or double chain can be created depending on the paramters given in the
input for this script.
"""

import random
import math as m
import numpy as np
import sys

# Opens input file, read data and sets parameters
params = open(sys.argv[1], 'r')
lines = params.readlines()
for i in range(len(lines)):
	lines[i] = lines[i].strip()
	lines[i] = lines[i].split(" ")
	lines[i] = np.array(lines[i])

# Number of atoms per chain
n = int((lines[0])[1])
# Number of chains (1 or 2)
chains = int((lines[1])[1])
# This is half the box size, it goes from -(box_size) -> box_size
box_size = float((lines[2])[1]) 
# Number of types of atoms. For homogeneous chains, this should be
# same a number of chains to allow VMD to label chains differently
atom_types = int((lines[3])[1])
# 1 if single chain (only covalent), 2 if double chain (covalent and hydrogen)
bond_types = int((lines[4])[1])
# Should always be 1
angle_types = int((lines[5])[1])
# Relative masses. Should both be 1 for homogeneous chain
mass1 = float((lines[6])[1])
mass2 = float((lines[7])[1])
# Equilibrium length of covalent bond
chain_bond_length = float((lines[8])[1])
# Equilibrium length of hydrogen bond
H_bond_length = float((lines[9])[1])

# Writes initiation info to top of file
output = open("lammps_input", "w")
if chains == 1:
	output.write("{} {} {}\n{} {}\n{} {}\n{} {}\n\n".format
				("LAMMPS data file from restart file:","timestep = 0,",
				"procs = 1",n,"atoms",n-1,"bonds",n-2,"angles"))
elif chains == 2:
	output.write("{} {} {}\n{} {}\n{} {}\n{} {}\n\n".format
				("LAMMPS data file from restart file:","timestep = 0,",
				"procs = 1",2*n,"atoms",3*n-2,"bonds",2*n-4,"angles"))
output.write("{} {}\n{} {}\n{} {}\n\n".format(atom_types,"atom types",
			 bond_types,"bond types",angle_types,"angle types"))
output.write("{} {}{} {}\n{} {}{} {}\n{} {}{} {}\n\n".format
			(-1*box_size,box_size,"xlo","xhi",-1*box_size,box_size,"ylo","yhi",
			-1*box_size,box_size,"zlo","zhi"))

if atom_types == 1:
	output.write("{}\n\n{} {}\n".format("Masses",mass1,mass1))
elif atom_types == 2:
	output.write("{}\n\n{} {}\n{} {}\n".format("Masses",1,mass1,2,mass2))

# Sets positions of first chain atoms to 0
xpos = 0
ypos = 0
zpos = 0

# Creates empty array for atom positions
positions = np.zeros((2*n,3))

for i in range(1,n+1):
	# Prints current positions to array
	positions[i-1,0] = xpos
	positions[i-1,1] = ypos
	positions[i-1,2] = zpos
	# Generates randomly x-distance to next atom
	xrand = random.uniform(-1*(chain_bond_length),chain_bond_length)
	# Generates randomly y-distance to next atom but only accepts if 
	# bond length is =< equil. bond length
	yrand = random.uniform(-1*(chain_bond_length),chain_bond_length)
	while m.sqrt(xrand**2 + yrand**2) >= chain_bond_length:
		yrand = random.uniform(-1*(chain_bond_length),chain_bond_length)
	# Calculates z-distance to next atom to give chain correct length
	# and randomly sets positive or negative direction
	sign = 0
	while sign == 0:
		sign = random.randint(-1,1)
	zrand = (m.sqrt(chain_bond_length**2- xrand**2 - yrand**2))*sign
	# Position of next atom is position of previous atom + distances above
	xpos = xpos + xrand
	ypos = ypos + yrand
	zpos = zpos + zrand

# If 2 chains, initial positions of second chain are same in y- and z-
# directions, and x-position of 1st chain atom + H-bond length in x-direction	
if chains == 2:
	for i in range(n+1,2*n+1):
		positions[i-1,0] = positions[i-(n+1),0] + H_bond_length
		positions[i-1,1] = positions[i-(n+1),1]
		positions[i-1,2] = positions[i-(n+1),2]

# Writes atoms with labels and positions, velocities of atoms (all 0), 
# bonds with labels, and angles with labels to output file
if chains == 1:
	output.write("\nAtoms \n\n")
	atom_type = 1
	for i in range(1, n+1):
		output.write("{} {} {} {} {} {}\n".format(i,1,atom_type,positions[i-1,0],positions[i-1,1],positions[i-1,2]))
	output.write("\nVelocities \n\n")	
	for i in range(1, n+1):
		output.write("{} {} {} {}\n".format(i,0,0,0))
	# Bonds defined by two atoms bonded
	output.write("\nBonds \n\n")
	for i in range(1,n):
		output.write("{} {} {} {}\n".format(i,1,i,i+1))
	# Angles defined by three atoms anglle is between 
	output.write("\nAngles \n\n")
	for i in range(1,n-1):
		output.write("{} {} {} {} {}\n".format(i,1,i,i+1,i+2))
		
# Same for case of two chains
elif chains == 2:
	output.write("\nAtoms \n\n")
	if atom_types == 1:
		for i in range(1,2*n+1):
			atom_type = 1
			output.write("{} {} {} {} {} {}\n".format(i,1,atom_type,positions[i-1,0],positions[i-1,1],positions[i-1,2]))

	elif atom_types == 2:	
		for i in range(1,2*n+1):
			# If 2 atom types, each chain set as different type
			if i < n:
				atom_type = 1
			elif i > n:
				atom_type = 2
			output.write("{} {} {} {} {} {}\n".format(i,1,atom_type,positions[i-1,0],positions[i-1,1],positions[i-1,2]))
	
	output.write("\nVelocities \n\n")
	for i in range(1,2*n+1):
		output.write("{} {} {} {}\n".format(i,0,0,0))
	
	output.write("\nBonds \n\n")
	# Covalent bonds between atoms in first chain
	for i in range(1,n):
		output.write("{} {} {} {}\n".format(i,1,i,i+1))
	# Covalent bonds between atoms in second chain
	for i in range(n,2*n-1):
		output.write("{} {} {} {}\n".format(i,1,i+1,i+2))
	# Hydrogen bonds between corresponding atoms in opposite chains
	for i in range(2*n-1,3*n-1):
		output.write("{} {} {} {}\n".format(i,2,i-(2*n-2),i-(n-2)))
		
	output.write("\nAngles \n\n")
	# Angles only defined as angle between two covalent bonds
	for i in range(1,n-1):
		output.write("{} {} {} {} {}\n".format(i,1,i,i+1,i+2))
	for i in range(n-1,2*n-3):
		output.write("{} {} {} {} {}\n".format(i,1,i+2,i+3,i+4))

