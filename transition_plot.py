"""
transition_plot.py - Cameron Johnston
This script:
- Takes an input file containing average fractions of remaining bonds 
  between two DNA chains at different values of Umin and plots this with
  fixed chain length N and bond stiffness K_B.
- Calculates the Umin value at the transition where half of bonds remain
  and the error on this value, and writes them to a file.
- Plots the Umin transition value against K_B for each chain length.
"""

import numpy as np
import sys
import matplotlib.pyplot as pl

def transition_calc(X, Y, Yerr):
	# Method to caculate transition Umin and error
	# Finds where (fraction of remaining bonds) - 0.5 becomes positive
	diff = Y - 0.5
	for i in range(diff.size):
		if diff[i] < 0:
			# 1000 is arbitrarily large value so smallest value marks transition
			diff[i] = 1000
	# Finds smallest value and takes its index
	index = np.argmin(diff)
	# Finds values which transition falls between
	Xmin, Xmax = X[index-1], X[index]
	Ymin, Ymax = Y[index-1], Y[index]
	Yminerr, Ymaxerr = Yerr[index-1], Yerr[index]
	# Uses gradient of staight line between points to find where this
	# line crosses Y = 0.5
	Umin_transition = Xmin + ((0.5-Ymin)*(Xmax-Xmin)/(Ymax-Ymin))
	# Calculates error by taking partial derivatives wrt to each 
	# value in calculation which has an error
	dXdYmin = ((Xmax-Xmin)/(Ymax-Ymin))*(((0.5-Ymin)/(Ymax-Ymin))-1)
	dXdYmax = ((Xmax-Xmin)/(Ymax-Ymin))*((Ymin-0.5)/(Ymax-Ymin))
	Uminerr = np.sqrt((Yminerr**2)*(dXdYmin**2) + (Ymaxerr**2)*(dXdYmax**2))
	
	return Umin_transition, Uminerr
	

def plot(num):
	# This method plots average remaining bonds (theta) against Umin
	num = int(num)
	# Method works for any number of input files
	infile = open(sys.argv[num], 'r')
	data = np.loadtxt(infile)
	# Sorts data by Umin value
	data = data[data[:,0].argsort()]
	# Creates arrays of relevant data
	Umin_full, ave_bonds_full, sd_bonds_full = data[:,0], data[:,1], data[:,2]
	# Empty arrays created for averages
	Umin = np.array([])
	ave_bonds = np.array([])
	errors = np.array([])
	# Counts how many values there are for each Umin (assumes all Umins
	# have same number of repetitions)
	repetitions = 1
	for i in range(1,Umin_full.size):
		if Umin_full[i] == Umin_full[0]:
			repetitions = repetitions + 1
	
	# Averages theta and error for each value of Umin	
	for i in range(Umin_full.size):
		if (i%repetitions) == 0:
			Umin = np.append(Umin, Umin_full[i])
			bonds = np.array([])
			sd = np.array([])
			for j in range(repetitions):
				bonds = np.append(bonds, ave_bonds_full[i+j])
				sd = np.append(sd, sd_bonds_full[i+j])
			mean = np.average(bonds)
			sd_sq = sd**2
			std_err_mean = (np.sqrt(np.sum(sd_sq)))/(repetitions*np.sqrt(repetitions))
			ave_bonds = np.append(ave_bonds, mean)
			errors = np.append(errors, std_err_mean)
	# Standard error on mean found as (average std dev)/sqrt(number of repetitions)
	# Finds values of N and K_B from input file name
	# so input file names MUST have format bonds-<N>-<K_B>.dat
	name = infile.name
	name = name.strip("bonds-")
	name = name.strip(".dat")
	values = name.split("-")
	length, K = int(values[0]), int(values[1])
	# Plots data with error bars and assigns label
	pl.errorbar(Umin, ave_bonds, xerr = None, yerr = errors, capsize = 2, label = ("N = "+str(length)+", "+r"$K_B=$"+str(K)))
	
	return length, K, Umin, ave_bonds, errors
	

def main():
	
	# Creates output file for calculated transition values
	outfile = open('transition_points.dat', 'w')
	outfile.write("{}	{}	{}	{}\n".format("N", "K", "Umin at transition", "Error on Umin"))
	
	# Creates empty lists to later plot transition data
	K100, K200, K500 = [], [], []
	U100, U200, U500 = [], [], []
	Uerr100, Uerr200, Uerr500 = [], [], []

	# Plots Y = 0.5 to show transition points
	pl.axhline(y=0.5, color = 'r', linestyle = '--', label = r"$\theta = 0.5$")
	# Performs plot() and transition_calc() methods for every input file
	for i in range(1,(len(sys.argv))):
		length, K, Umin, ave_bond, sd_err_mean = plot(i)
		Umin_transition, Uminerr = transition_calc(Umin, ave_bond, sd_err_mean)
		# Writes transition data to output file
		outfile.write("{}	{}	{}	{}\n".format(length, K, float(Umin_transition), float(Uminerr)))
		# Appends transition data to relevant lists
		if length == 100:
			K100.append(K)
			U100.append(Umin_transition)
			Uerr100.append(Uminerr)
		if length == 200:
			K200.append(K)
			U200.append(Umin_transition)
			Uerr200.append(Uminerr)
		if length == 500:
			K500.append(K)
			U500.append(Umin_transition)
			Uerr500.append(Uminerr)
	pl.xlabel(r"$U_{min}$ ($\epsilon$)")
	pl.ylabel(r"$\theta$")
	pl.legend()
	# Shows plot with linear axes
	pl.show()

	# Creates same plot again with logarithmic x-axis
	pl.axhline(y=0.5, color = 'r', linestyle = '--', label = r"$\theta = 0.5$")
	for i in range(1,(len(sys.argv))):
		plot(i)
	pl.xlabel(r"$U_{min}$ ($\epsilon$)")
	pl.ylabel(r"$\theta$")
	pl.xscale("log")
	pl.legend()
	pl.show()
	
	#Plots transition data with line for each value of N
	pl.errorbar(K100, U100, yerr = Uerr100, capsize = 2, label = r"$N=100$")
	pl.errorbar(K200, U200, yerr = Uerr200, capsize = 2, label = r"$N=200$")
	pl.errorbar(K500, U500, yerr = Uerr500, capsize = 2, label = r"$N=500$")
	pl.xlabel(r"$K_B$ ($\epsilon$)")
	pl.ylabel(r"$U_{min}^{transition}$ ($\epsilon$)")
	pl.legend()
	pl.show()

main()
