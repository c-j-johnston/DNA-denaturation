# DNA-denaturation
Python scripts written and modified LAMMPS script for DNA Denaturation project.

The repository contains the following:

1. input.txt - a text file with input parameters for input_file_create.py.
2. input_file_create.py - a Python script which randomly generates an input for the LAMMPS script using the parameters in input.txt. Initial positions of monomers in a single or double chain are randomly chosen.
3. lammps_input - an example input file generated by input_file_create.py. File must have this name to be read by denaturation.lam.
4. denaturation.lam - LAMMPS script used to run simulations of DNA denaturation. Bond stiffness and Umin should be changed to run different experiments. Parameters to change are shown in comments in code. LAMMPS software must be downloaded to run this script.
5. separation_count.py - Python script taking dump file from denaturation.lam which counts the number of remainng bonds at each simulation timestep and averages these values past 1 million timesteps. Also finds standard deviation of this value and plots fraction of remaining bonds over time.
6. transition_plot.py - takes data file containing average remaining bond fraction (theta) over many Umin values for a fixed chain length and plots theta vs Umin. Finds value of Umin for transition from melted to unmelted phase and errors. Plots bond stiffness against transition Umin for different chain lengths.
