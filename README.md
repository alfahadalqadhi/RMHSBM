# RMHSBM
Code that accompanies a paper on RMHSBM, the code includes the simulation and an application on the BNU1 dataset.
Data can be downloaded from: https://www.cis.jhu.edu/~parky/Microsoft/JHU-MSR/ZMx2/BNU1/DS16784.

# Full Simulation
Run runsim.R and runsim_nocorr.R to create the simulation R files. Run Make_Models.R to create the models data files, you can edit the nDs vector to change the number of mismatched blocks. Now you can the files you created with runsim.R, they will produce the data from the simulations. Running MakePlotData.R makes a data file that can be used by CreatePlots.R to make plots for different numbers of mismatched blocks. RMHSBM_functions.R includes helper functions that you need.

# Reproducing the simulation plots from the paper
You can reproduce plots 4-9 from the paper by running Make_Models.R, runsim18.R and runsim18_nocorr.R to creates the data files. Then run MakePlotData.R and CreatePlots.R; change line 68 in MakePlotData.R to for(simnum in c(18)){ .

# For the BNU1 data analysis
Run BNU1Analysis.R then BNU1plots.R  after downloading the data.
