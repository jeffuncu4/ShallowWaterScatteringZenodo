This code uses Dedalus to solve the rotating shallow water equations in a 2d horizontal domain. A plane wave is forced from the boundary to propagate through a vortex in cyclogeostrophic balance in the centre of the domain. The amount of wave scattering is quantified with the scattering ratio S for a variety of vortex/wave parameters. 3d fits which relate the vortex/wave characteristics to S are then calculated. The related article can be found at EarthArXiv with https://doi.org/10.31223/X5639Q and is under review in the Journal of Fluid Mechanics as of the publication of this repository.

We use anaconda to build our environment. You can download the required packages with the included environment file with 'conda env create -f env.yml". The environment can then be activated with "conda activate dedalus".

The steps to run experiments are as follows.
1. In the code files SimulationClass.py and AnalysisClass.py, change the value for "home_dir" with the folder that you intend to place this repository in, and change the value for "data_dir" with the folder which will hold experiment data. Within the folder you define as "data_dir", make the folders "experiments", "vortices" and "analytics".
2. Use run_experiment.py to run an experiment by simply modifying the following parameters: the Rossby number Ro, the Burger number Bu, the ratio of the vortex length scale to wavelength Lr, and the ratio of the vortex velocity to wave velocity Ur.
3. Because there is a cyclogeostrophic adjustment between the initial parameters and the vortex parameters, one can use measured_parameters.py to obtain the adjusted parameters for analysis, however, this step can be skipped if you conduct the analysis following our codes as outlined below. 

The steps to obtain the relation between scattering ratio S and the non-dimensional parameters are described below.
1. Run save_experiment_statistics.py to calculate S and measure the adjusted simulation parameters for a given set of experiments. The outputs will be saved to a .csv file. If you wish to calculate S for a single experiment use calculate_fluxratio.py.
2. Use RoBuK_fitting.py to do a 3D fit of the data and plot the fits.

We have various other codes which plot the simulation setup, process data with pandas, and calculate the energy spectrum but they do not need to be included for running the experiments and calculating S.
