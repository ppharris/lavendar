#######################################
# File to specify setup of experiment #
#######################################
# core py modules:
import os
# local modules:
import observations
import plot


# set directory for JULES output
output_directory = os.getcwd()+'/output'
# set directory containing input JULES NML files
nml_directory = os.getcwd()+'/example_nml'
# set directory for temporary JULES nml files
work_directory = "/scratch/ppha/lavendar_test/work"
# set model executable
model_exe = '/home/users/ewanp82/models/jules4.9/build/bin/jules.exe'
# set function to extract JULES modelled observations for prior JULES
jules_hxb = observations.extract_jules_hxb
# set function to extract prior ensemble of modelled observations
jules_hxb_ens = observations.extract_jules_hxb_ens
# set function to extract observations to be assimilated
obs_fn = observations.extract_twin_data
# set JULES parameters to optimised during data assimilation
opt_params = {'pft_params': {
                  'jules_pftparm': {
                      'neff_io': [7, 6.24155040e-04, (5e-05, 0.0015)],
                      'alpha_io': [7, 6.73249126e-02, (0, 1.0)],
                      'fd_io': [7, 8.66181324e-03, (0.0001, 0.1)]}},
              'crop_params': {
                  'jules_cropparm': {
                      'gamma_io': [2, 2.07047321e+01, (0.0, 40.0)],
                      'delta_io': [2, -2.97701647e-01, (-2.0, 0.0)],
                      'mu_io': [2, 2.37351160e-02, (0.0, 1.0)],
                      'nu_io': [2, 4.16006288e+00, (0.0, 20.0)]}}}
# set error on prior parameter estimates
prior_err = 0.25
# set size of ensemble to be used in data assimilation experiments
ensemble_size = 50
# set number of processors to use in parallel runs of JULES ensemble
num_processes = 100
# plotting save function
save_plots = plot.save_plots
