'''
Author: Dan-Mircea Mirea
Date: 07/18/22 (mod 11/17/22)
Scriptescription: Python script that performs PySTAN fitting for Gaussian-likelihood-CRP-based latent cause models with variable (1) pystan sampling function specification, (2) seed and (3) model type; corrected subject number bug;

CmdStanPy-based!!!

'''

##############################
###### Import libraries ######
##############################

#import stan

import numpy as np

import datetime
import argparse
import os

#import packages
from cmdstanpy import cmdstan_path, CmdStanModel, rebuild_cmdstan

###################################################
###### Arguments, functions and data loading ######
###################################################

parser = argparse.ArgumentParser()

parser.add_argument("-ni", "--n_iter", help="Number of iterations",
                    type = int, required = True, dest = 'n_iter')
parser.add_argument("-nc", "--n_chains", help="Number of chains",
                    type = int, required = True, dest = 'n_chains')
parser.add_argument("-s","--seed_val", help="Desired seed",
                    type = int, required = True, dest = "seed_val")
parser.add_argument("-m", "--model_type", help="Type of model to be fit",
                    type = str, required = True, dest = 'model_type')
parser.add_argument("-si", "--start_index", help="First JSON file to model",
                    type = str, required = True, dest = "start_index")
parser.add_argument("-ei", "--end_index", help="Last JSON file to model",
                    type = str, required = True, dest = "end_index")
parser.add_argument("-f", "--folder_path", help="Folder path for indiv. jsons",
                    type = str, required = True, dest = "folder_path")


args = parser.parse_args()
argv=vars(args)

n_iter = argv['n_iter']
n_chains = argv['n_chains']
seed_val = argv["seed_val"]
model_type = argv['model_type']
start_index = argv['start_index']
end_index = argv['end_index']
folder_path = argv['folder_path']

n_warmup = int(n_iter/5) #1/5 warmup iterations


files = os.listdir('Jsons_extra/'+folder_path)
file_range = files[int(start_index):int(end_index)]

folder_loc = '/jukebox/niv/dannie/QIP_fitting/'
print(folder_loc)

for file_ in file_range:
    
    data_file = folder_loc+'Jsons_extra/'+folder_path+file_

    model_file = folder_loc+model_type+'.stan' #add file extension

    ###rebuild_cmdstan()


    ### Make folder for fitting results if it doesn't exist
    if not os.path.exists('Fits'):
        os.makedirs('Fits')

    ### Define folder name 
    folder_name = 'Fits/' + folder_path+ model_type + '_n_iter=' + str(
            n_iter) + '_n_chains=' + str(n_chains) + '_seed=' + str(seed_val)

    if not os.path.exists(folder_name):
        os.makedirs(folder_name)

    model_name = file_.split('.')[0]

    # print(data)

    #####################
    ###### Fitting ######
    #####################

    ### Print stuff

    print('Starting fitting from file',data_file)
    e = datetime.datetime.now()
    print("The time at the beginning of fitting: = %s:%s:%s" % (e.hour, e.minute, e.second))



    # Instantiate the Stan model, assemble the data ---------------------
    # specify Stan program file
  

    print(cmdstan_path())

    # instantiate the model; compiles the Stan program as needed.
    model = CmdStanModel(stan_file=model_file, model_name = model_name)

    # inspect model object
    print(model)

    # Run the HMC-NUTS sampler ------------------------------------------
    # specify data file

    
    # fit the model
    fit = model.sample(data=data_file, chains = n_chains,
                    iter_sampling = n_iter, iter_warmup = n_warmup,
                    seed = seed_val, output_dir=folder_name,
                       save_latent_dynamics = True)

    # printing the object reports sampler commands, output files
    print(fit)

    # Access the sample -------------------------------------------------
    fit.draws().shape
    fit.draws(concat_chains=True).shape

    # Summarize the results ---------------------------------------------
    fit.summary()
    fit.diagnose()

    # Save the Stan CSV files
    #fit.save_csvfiles(dir='Fits/')
