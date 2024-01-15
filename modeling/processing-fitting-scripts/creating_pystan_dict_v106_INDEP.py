##############################
###### Import libraries ######
##############################

import numpy as np
import pandas as pd

import datetime
import argparse

import json


###################################################
###### Arguments, functions and data loading ######
###################################################

parser = argparse.ArgumentParser()

parser.add_argument("-fi", "--file_in", help="Input file (CSV)",
                    type = str, required = True, dest = "file_in")
parser.add_argument("-fo", "--folder_out", help="Output folder (JSON)",
                    type = str, required = True, dest = "folder_out")



args = parser.parse_args()
argv=vars(args)

file_in = argv['file_in']
folder_out = argv['folder_out']


### Loading data file
sub_data_full = pd.read_csv(file_in, header=0) #eg: 'data_filt_block.csv'

### Data needs to have these 8 columns!
sub_data_full = sub_data_full[['subject','block','dim1_num','dim2_len','trial','catresp','class','cat_max']]

for sub_name in np.unique(sub_data_full['subject']):

    sub_data = sub_data_full[sub_data_full['subject'] == sub_name]


    #########################################
    ###### Data processing for fitting ######
    #########################################


    ### Print stuff
    print('Starting data processing for subject',sub_name)
    e = datetime.datetime.now()
    print("The time at the beginning of processing: = %s:%s:%s" % (e.hour, e.minute, e.second))
    

    ### Modify trial list so that it only includes valid trials
    ind = np.where(np.isnan(sub_data['catresp']))[0]

    if len(ind) != 0:
        for i in ind:
            #print(i)
            trial = sub_data.iloc[i,4]
            sub = sub_data.iloc[i,0]
            bl = sub_data.iloc[i,1]
            for k in range(
                    i+1,i+max(sub_data[(sub_data['subject'] == sub) & (sub_data['block'] == bl)].iloc[:,4])-trial+1):
                #print(k)
                sub_data.iloc[k,4] -=1

    ### Drop missing trials - where participants didn't make a response          
    sub_data = sub_data.dropna(subset = ['catresp'])

    ### Sanity check - trial dropping & trial succession
    old_trial = 0
    for row in sub_data.iterrows():
        sub = row[1][0]
        bl = row[1][1]
        trial = row[1][4]
        #print(trial)
    
        if (trial != old_trial + 1) & (trial != 1):
            print('Missing trial in succession')
            print(sub, bl, trial)
        
        old_trial = trial


    ### Define lists for dictionary
    responses = np.array(list(sub_data['catresp'])).astype(int).tolist()
    #print(type(responses[0]))
    cat_max = np.array(list(sub_data['cat_max'])).astype(int).tolist()
    dim1_num = np.array(list(sub_data['dim1_num'])).astype(int).tolist() 
    dim2_len = np.array(list(sub_data['dim2_len'])).astype(int).tolist() 
    block_list = np.array(list(sub_data['block'])).astype(int).tolist()
    trial_list = np.array(list(sub_data['trial'])).astype(int).tolist()
    subject_list = np.array(list(sub_data['subject']))

    N = len(sub_data)
    NS = len(set(sub_data['subject']))
    K = max(responses)+1

    ### Create subject ID -> subject number dictionary

    subject_dict = {}
    sub_IDs = np.unique(subject_list)
    for i in range(len(sub_IDs)):
        subject_dict[sub_IDs[i]] = i+1


    #print(subject_dict)


    ### Create arrays
    nk = np.zeros((N,K)) # observation counts for each cluster

    sq_dist_num = np.zeros((N,K)) # density square distance

    sq_dist_len = np.zeros((N,K)) # length square distance

    sub_nos = np.zeros(N) # subject number

    for i in range(N):
        sub = subject_list[i]

        # count how many observations have been made of each cluster
        if trial_list[i] == 1 :
            nk[i, :] = 0 # if new trial, count is zero
        else:
            nk[i, :] = nk[i-1, :] # initialise counter from timepoint t-1
            nk[i, int(responses[i-1])-1] = nk[i, int(responses[i-1])-1] + 1 # add previous response to counter

        # calculate square distances between each stimulus value and the cluster mean values, for all existing clusters
        if cat_max[i] > 0 :
            for k in range(1,cat_max[i]+1):
                num_mean = np.mean(sub_data[(sub_data['subject']==sub)&(
                    sub_data['block'] == block_list[i])&(sub_data['trial'] < trial_list[i])&(sub_data['catresp'] == k)]['dim1_num'])
                sq_dist_num[i,k-1] = np.square(dim1_num[i] - num_mean)
    

                len_mean = np.mean(sub_data[(sub_data['subject']==sub)&(
                    sub_data['block'] == block_list[i])&(sub_data['trial'] < trial_list[i])&(sub_data['catresp'] == k)]['dim2_len'])
                sq_dist_len[i,k-1] = np.square(dim2_len[i] - len_mean)
    
        sub_nos[i] = subject_dict[sub] # store subject number at each timepoint


        ### Put data in a dictionary + convert elements in array to into & arrays to lists
    data = {'N':N,'NS':NS,'K':K, 'trial': trial_list, 'resp': responses,
            'nk':nk.astype(int).tolist(), 'kmax':cat_max, 'block_':block_list,
            'sub_nos': sub_nos.astype(int).tolist(),
        'num_sq_dist':sq_dist_num.tolist(),'len_sq_dist':sq_dist_len.tolist(),
            'dim1_num':dim1_num, 'dim2_len':dim2_len}

    #filename = 'Jsons/Individual_jsons/'+sub_name + '.json'
    filename = folder_out+sub_name + '.json'

    ### Dump dictionary into JSON file
    with open(filename, "w") as f:
        json.dump(data, f)
