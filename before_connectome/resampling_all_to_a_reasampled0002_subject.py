#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 26 11:24:35 2023

@author: ali
"""

import os
import nibabel as nib
from nibabel import load, save, Nifti1Image, squeeze_image
import multiprocessing
import numpy as np
import pandas as pd
import shutil
import sys
from scipy.cluster.vq import kmeans, vq
from nibabel import load, save, Nifti1Image, squeeze_image
import subprocess

ref_image = '/Users/ali/Desktop/Jun23/ADRC_mrtrix_dwifslpreproc/before_connectome/refrence_resampling_ADRC0002/HCP_DTI_reverse_phase.nii.gz'

dusom__orig_data_folders_path = '/Volumes/dusom_mousebrains/All_Staff/Projects/ADRC/Data/raw_ADRC_data/ADRC-20230511/'
list_folders_path = os.listdir(dusom__orig_data_folders_path)
list_of_subjs_long = [i for i in list_folders_path if 'ADRC' in i]
list_of_subjs = list_of_subjs_long
list_of_subjs = [str for str in list_of_subjs if len(str) < 20]

Munin_destination_path ='/Volumes/Data/Badea/Lab/ADRC-20230511/'




#### empty_folder_maker 
for subj in list_of_subjs:
    Munin_fodler = Munin_destination_path+ subj 
    if not os.path.isdir(Munin_fodler) : os.mkdir(Munin_fodler)
    Munin_fodler_visit =  Munin_destination_path+ subj  + '/visit1'
    if not os.path.isdir(Munin_fodler_visit) : os.mkdir(Munin_fodler_visit)
    

    


#### DTI 
for subj in list_of_subjs:
    subj_DTI_orig_path = dusom__orig_data_folders_path+ subj + "/visit1/HCP_DTI.nii.gz"
        #a= int(str(subprocess.check_output("mrinfo -size "+ subj_T11_orig_path, shell=True))[2:5])
        #if a==256: 
        #DTI_temp = Munin_destination_path+subj+'/visit1/HCP_DTI_temp.nii.gz'
        #os.system( 'mrgrid '+ subj_DTI_orig_path + ' regrid -voxel 1.5 '+ DTI_temp +' -force  >/dev/null 2>&1')
    DTI_munin_destination =  Munin_destination_path+subj+'/visit1/HCP_DTI.nii.gz'
    command =  'mrgrid '+subj_DTI_orig_path+' regrid -template ' + ref_image  + ' -scale 1,1,1 '    + DTI_munin_destination + ' -force'

    os.system( command + '>/dev/null 2>&1')
    
    
#### reverse_DTI 
for subj in list_of_subjs:
    subj_rev_DTI_orig_path = dusom__orig_data_folders_path+ subj + "/visit1/HCP_DTI_reverse_phase.nii.gz"
        #a= int(str(subprocess.check_output("mrinfo -size "+ subj_T11_orig_path, shell=True))[2:5])
        #if a==256: 
        #DTI_temp = Munin_destination_path+subj+'/visit1/HCP_DTI_temp.nii.gz'
        #os.system( 'mrgrid '+ subj_DTI_orig_path + ' regrid -voxel 1.5 '+ DTI_temp +' -force  >/dev/null 2>&1')
    DTI_rev_munin_destination =  Munin_destination_path+subj+'/visit1/HCP_DTI_reverse_phase.nii.gz'
    command =  'mrgrid '+subj_rev_DTI_orig_path+' regrid -template ' + ref_image  + ' -scale 1,1,1 '    + DTI_rev_munin_destination + ' -force'

    os.system( command + '>/dev/null 2>&1')    



    
#### T1s 
for subj in list_of_subjs:
    subj_T1_orig_path = dusom__orig_data_folders_path+ subj + "/visit1/T1.nii.gz"
        #a= int(str(subprocess.check_output("mrinfo -size "+ subj_T11_orig_path, shell=True))[2:5])
        #if a==256: 
        #DTI_temp = Munin_destination_path+subj+'/visit1/HCP_DTI_temp.nii.gz'
        #os.system( 'mrgrid '+ subj_DTI_orig_path + ' regrid -voxel 1.5 '+ DTI_temp +' -force  >/dev/null 2>&1')
    T1_munin_destination =  Munin_destination_path+subj+'/visit1/T1.nii.gz'
    command =  'mrgrid '+subj_T1_orig_path+' regrid -template ' + ref_image  + ' -scale 1,1,1 '    + T1_munin_destination + ' -force'
    print(subj)
    os.system( command + '>/dev/null 2>&1')    

#### all other files to copy




for subj in list_of_subjs:
    other_files_orig_path = dusom__orig_data_folders_path + subj + '/visit1'
    other_files_orig_names = os.listdir(other_files_orig_path)
    other_files_orig_names = [ i for i in other_files_orig_names if not "T1.nii.gz" in i]
    other_files_orig_names = [ i for i in other_files_orig_names if not "HCP_DTI.nii.gz" in i]
    other_files_orig_names = [ i for i in other_files_orig_names if not "HCP_DTI_reverse_phase.nii.gz" in i]
    for file_names in other_files_orig_names:        
        orig_file =   other_files_orig_path +  '/' +   file_names
        destination_folder = Munin_destination_path+subj+'/visit1/'
        command =  'cp '+orig_file+' ' + destination_folder
        os.system( command )    
    print(subj)


