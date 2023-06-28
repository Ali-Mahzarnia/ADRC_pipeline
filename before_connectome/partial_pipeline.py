#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:21:13 2023

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
#nthreads = multiprocessing.cpu_count()
nthreads = 8;
#n_threds = str(multiprocessing.cpu_count())

subj = sys.argv[1] #reads subj number with s... from input of python file 
#subj = "ADRC0001"

#subj_ref = subj



#subj = "H21593" #reads subj number with s... from input of python file 

index_gz = ".gz" 
#index_gz = "" 


root= '###/ADRC_mrtrix_dwifsl/'



orig_subj_path = '###/ADRC-20230511/'+ subj + "/visit1/"



for file in os.listdir(orig_subj_path):     #rename files so there's no space
    os.rename(orig_subj_path+ file ,orig_subj_path +file.replace(" ", "_"))



DTI_forward_nii_gz =  orig_subj_path  + 'HCP_DTI.nii.gz'
if not os.path.isfile(DTI_forward_nii_gz) : print('where is original 4d?')

DTI_reverse_phase_nii_gz = orig_subj_path + "HCP_DTI_reverse_phase.nii.gz"
if not os.path.isfile(DTI_reverse_phase_nii_gz) : print('where is original reverse DTI?')

path_perm = root + 'perm_files/'
if not os.path.isdir(path_perm) : os.mkdir(path_perm)
#bval_path = path_perm  + subj + '_bvals_RAS.txt'
#old_bval = np.loadtxt(bval_path_orig)
#new_bval = np.round(old_bval)
#new_bval.shape
#np.savetxt(bval_path, new_bval, newline=" ", fmt='%f') # saving the RAS bvec

#bval_path = '###/Downloads/N59066_bval.txt'



bvec = [] # bvec extraction
with open(orig_subj_path + "HCP_DTI.bxh") as file:
    while line := file.readline():
        if "<value>" in line:
            temp_loop = line
            temp_loop_split = temp_loop.split()
            temp_loop_2  = [ i.replace("<value>" , "") for i in temp_loop_split ]
            temp_loop_2  = [ i.replace("</value>" , "") for i in temp_loop_2 ]
            temp_loop_2 = [float(i) for i in temp_loop_2]
            bvec.append(temp_loop_2)
            #print( temp_loop_2)
            
            
norms  = np.linalg.norm(bvec, axis= 1)
norms [norms ==0 ] = 1 
bvec = np.array(bvec)
bvec = bvec /   norms.reshape(len(norms),1)   
bvec = bvec.transpose()    

 
bvec_path_orig = path_perm + subj +'_bvec.txt' 
np.savetxt(bvec_path_orig ,bvec,fmt='%.2f')
#os.system("awk '{print NF; exit}' "+bvec_path_orig) #print the lenght of bvec : must match the previous number


bvals = [] # bval extraction
with open(orig_subj_path + "HCP_DTI.bxh") as file:
    while line := file.readline():
        if "bvalues" in line:
            temp_loop = line
            temp_loop_split = temp_loop.split()
            temp_loop_2  = [ i.replace("<bvalues>" , "") for i in temp_loop_split ]
            temp_loop_2  = [ i.replace('</bvalues>' , "") for i in temp_loop_2 ]

temp_loop_2 = [float(i) for i in temp_loop_2]            
bvals =  temp_loop_2
bvals  =  bvals *(norms**2)

codebook, _ = kmeans(bvals, 3) 
ccodebook2 = codebook.round()
cluster_indices, _ = vq(bvals, codebook)
rnd_bvals= ccodebook2[cluster_indices]
bvals = rnd_bvals
bvals = bvals.reshape(1,len(bvals))
new_bval = bvals
bval_path_orig = path_perm + subj +'_bvals.txt'
np.savetxt(bval_path_orig ,bvals,fmt='%.2f')
#os.system("awk '{print NF; exit}' "+bval_path_orig) #print the lenght of bval : must match the previous number



#os.system( 'bash co_reg_4d_stack_tmpnew.bash $diff_path $subj_name 0 $outpath 0') 


if not os.path.isdir(root +  'temp/' ) : os.mkdir(root +  'temp/' )
subj_path = root +  'temp/' + subj + '/'
if not os.path.isdir(subj_path) : os.mkdir(subj_path)
#timeseries_input_fldr = root +  'temp/' + subj + '/input_series/'
#if not os.path.isdir(timeseries_input_fldr) : os.mkdir(timeseries_input_fldr)
#timeseries_coreg_fldr = root +  'temp/' + subj + '/coreg_series/'
#if not os.path.isdir(timeseries_coreg_fldr) : os.mkdir(timeseries_coreg_fldr)



#os.system('/Applications/ANTS/ImageMath 4 ' + timeseries_input_fldr+subj+'_.nii.gz' + ' TimeSeriesDisassemble '+ DTI_forward_nii_gz) #dissassmple time series
#ref_index   = [index for index, value in enumerate(bvals) if value == 0] # find index of b0  : wherever bval is 0
#ref_index = np.where(bvals == 0)[1]

#b_0_ref_nii_gz =  timeseries_input_fldr + subj + '_' +  str(ref_index[0]  + 1000)+ '.nii.gz'




#for i in range(bvals.shape[1]): # register all volumes to b0 volume
#    b_i_ref_nii_gz =  timeseries_input_fldr + subj + '_' +  str(i + 1000)+ '.nii.gz'
#    out_prefix = timeseries_coreg_fldr + subj + '_' +  str(i + 1000) + "_"
#    os.system("/Applications/ANTS//antsRegistration -v 1 -d 3 -m Mattes["+b_0_ref_nii_gz+" ,"+b_i_ref_nii_gz+",1,32,None] -r ["+b_0_ref_nii_gz+" ,"+b_i_ref_nii_gz+",1] -t affine[0.1] -c [300x300x0x0,1e-8,20] -s 4x2x1x0.5vox -f 6x4x2x1 -u 1 -z 1 -l 1 -o "+out_prefix + " >/dev/null 2>&1")  
#    os.system("/Applications/ANTS/antsApplyTransforms -d 3 -e 0 --float  -u float -i "+b_i_ref_nii_gz+" -r "+b_0_ref_nii_gz+" -o "+out_prefix+"coregs.nii.gz -t "+out_prefix+"0GenericAffine.mat" + " >/dev/null 2>&1") 
    #print(i)




#os.system( "cp " + b_0_ref_nii_gz + " "+ timeseries_coreg_fldr + subj + '_' +  str(0 + 1000)+ '_coregs.nii.gz')
#shutil.rmtree(timeseries_input_fldr )

#DTI_coreg_nii_gz =  path_perm + subj + 'coreg_DTI.nii.gz'
#os.system('/Applications/ANTS/ImageMath 4 ' + DTI_coreg_nii_gz + ' TimeSeriesAssemble 1 0 '+ timeseries_coreg_fldr+subj+'*.nii.gz') #assmple time series
#shutil.rmtree(timeseries_coreg_fldr )


#if not os.path.isfile(orig_subj_path + "HCP_DTI.bxh") : print('where is original bvec, bval text file?')







T1_orig =  orig_subj_path+'T1.nii.gz'
if not os.path.isfile(T1_orig) : print('where is original DWI?')






#b0_orig =  orig_subj_path+subj+'_subjspace_b0.nii.gz'
#if not os.path.isfile(b0_orig) : print('where is original b0?')











T1 =T1_orig

#os.system('/Applications/Convert3DGUI.app/Contents/bin/c3d ' + T1_orig +" -orient RAS -o "+T1) # RAS rotation T1

#if not os.path.isdir(root +  'temp/' ) : os.mkdir(root +  'temp/' )
#subj_path = root +  'temp/' + subj + '/'
#if not os.path.isdir(subj_path) : os.mkdir(subj_path)

#diff=nib.load(nii_gz_path_orig) # read the data of this volumetrtic file as nib object
#diff_data=diff.get_fdata() #read data as array 

#for i in range(int(diff_data.shape[3])):
#    diff_data_3d=diff_data[:,:,:,i]
#    squuezed=squeeze_image(nib.Nifti1Image(diff_data_3d,diff.affine)) #squeeze the last dimension
#    vol_ith_path = subj_path + subj + '_'+str(i)+'.nii.gz'
#    nib.save(squuezed, vol_ith_path) 
#    vol_ith_RAS_path = subj_path + subj + '_'+str(i)+'_RAS.nii.gz'
#    os.system("/Applications/Convert3DGUI.app/Contents/bin/c3d "+vol_ith_path+" -orient RAS -o "+vol_ith_RAS_path)

#volS_RAS_path = subj_path + subj + '_' 
#nii_gz_path =  path_perm  + subj +'_RAS_coreg.nii.gz'
#os.system(f'###/Downloads/ANTsR/install/bin/ImageMath 4 {nii_gz_path} TimeSeriesAssemble 1 0 {volS_RAS_path}*RAS.nii.gz') # concatenate volumes of fmri
#nii_gz_path = nii_gz_path_orig

#bvec_path = path_perm+subj+'_bvecs_RAS.txt' 
old_bvec = np.loadtxt(bvec_path_orig)
new_bvec = old_bvec
#new_bvec = old_bvec [:, [2,1,0] ] # swap x and y
new_bvec[1, 1:] = -new_bvec[1, 1:] # flip y sign
#new_bvec[1:,1] = -new_bvec[1:,1] # flip x sign
#new_bvec[1:,2] = -new_bvec[1:,2] # flip z sign
#new_bvec=new_bvec.transpose()
np.savetxt(bvec_path_orig, new_bvec, fmt='%f') # saving the RAS bvec
#bvec_path = bvec_path_orig
#bvec_path  = '###/Downloads/N59066_bvecs.txt'


###
bvec_path = bvec_path_orig
bval_path = bval_path_orig
####



#changng to mif format
T1_mif = subj_path+subj+'_T1.mif'+index_gz
os.system('mrconvert ' +T1+ ' '+T1_mif + ' -force' )

out_mif = subj_path + subj+'_subjspace_dwi.mif'+index_gz
os.system('mrconvert '+DTI_forward_nii_gz+ ' ' +out_mif+' -fslgrad '+bvec_path_orig+ ' '+ bval_path_orig+' -bvalue_scaling false -force') #turn off the scaling otherwise bvals becomes 0 4000 1000 instead of 2000 
    
#os.system('mrinfo '+out_mif+ ' -export_grad_fsl ' + '###/Desktop/Feb23/mrtrix_pipeline/temp/N59141/test.bvecs ###/Desktop/Feb23/mrtrix_pipeline/temp/N59141/test.bvals -force'  ) #print the 4th dimension



#os.system('mrinfo '+out_mif) #info





#os.system('mrinfo -size '+out_mif+" | awk '{print $4}' ") #print the 4th dimension 
#os.system("awk '{print NF; exit}' "+bvec_path) #print the lenght of bvec : must match the previous number
#os.system("awk '{print NF; exit}' "+bval_path) #print the lenght of bval : must match the previous number
#os.system("mrview "+ out_mif) #viewer 
 
#preprocessing
#denoise:

output_denoise =   subj_path+subj+'_den.mif'
os.system('dwidenoise '+out_mif + ' ' + output_denoise+' -force') #denoising
#compute residual to check if any resdiual is loaded on anatomy
output_residual = subj_path+subj+'residual.mif'
os.system('mrcalc '+out_mif + ' ' + output_denoise+ ' -subtract '+ output_residual+ ' -force') #compute residual
#os.system('mrview '+ output_residual) #inspect residual







#distortion correction
nii_gz_path_PA = orig_subj_path +'HCP_DTI_reverse_phase.nii.gz'


# read bvec bval of reverse phase:


bvec_rvrs = [] # bvec extraction
with open(orig_subj_path + "HCP_DTI_reverse_phase.bxh") as file:
    while line := file.readline():
        if "<value>" in line:
            temp_loop = line
            temp_loop_split = temp_loop.split()
            temp_loop_2  = [ i.replace("<value>" , "") for i in temp_loop_split ]
            temp_loop_2  = [ i.replace("</value>" , "") for i in temp_loop_2 ]
            temp_loop_2 = [float(i) for i in temp_loop_2]
            bvec_rvrs.append(temp_loop_2)
            #print( temp_loop_2)
            
            
            

norms  = np.linalg.norm(bvec_rvrs, axis= 1)
norms [norms ==0 ] = 1 
bvec_rvrs  = np.array(bvec_rvrs )
bvec_rvrs = bvec_rvrs /   norms.reshape(len(norms),1)   
bvec_rvrs  = bvec_rvrs .transpose()   

bvec_rvrs[1, 1:] = -bvec_rvrs[1, 1:]



bvec_rvrs_path = path_perm + subj +'_bvec_rvrs.txt'


 
np.savetxt(bvec_rvrs_path ,bvec_rvrs,fmt='%.2f')


bvals_rvrs = [] # bval extraction
with open(orig_subj_path + "HCP_DTI_reverse_phase.bxh") as file:
    while line := file.readline():
        if "bvalues" in line:
            temp_loop = line
            temp_loop_split = temp_loop.split()
            temp_loop_2  = [ i.replace("<bvalues>" , "") for i in temp_loop_split ]
            temp_loop_2  = [ i.replace('</bvalues>' , "") for i in temp_loop_2 ]

temp_loop_2 = [float(i) for i in temp_loop_2]            
bvals_rvrs =  temp_loop_2
bvals_rvrs  =  bvals_rvrs *(norms**2)

codebook, _ = kmeans(bvals_rvrs, 3) 
ccodebook2 = codebook.round()
cluster_indices, _ = vq(bvals_rvrs, codebook)
rnd_bvals= ccodebook2[cluster_indices]
bvals = rnd_bvals
bvals = bvals.reshape(1,len(bvals_rvrs))

bval_path_rvrs = path_perm + subj +'_bvals_rvrs.txt'
np.savetxt(bval_path_rvrs ,bvals,fmt='%.2f')


bvec_path_PA = bvec_rvrs_path
bval_path_PA = bval_path_rvrs







############


PA_mif = subj_path+subj+'_PA.mif'
os.system('mrconvert '+nii_gz_path_PA+ ' ' +PA_mif + ' -force') #PA to mif    


mean_b0_PA_mif = subj_path+subj+'mean_b0_PA.mif'
    #take mean of PA and AP to unwarp

os.system('mrconvert '+PA_mif+ ' -fslgrad '+bvec_path_PA+ ' '+ bval_path_PA + ' - | mrmath - mean '+ mean_b0_PA_mif+' -axis 3 -force')

#Extracting b0 images from the AP dataset, and concatenating the b0 images across both AP and PA images:
mean_b0_AP_mif = subj_path+subj+'mean_b0_AP.mif'
os.system('dwiextract '+output_denoise+ ' - -bzero | mrmath - mean '+ mean_b0_AP_mif+' -axis 3 -force')
b0_pair_mif = subj_path+subj+'b0_pair.mif'
os.system('mrcat '+mean_b0_AP_mif+ ' ' +mean_b0_PA_mif+' -axis 3 '+ b0_pair_mif+' -force')





den_preproc_mif = subj_path+subj+'_den_preproc.mif' 



############################# implementing dwifslpreproc
#dwifslpreproc: Generated scratch directory: ###/ADRC_mrtrix/temp/ADRC0001/dwifslpreproc-tmp-CJQZIR/
scratch_path = subj_path + 'dwifsl_scratch/'
if not os.path.isdir(scratch_path) : os.mkdir(scratch_path)


#Command:  mrconvert ###/ADRC_mrtrix/temp/ADRC0001/ADRC0001_den.mif ###/ADRC_mrtrix/temp/ADRC0001/dwifslpreproc-tmp-CJQZIR/dwi.mif -fslgrad ###/ADRC_mrtrix/perm_files/ADRC0001_bvec.txt ###/ADRC_mrtrix/perm_files/ADRC0001_bvals.txt -json_export ###/ADRC_mrtrix/temp/ADRC0001/dwifslpreproc-tmp-CJQZIR/dwi.json
#          mrconvert ###/Desktop/Jun23/adrc_dwifsl/temp/ADRC0001/ADRC0001_den.mif ###/Desktop/Jun23/adrc_dwifsl/temp/ADRC0001/dwifsl_scratch/dwi.mif -fslgrad ###/Desktop/Jun23/adrc_dwifsl/perm_files/ADRC0001_bvec.txt ###/Desktop/Jun23/adrc_dwifsl/perm_files/ADRC0001_bvals.txt -json_export ###/Desktop/Jun23/adrc_dwifsl/temp/ADRC0001/dwifsl_scratch/dwi.json -force
dwi_mif = scratch_path + 'dwi.mif'
dwi_json = scratch_path + 'dwi.json'
command = 'mrconvert ' + output_denoise+ " " + dwi_mif+' -fslgrad ' + bvec_path + ' '+ bval_path + ' -json_export '+ dwi_json  + ' -force'
print(command)
os.system(command)


#Command:  mrconvert ###/ADRC_mrtrix/temp/ADRC0001/ADRC0001b0_pair.mif ###/ADRC_mrtrix/temp/ADRC0001/dwifslpreproc-tmp-CJQZIR/se_epi.mif
se_epi_mif  =  scratch_path + 'se_epi.mif'
command = 'mrconvert ' + b0_pair_mif+ " " + se_epi_mif+' -force'
print(command)
os.system(command)

#dwifslpreproc: Changing to scratch directory (###/ADRC_mrtrix/temp/ADRC0001/dwifslpreproc-tmp-CJQZIR/)
#dwifslpreproc: Total readout time not provided at command-line; assuming sane default of 0.1

#Command:  mrinfo dwi.mif -export_grad_mrtrix grad.b
grad_b  = scratch_path+'grad.b'
command = 'mrinfo ' + dwi_mif+ ' -export_grad_mrtrix '+grad_b  +' -force'
print(command)
os.system( command )

#dwifslpreproc: 1 spatial axis of DWIs has non-even size; this will be automatically padded for compatibility with topup, and the extra slice erased afterwards

#Command:  mrconvert se_epi.mif -coord 2 76 - | mrcat se_epi.mif - se_epi_pad2.mif -axis 2
se_epi_pad2_mif = scratch_path + 'se_epi_pad2.mif'
end_vol = subprocess.check_output('mrinfo -size '+out_mif+" | awk '{print $3}' ", shell=True)
end_vol = end_vol.decode()
end_vol = str(int(end_vol[0:2])-1)
print(end_vol)
command = 'mrconvert ' + se_epi_mif + ' -coord 2 '+end_vol+' - | mrcat '+ se_epi_mif + ' - ' + se_epi_pad2_mif + ' -axis 2 -force'
print(command)
os.system(command)

#Command:  mrconvert dwi.mif -coord 2 76 -clear dw_scheme - | mrcat dwi.mif - dwi_pad2.mif -axis 2
dwi_pad2_mif = scratch_path + 'dwi_pad2.mif'
command = 'mrconvert '+ dwi_mif + ' -coord 2 '+end_vol+'  -clear dw_scheme - | mrcat ' + dwi_mif + ' - ' + dwi_pad2_mif+ ' -axis 2 -force'
print(command)
os.system(command)

#Command:  mrconvert se_epi_pad2.mif topup_in.nii -import_pe_table se_epi_manual_pe_scheme.txt -strides -1,+2,+3,+4 -export_pe_table topup_datain.txt
topup_in_nii = scratch_path + 'topup_in.nii'
topup_datain_txt =  scratch_path + 'topup_datain.txt'
se_epi_manual_pe_scheme_txt = root + 'se_epi_manual_pe_scheme.txt'
command = 'mrconvert '+ se_epi_pad2_mif +' ' + topup_in_nii + ' -import_pe_table ' + se_epi_manual_pe_scheme_txt + ' -strides -1,+2,+3,+4 -export_pe_table '+ topup_datain_txt + ' -force'
print(command)
os.system(command)

#Command:  topup --imain=topup_in.nii --datain=topup_datain.txt --out=field --fout=field_map.nii.gz --config=/usr/local/packages/fsl-6.0.3/etc/flirtsch/b02b0.cnf --verbose
field_map_nii_gz = scratch_path + 'field_map.nii'
dwi_topup =  scratch_path + 'dwi_topup'
command = 'topup --imain=' +topup_in_nii +' --datain=' + topup_datain_txt + ' --out=field --fout='+field_map_nii_gz + ' --config=$FSLDIR/src/topup/flirtsch/b02b0.cnf --out='+dwi_topup +' --verbose'
print(command)
os.system(command)
#cannot be done here so echo and do in terminal

#Command:  mrconvert dwi_pad2.mif -import_pe_table dwi_manual_pe_scheme.txt - | mrinfo - -export_pe_eddy applytopup_config.txt applytopup_indices.txt
applytopup_config_txt  = scratch_path + 'applytopup_config.txt'
applytopup_indices_txt = scratch_path + 'applytopup_indices.txt'
dwi_manual_pe_scheme_txt = root + 'dwi_manual_pe_scheme.txt'
command = 'mrconvert '+ dwi_pad2_mif+ ' -import_pe_table '+ dwi_manual_pe_scheme_txt + ' - | mrinfo - -export_pe_eddy ' + applytopup_config_txt + ' '+ applytopup_indices_txt + ' -force'
print(command)
os.system(command)



#Command:  mrconvert dwi_pad2.mif dwi_pad2_pe_0.nii -coord 3 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91 -strides -1,+2,+3,+4 -json_export dwi_pad2_pe_0.json
dwi_pad2_pe_0_nii = scratch_path + 'dwi_pad2_pe_0.nii'
dwi_pad2_pe_0_json = scratch_path + 'dwi_pad2_pe_0.json'
command = 'mrconvert ' + dwi_pad2_mif + ' ' +dwi_pad2_pe_0_nii+' -strides -1,+2,+3,+4 -json_export ' + dwi_pad2_pe_0_json+' -force'
print(command)
os.system(command)


#Command:  applytopup --imain=dwi_pad2_pe_0.nii --datain=applytopup_config.txt --inindex=1 --topup=field --out=dwi_pad2_pe_0_applytopup.nii --method=jac
applytopup_config_txt = scratch_path + 'applytopup_config.txt'
dwi_pad2_pe_0_applytopup_nii = scratch_path+ 'dwi_pad2_pe_0_applytopup.nii'
command = 'applytopup --imain='+dwi_pad2_pe_0_nii+' --datain='+applytopup_config_txt+' --inindex=1 --topup='+ dwi_topup+' --out='+dwi_pad2_pe_0_applytopup_nii  +' --method=jac'
print(command)
os.system(command)
#cannot be done here so echo and do in terminal

#Command:  mrconvert dwi_pad2_pe_0_applytopup.nii.gz dwi_pad2_pe_0_applytopup.mif -json_import dwi_pad2_pe_0.json
dwi_pad2_pe_0_applytopup_nii = dwi_pad2_pe_0_applytopup_nii +'.gz'
dwi_pad2_pe_0_applytopup_mif = scratch_path + 'dwi_pad2_pe_0_applytopup.mif'
command = 'mrconvert ' + dwi_pad2_pe_0_applytopup_nii+ ' ' + dwi_pad2_pe_0_applytopup_mif+  ' -json_import '+ dwi_pad2_pe_0_json +' -force'
print(command)
os.system(command)


#Command:  dwi2mask dwi_pad2_pe_0_applytopup.mif - | maskfilter - dilate - | mrconvert - eddy_mask.nii -datatype float32 -strides -1,+2,+3
eddy_mask_nii = scratch_path +'eddy_mask.nii'
command = 'dwi2mask '+ dwi_pad2_pe_0_applytopup_mif + ' - | maskfilter - dilate - | mrconvert - ' + eddy_mask_nii +' -datatype float32 -strides -1,+2,+3 -force'
print(command)
os.system(command)


#Command:  mrconvert dwi_pad2.mif -import_pe_table dwi_manual_pe_scheme.txt eddy_in.nii -strides -1,+2,+3,+4 -export_grad_fsl bvecs bvals -export_pe_eddy eddy_config.txt eddy_indices.txt
bvecs_eddy  = path_perm + subj+ "_bvecs_eddy.txt"
bvals_eddy = path_perm +subj +"_bvals_eddy.txt"
eddy_config_txt = scratch_path +'eddy_config.txt'
eddy_indices_txt =scratch_path + 'eddy_indices.txt'
eddy_in_nii = scratch_path + 'eddy_in.nii'
command = 'mrconvert '+ dwi_pad2_mif+ ' -import_pe_table ' +dwi_manual_pe_scheme_txt + ' ' + eddy_in_nii + ' -strides -1,+2,+3,+4 -export_grad_fsl '+ bvecs_eddy +  ' ' + bvals_eddy+ ' -export_pe_eddy ' +eddy_config_txt + " " + eddy_indices_txt + ' -force'
print(command)
os.system(command)


#Command:  eddy_openmp --imain=eddy_in.nii --mask=eddy_mask.nii --acqp=eddy_config.txt --index=eddy_indices.txt --bvecs=bvecs --bvals=bvals --topup=field --slm=linear --data_is_shelled --out=dwi_post_eddy --verbose
dwi_post_eddy = scratch_path + "dwi_post_eddy.nii"
command = 'eddy_openmp --imain='+eddy_in_nii+ ' --mask='+eddy_mask_nii + ' --acqp='+eddy_config_txt + ' --index='+ eddy_indices_txt + ' --bvecs=' +bvecs_eddy+' --bvals='+ bvals_eddy+ ' --topup='+dwi_topup+' --slm=linear --data_is_shelled --out='+dwi_post_eddy + ' --verbose'
print(command)
os.system(command)
####does not exist on mac so switching to munin


#Command:  mrconvert dwi_post_eddy.nii.gz result.mif -coord 2 0:76 -strides -1,-2,3,4 -fslgrad dwi_post_eddy.eddy_rotated_bvecs bvals
dwi_post_eddy = dwi_post_eddy + ".gz"
result_mif = scratch_path + subj + '_coreg_result.mif'
dwi_post_eddy_eddy_rotated_bvecs = scratch_path + 'dwi_post_eddy.nii.eddy_rotated_bvecs'
command = 'mrconvert '+ dwi_post_eddy+ ' '+ result_mif+ ' -coord 2 0:'+end_vol+' -strides -1,-2,3,4 -fslgrad ' + dwi_post_eddy_eddy_rotated_bvecs + ' ' +bvals_eddy+ ' -force'
print(command)
os.system(command)


#Command:  mrconvert result.mif ###/ADRC_mrtrix/temp/ADRC0001/ADRC0001_den_preproc.mif
command = 'mrconvert '+ result_mif+ ' '+ den_preproc_mif+ ' -force'
print(command)
os.system(command)




#import nipype.interfaces.mrtrix3 as mrt
#preproc = mrt.DWIPreproc()
#preproc.inputs.in_file = output_denoise
#preproc.inputs.rpe_options = "pair"
#preproc.inputs.in_epi = b0_pair_mif
#preproc.inputs.out_file = den_preproc_mif
#preproc.grad_fsl = (bvec_path , bval_path )
#preproc.inputs.eddy_options = '--slm=linear --data_is_shelled'
#preproc.inputs.ro_time = 40 #YOUR READOUT TIME HERE
#preproc.inputs.pe_dir = 'AP'
#preproc.inputs.nthreads = 8
#preproc.inputs.args = "-force"
#preproc.inputs.export_grad_mrtrix = False
#preproc.inputs.export_grad_fsl = False
#preproc.outputs.out_grad_mrtrix = "###/ADRC_mrtrix/temp/ADRC0001/gradgrad.b"
#print(preproc.cmdline)
#res = preproc.run()

#os.environ['FSLDIR'] = '/usr/local/fsl' # only for spider
#os.system('dwifslpreproc '+output_denoise+ ' ' +den_preproc_mif+' -pe_dir AP -rpe_pair -se_epi '+ b0_pair_mif +' -eddy_options " --slm=linear --data_is_shelled" -nthreads '+ str(nthreads) +' -scratch ' + subj_path+' -fslgrad ' +bvec_path + ' '+ bval_path  +' -force')
#if by echo does not work copy and paste
#os.system('mrview '+den_preproc_mif+' -overlay.load '+out_mif)
#den_preproc_mif = output_denoise # bypass this it takes a lot of time

#os.system('mrview '+den_preproc_mif+' -overlay.load '+out_mif)


#createing mask after bias correction:
den_unbiased_mif = subj_path+subj+'_den_preproc_unbiased.mif'
bias_mif = subj_path+subj+'_bias.mif'
command = 'dwibiascorrect ants '+den_preproc_mif+' '+den_unbiased_mif+ ' -scratch ' +subj_path +' -bias '+ bias_mif + ' -force'
print(command)
os.system(command)
#cannot be done here go on on terminal after echoing and python it
#den_unbiased_mif = den_preproc_mif  # bypassing
#dwibiascorrect ants ###/ADRC_mrtrix_dwifsl/temp/ADRC0001/ADRC0001_den_preproc.mif ###/ADRC_mrtrix_dwifsl/temp/ADRC0001/ADRC0001_den_preproc_unbiased.mif  -bias ###/ADRC_mrtrix_dwifsl/temp/ADRC0001/ADRC0001__bias.mif -scratch ###/ADRC_mrtrix_dwifsl/temp/ADRC0001/ -force


dwi_mif = subj_path+subj+'_dwi.mif'
command = 'dwiextract '+den_unbiased_mif+ ' - -no_bzero | mrmath - mean '+ dwi_mif+' -axis 3 -force'
print(command)
os.system(command)

dwi_nii_gz = path_perm+subj+'_dwi.nii.gz'
command = 'mrconvert ' +dwi_mif + ' '+ dwi_nii_gz + ' -force' 
print(command)
os.system(command)

coreg_bvecs = path_perm+subj+'_coreg_bvecs.txt'
coreg_bvals = path_perm+subj+'_coreg_bvals.txt'
coreg_nii_gz = path_perm+subj+'_coreg.nii.gz'
command = 'mrconvert ' +den_unbiased_mif + ' '+ coreg_nii_gz + ' -export_grad_fsl '+ coreg_bvecs + ' ' +coreg_bvals + ' -force'
print(command)
os.system(command)





mask_mif  =  path_perm+subj+'_mask.mif'
command = 'dwi2mask '+den_unbiased_mif+  ' '+ mask_mif + ' -force'
print(command)
os.system(command) 

#os.system('mrview '+fa_mif + ' -overlay.load '+ mask_mif )  
mask_mrtrix_nii = path_perm +subj+'_mask.nii.gz'
command = 'mrconvert ' +mask_mif+ ' '+mask_mrtrix_nii + ' -force' 
print(command)
os.system(command)





########### making a mask out of labels

#label_path_orig = orig_subj_path +subj+'_labels.nii.gz'
#label_path = path_perm +subj+'_labels.nii.gz'
#os.system("/Applications/Convert3DGUI.app/Contents/bin/c3d "+label_path_orig+" -orient RAS -o "+label_path)
#label_path = label_path_orig
#mask_output = subj_path +subj+'_mask_of_label.nii.gz'
#label_nii=nib.load(label_path)
#mask_labels_data = label_nii.get_fdata()
#mask_labels = np.unique(mask_labels_data)
#mask_labels=np.delete(mask_labels, 0)
#mask_of_label =label_nii.get_fdata()*0



#path_atlas_legend = root+ 'IIT/IITmean_RPI_index.xlsx'
#legend  = pd.read_excel(path_atlas_legend)
#index_csf = legend [ 'Subdivisions_7' ] == '8_CSF'
#index_wm = legend [ 'Subdivisions_7' ] == '7_whitematter'

#vol_index_csf = legend[index_csf]
#vol_index_csf = vol_index_csf['index2']

#mask_labels_no_csf= set(mask_labels)-set(vol_index_csf )
#mask_labels_no_csf = mask_labels
#for vol in mask_labels_no_csf: mask_of_label[  mask_labels_data == int(vol)] = int(1)
#mask_of_label= mask_of_label.astype(int)

#file_result= nib.Nifti1Image(mask_of_label, label_nii.affine, label_nii.header)
#nib.save(file_result,mask_output  )  
#mask_labels_mif   = subj_path +subj+'mask_of_label.mif'+index_gz
#os.system('mrconvert '+mask_output+ ' ' +mask_labels_mif+ ' -datatype uint16 -force')
#os.system('mrview '+fa_mif + ' -overlay.load '+ mask_labels_mif ) 


#new_bval_path = path_perm+subj+'_new_bvals.txt' 
#new_bvec_path = path_perm+subj+'_new_bvecs.txt' 
#os.system('dwigradcheck ' + out_mif +  ' -fslgrad '+bvec_path_orig + ' '+ bval_path_orig +' -mask '+ mask_mif + ' -number 100000 -export_grad_fsl '+ new_bvec_path + ' '  + new_bval_path  +  ' -force' )
#bvec_temp=np.loadtxt(new_bvec_path)

    
#Estimating the Basis Functions:
wm_txt =   subj_path+subj+'_wm.txt' 
gm_txt =  subj_path+subj+'_gm.txt' 
csf_txt = subj_path+subj+'_csf.txt'
voxels_mif =  subj_path+subj+'_voxels.mif'+index_gz
command = 'dwi2response dhollander '+den_unbiased_mif+ ' ' +wm_txt+ ' ' + gm_txt + ' ' + csf_txt + ' -voxels ' + voxels_mif+' -mask '+ mask_mif + ' -scratch ' +subj_path + ' -fslgrad ' +bvec_path + ' '+ bval_path   +'  -force' 
print(command)
os.system(command)

#not on spider but on terminal
#Viewing the Basis Functions:
#os.system('mrview '+den_unbiased_mif+ ' -overlay.load '+ voxels_mif)
#os.system('shview '+wm_txt)
#os.system('shview '+gm_txt)
#os.system('shview '+csf_txt)

#Applying the basis functions to the diffusion data:
wmfod_mif =  subj_path+subj+'_wmfod.mif'+index_gz
gmfod_mif = subj_path+subj+'_gmfod.mif'+index_gz
csffod_mif = subj_path+subj+'_csffod.mif'+index_gz

#os.system('dwi2fod msmt_csd ' +den_unbiased_mif+ ' -mask '+mask_mif+ ' ' +wm_txt+ ' ' + wmfod_mif+ ' ' +gm_txt+ ' ' + gmfod_mif+ ' ' +csf_txt+ ' ' + csffod_mif + ' -force' )
command = 'dwi2fod msmt_csd ' +den_unbiased_mif+ ' -mask '+mask_mif+ ' ' +wm_txt+ ' ' + wmfod_mif+ ' ' +gm_txt+ ' ' + gmfod_mif+ ' ' +csf_txt+ ' ' + csffod_mif + ' -force' 
print(command)
os.system(command)

#combine to single image to view them
#Concatenating the FODs:
vf_mif =   subj_path+subj+'_vf.mif' 
command = 'mrconvert -coord 3 0 ' +wmfod_mif+ ' -| mrcat '+csffod_mif+ ' ' +gmfod_mif+ ' - ' + vf_mif+' -force' 
print(command)
os.system(command)
#os.system('mrconvert -coord 3 0 ' +wmfod_mif+ ' -| mrcat ' +gmfod_mif+ ' - ' + vf_mif+' -force' ) # without csf

#Viewing the FODs:
#os.system('mrview ' +fa_mif+ ' -odf.load_sh '+wmfod_mif )

#Normalizing the FODs:
wmfod_norm_mif =  subj_path+subj+'_wmfod_norm.mif'+index_gz
gmfod_norm_mif = subj_path+subj+'_gmfod_norm.mif'
csffod_norm_mif = subj_path+subj+'_csffod_norm.mif'  
command = 'mtnormalise ' +wmfod_mif+ ' '+wmfod_norm_mif+ ' ' +gmfod_mif+ ' '+gmfod_norm_mif + ' ' +csffod_mif+ ' '+csffod_norm_mif +' -mask ' + mask_mif + '  -force'
print(command)
os.system(command)







#making fa and Kurt:
    
dt_mif = path_perm+subj+'_dt.mif'+index_gz
fa_mif = path_perm+subj+'_fa.mif'+index_gz
dk_mif = path_perm+subj+'_dk.mif'+index_gz
mk_mif = path_perm+subj+'_mk.mif'+index_gz
md_mif = path_perm+subj+'_md.mif'+index_gz
ad_mif = path_perm+subj+'_ad.mif'+index_gz
rd_mif = path_perm+subj+'_rd.mif'+index_gz

#output_denoise = '###/Desktop/Feb23/mrtrix_pipeline/temp/N59141/N59141_subjspace_dwi_copy.mif.gz'#



if np.unique(new_bval).shape[0] > 2 :
    os.system('dwi2tensor ' + den_unbiased_mif + ' ' + dt_mif + ' -dkt ' +  dk_mif +' -fslgrad ' +  bvec_path + ' ' + bval_path + ' -force'  )
    os.system('tensor2metric  -fa ' + fa_mif  + ' '+ dt_mif + ' -adc '  + md_mif+' -ad '  + ad_mif + ' -rd '  + rd_mif   + ' -force' ) 

    #os.system('mrview '+ fa_mif) #inspect residual
else: 
    os.system('dwi2tensor ' + den_unbiased_mif + ' ' + dt_mif  +' -fslgrad ' +  bvec_path + ' ' + bval_path + ' -force'  )
    os.system('tensor2metric  -fa ' + fa_mif  + ' '+ dt_mif + ' -force' ) 
    os.system('tensor2metric  -rd ' + rd_mif  + ' '+ dt_mif + ' -force' ) # if doesn't work take this out :(
    os.system('tensor2metric  -ad ' + ad_mif  + ' '+ dt_mif + ' -force' ) # if doesn't work take this out :(
    os.system('tensor2metric  -adc ' + md_mif  + ' '+ dt_mif + ' -force' ) # if doesn't work take this out :(
    #os.system('mrview '+ fa_mif) #inspect residual

fa_nii = path_perm+subj+'_fa.nii'+index_gz
command = 'mrconvert ' +fa_mif+ ' '+fa_nii + ' -force' 
print(command)
os.system(command) 


#os.system('mrinfo '+den_unbiased_mif) #info





#os.system('mrinfo -size '+out_mif+" | awk '{print $4}' ") #print the 4th dimension 
#os.system("awk '{print NF; exit}' "+bvec_path) #print the lenght of bvec : must match the previous number
#os.system("awk '{print NF; exit}' "+bval_path) #print the lenght of bval : must match the previous number
#os.system("mrview "+ out_mif) #viewer 










#Viewing the normalise FODs:
#os.system('mrview ' +fa_mif+ ' -odf.load_sh '+wmfod_norm_mif )

###############







#boundry of white and grey matter used for streamlines
#Converting the anatomical image to MRtrix format:

#Segmenting the anatomical image with FSL's FAST to 5 different classes
#fivett_nocoreg_mif  = subj_path+subj+'5tt_nocoreg.mif'
#os.system('5ttgen fsl '  +T1_mif+ ' '+fivett_nocoreg_mif + ' -force')
#cannot be done here go on on terminal after echoing and python it

#os.system('mrview ' +fivett_nocoreg_mif  )

#mean_b0_mif = subj_path+subj+'_mean_b0.mif'
#Extracting the b0 images: for Coregistering the anatomical and diffusion datasets:
#os.system('dwiextract '+ den_unbiased_mif+' - -bzero | mrmath - mean '+ mean_b0_mif +' -axis 3 -force')

#Converting the b0 and 5tt images bc we wanna use fsl this part and fsl does not accept mif:
#mean_b0_nii_gz    = subj_path+subj+'_mean_b0.nii.gz'
#fivett_nocoreg_nii_gz = subj_path+subj+'_5tt_nocoreg.nii.gz'
#os.system('mrconvert ' +mean_b0_mif + ' '+ mean_b0_nii_gz + ' -force')
#os.system('mrconvert ' + fivett_nocoreg_mif + ' ' + fivett_nocoreg_nii_gz + ' -force')

# now Extracting the grey matter segmentation with fsl:
#fivett_vol0_nii_gz    =  subj_path+subj+'_5tt_vol0.nii.gz'
#os.system('fslroi '+ fivett_nocoreg_nii_gz+ ' '+ fivett_vol0_nii_gz + ' 0 1')
#if not working here, works on terminal after echoing and copy\pasting


#Coregistering the anatomical and diffusion datasets:
#diff2struct_fsl_mat =subj_path+subj+'_diff2struct_fsl.mat'
#os.system('flirt -in '+ mean_b0_nii_gz + ' -ref ' + fivett_vol0_nii_gz + ' -interp nearestneighbour -dof 6 -omat ' + diff2struct_fsl_mat  )
#if not working here, works on terminal after echoing and copy\pasting

#Converting the transformation matrix to MRtrix format: 
#diff2struct_mrtrix_txt = subj_path+subj+'_diff2struct_mrtrix.txt'
#os.system('transformconvert ' + diff2struct_fsl_mat + ' '+ mean_b0_nii_gz+ ' '+ fivett_nocoreg_nii_gz + ' flirt_import '+ diff2struct_mrtrix_txt + ' -force' )

#Applying the transformation matrix to the non-coregistered segmentation data:
#using the iverse transfomration coregsiter anatomiacl to dwi
#fivett_coreg_mif   = subj_path+subj+'_fivett_coreg.mif'
#os.system('mrtransform ' + fivett_nocoreg_mif + ' -linear ' + diff2struct_mrtrix_txt + ' -inverse '+ fivett_coreg_mif + ' -force')

#Viewing the coregistration in mrview:
#os.system( 'mrview '+ den_unbiased_mif +' -overlay.load ' + fivett_nocoreg_mif + ' -overlay.colourmap 2 -overlay.load ' + fivett_coreg_mif + ' -overlay.colourmap 1 ')
#os.system( 'mrview '+ T1_mif +' -overlay.load ' + fivett_nocoreg_mif + ' -overlay.colourmap 2 -overlay.load ' + fivett_coreg_mif + ' -overlay.colourmap 1 ')
#fivett_coreg_mif = fivett_nocoreg_mif

#Creating the grey matter / white matter boundary: seed boundery bc they're used to create seeds for streamlines
#gmwmSeed_coreg_mif = subj_path+subj+'_gmwmSeed_coreg.mif'
#os.system( '5tt2gmwmi ' +  fivett_coreg_mif+ ' '+ gmwmSeed_coreg_mif + ' -force')


#os.system('mrview ' +fivett_nocoreg_mif + ' -overlay.load ' + mean_b0_mif )

# Viewing the GM/WM boundary:
#os.system('mrview ' + den_unbiased_mif + ' -overlay.load ' + gmwmSeed_coreg_mif)

shutil.rmtree(scratch_path)



gmwmSeed_coreg_mif  = mask_mif  


###############


####read to creat streamlines
#Creating streamlines with tckgen: be carefull about number of threads on server
tracks_10M_tck  = subj_path +subj+'_tracks_10M.tck' 

#os.system('tckgen -act ' + fivett_coreg_mif + '  -backtrack -seed_gmwmi '+ gmwmSeed_coreg_mif + ' -maxlength 250 -cutoff 0.06 -select 10000000 ' + wmfod_norm_mif + ' ' + tracks_10M_tck + ' -force')
#seconds1 = time.time()
command ='tckgen -backtrack -seed_image '+ gmwmSeed_coreg_mif + ' -maxlength 250 -cutoff 0.1 -select 10000000 ' + wmfod_norm_mif + ' ' + tracks_10M_tck + ' -force'
print(command)
os.system( command )


#os.system('tckgen -backtrack -seed_image '+ gmwmSeed_coreg_mif + ' -maxlength 1000 -cutoff 0.3 -select 50k ' + wmfod_norm_mif + ' ' + tracks_10M_tck + ' -force')
#seconds2 = time.time()
#(seconds2 - seconds1)/360 # a million track in hippo takes 12.6 mins


#Extracting a subset of tracks:
smallerTracks = path_perm+subj+'_smallerTracks2mill.tck'
#os.system('echo tckedit '+ tracks_10M_tck + ' -number 2000000 -minlength 0.1 ' + smallerTracks + ' -force')

command = 'tckedit '+ tracks_10M_tck + ' -number 2000000 ' + smallerTracks + ' -force'
print(command)
os.system(command)

#os.system('tckedit '+ tracks_10M_tck + ' -number 2000000 -minlength 2 ' + smallerTracks + ' -force')
#os.system('mrview ' + den_unbiased_mif + ' -tractography.load '+ smallerTracks)
#os.system('mrview ' + den_unbiased_mif + ' -tractography.load '+ smallerTracks)







################################# upto here
'''

#Sifting the tracks with tcksift2: bc some wm tracks are over or underfitted
sift_mu_txt = subj_path+subj+'_sift_mu.txt'
sift_coeffs_txt = subj_path+subj+'_sift_coeffs.txt'
sift_1M_txt = subj_path+subj+'_sift_1M.txt'

os.system('echo tcksift2  -out_mu '+ sift_mu_txt + ' -out_coeffs ' + sift_coeffs_txt + ' ' + smallerTracks + ' ' + wmfod_norm_mif+ ' ' + sift_1M_txt  + ' -force')
os.system('tcksift2  -out_mu '+ sift_mu_txt + ' -out_coeffs ' + sift_coeffs_txt + ' ' + smallerTracks + ' ' + wmfod_norm_mif+ ' ' + sift_1M_txt  + ' -force')

#####connectome
##Running recon-all:

#os.system("SUBJECTS_DIR=`pwd`")
#sub_recon = subj_path+subj+'_recon3'
#os.system('recon-all -i '+ T1 +' -s '+ sub_recon +' -all -force')
# cant run here so do on command line


#Converting the labels:
#parcels_mif = subj_path+subj+'_parcels.mif'
#os.system('labelconvert '+ ' ###/sub-CON02_recon3/mri/aparc+aseg.mgz' + ' /Applications/freesurfer/7.3.2/FreeSurferColorLUT.txt ' +  '###/opt/anaconda3/pkgs/mrtrix3-3.0.3-ha664bf1_0/share/mrtrix3/labelconvert/fs_default.txt '+ parcels_mif)



#Coregistering the parcellation:
#diff2struct_mrtrix_txt = subj_path+subj+'_diff2struct_mrtrix.txt'
#parcels_coreg_mif = subj_path+subj+'_parcels_coreg.mif'
#os.system('mrtransform '+parcels_mif + ' -interp nearest -linear ' + diff2struct_mrtrix_txt + ' -inverse -datatype uint32 ' + parcels_coreg_mif )



#convert subj labels to mif

labels_data = label_nii.get_fdata()
labels = np.unique(labels_data)
labels=np.delete(labels, 0)
label_nii_order = labels_data*0.0

#sum(legend['index2'] == labels)
for i in labels: 
    leg_index =  np.where(legend['index2'] == i )
    leg_index = leg_index [0][0]
    ordered_num = legend['index'][leg_index]
    label3d_index = np.where( labels_data == i )
    label_nii_order [ label3d_index]  = ordered_num


file_result= nib.Nifti1Image(label_nii_order, label_nii.affine, label_nii.header)
new_label  = path_perm +subj+'_new_label.nii.gz'
nib.save(file_result, new_label ) 

parcels_mif = subj_path+subj+'_parcels.mif'+index_gz
#new_label = label_path
os.system('mrconvert '+new_label+ ' ' +parcels_mif + ' -force' )


#os.system('mrview '+ fa_mif + ' -overlay.load '+ new_label) 



#Creating the connectome without coregistration:
### connectome folders :
    
conn_folder = root + 'connectome/'
if not os.path.isdir(conn_folder) : os.mkdir(conn_folder)
    
    

distances_csv = conn_folder +subj+'_distances.csv'
os.system('tck2connectome ' + smallerTracks + ' ' + parcels_mif+ ' ' + distances_csv + ' -zero_diagonal -symmetric -scale_length -stat_edge  mean' + ' -force')
mean_FA_per_streamline =  subj_path+subj+'_per_strmline_mean_FA.csv'
mean_FA_connectome =  conn_folder+subj+'_mean_FA_connectome.csv'
os.system('tcksample '+ smallerTracks+ ' '+ fa_mif + ' ' + mean_FA_per_streamline + ' -stat_tck mean ' + ' -force')
os.system('tck2connectome '+ smallerTracks + ' ' + parcels_mif + ' '+ mean_FA_connectome + ' -zero_diagonal -symmetric -scale_file ' + mean_FA_per_streamline + ' -stat_edge mean '+ ' -force')

  




parcels_csv = conn_folder+subj+'_conn_sift_node.csv'
assignments_parcels_csv = path_perm +subj+'_assignments_con_sift_node.csv'
os.system('tck2connectome -symmetric -zero_diagonal -scale_invnodevol -tck_weights_in '+ sift_1M_txt+ ' '+ smallerTracks + ' '+ parcels_mif + ' '+ parcels_csv + ' -out_assignment ' + assignments_parcels_csv + ' -force')


parcels_csv_2 = conn_folder+subj+'_conn_plain.csv'
assignments_parcels_csv2 = path_perm +subj+'_assignments_con_plain.csv'
os.system('tck2connectome -symmetric -zero_diagonal '+ smallerTracks + ' '+ parcels_mif + ' '+ parcels_csv_2 + ' -out_assignment ' + assignments_parcels_csv2 + ' -force')

parcels_csv_3 = conn_folder+subj+'_conn_sift.csv'
assignments_parcels_csv3 = path_perm +subj+'_assignments_con_sift.csv'
os.system('tck2connectome -symmetric -zero_diagonal -tck_weights_in '+ sift_1M_txt+ ' '+ smallerTracks + ' '+ parcels_mif + ' '+ parcels_csv_3 + ' -out_assignment ' + assignments_parcels_csv3 + ' -force')
'''

shutil.rmtree(subj_path)


#scale_invnodevol scale connectome by the inverse of size of each node
#tck_weights_in weight each connectivity by sift
#out assignment helo converting connectome to tracks

#Viewing the connectome in Matlab:

#connectome = importdata('sub-CON02_parcels.csv');
#imagesc(connectome, [0 1])

#Viewing the lookup labels:

