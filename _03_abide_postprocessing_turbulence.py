# _03_abide_postprocessing_turbulence.py
"""
_03_abide_postprocessing_turbulence.py

Run postprocessing Turbulence pipeline on ABIDE data.

Example usage:

# Run on server and generate scripts
python _03_abide_postprocessing_turbulence.py server script

# Run on server and run batches
python _03_abide_postprocessing_turbulence.py server batch

--- Written by bfuhr and mvlombardo 16.12.2024
"""

# import modules ---------------------------------------------------------------
import os
import sys
import fnmatch
import glob
import time
import nibabel as nib
import pandas as pd
import numpy as np


# Parse input arguments --------------------------------------------------------
# First argument tells us whether we are running it on the lab's server or our own laptop
if sys.argv[1]=="server":
    run_on_server = True
elif sys.argv[1]=="laptop":
    run_on_server = False

# Second argument tells if you are running to generate scripts only or to run batches as well
if sys.argv[2]=="script":
    script_only = True
elif sys.argv[2]=="batch":
    script_only = False


# Directories ------------------------------------------------------------------
# specify my rootdir that all other directories follow from
if run_on_server:
    rootdir = "/media/DATA/RAW/abide"
    server_flag = "true"
    master_maskfile = "/usr/share/fsl/5.0/data/standard/MNI152_T1_1mm_brain_mask.nii.gz"
    fsl_matlab_path = "/usr/share/fsl/5.0/etc/matlab"
else:
    rootdir = "/Users/mlombardo/Dropbox/data/abide"
    server_flag = "false"
    master_maskfile = "/usr/local/fsl/data/standard/MNI152_T1_1mm_brain_mask.nii.gz"
    fsl_matlab_path = "/Users/mlombardo/fsl/etc/matlab"

# other directories I need
datadir = "%s/data" % rootdir
preprocdir = "%s/preproc" % datadir
postprocdir = "%s/postproc" % datadir
codedir = "%s/code" % rootdir
parcdir = "%s/parc" % rootdir
os.system("mkdir -p %s" % postprocdir)

# specify the metric you are computing (used later for naming directories and files)
metric_type = "turbulence"

# parcellation file to use
parcfile = "%s/Schaefer2018_1000Parcels_7Networks_order_FSLMNI152_2mm.nii.gz" % parcdir
cogfile =  "%s/schaefer_cog.mat" % parcdir
labelsfile =  "%s/schaefer_labels.mat" % parcdir

# load in MRI params and preproc QC file to get sublist to loop over
mriparamsfile = "%s/pheno/mriparams_summary_abide_I_abide_II.csv" % datadir
mri_params_data = pd.read_csv(mriparamsfile)

qcfile = "%s/pheno/qc_preproc_abideI_abideII.csv" % datadir
qc_data = pd.read_csv(qcfile)


# Fourth argument specifies number of computational threads
if len(sys.argv)>3:
    # specify number of computational threads for use n MATLAB
    maxNumCompThreads = int(sys.argv[3])
else:
    # specify number of computational threads for use n MATLAB
    maxNumCompThreads = 1


# generate my list of subjects from subjects already present in the raw data directory
sublist = os.listdir(preprocdir)

if run_on_server==False:
    # find .DS_Store if it exists and remove it from sublist
    pattern2find = ".DS_Store"
    matches = fnmatch.filter(sublist, pattern2find)
    if matches:
        sublist.remove(pattern2find)

if "tmp" in sublist:
    sublist.remove('tmp')


# loop over subjects
for sub in sublist:

    # subject-specific paths in raw and preproc directories
    subpath = "%s/%s" % (preprocdir,sub)
    postproc_subpath = "%s/%s" % (postprocdir,sub)

    # create output directory in the subject's postproc directory
    outdir = "%s/%s" % (postproc_subpath, metric_type)
    os.system("mkdir -p %s" % outdir)

    # start my bash script
    script = []
    script.append("#!/bin/bash")
    script.append("")
    script.append("start=`date +%s`")
    script.append("")

    # if you're running this on the server, 
    # you need to source these shell scripts 
    # to call AFNI and FSL programs
    if run_on_server:
        script.append("source /etc/afni/afni.sh")
        script.append("source /etc/fsl/fsl.sh")
        script.append("")

    # Step 1:  resample preproc data to 2mm (the voxel size of the Schaefer)
    resamp_ppfile = "%s/tmpdata2mm.nii.gz" % postproc_subpath
    orig_ppfile = "%s/rest_pp.nii.gz" % subpath
    resample_code2use = "3dresample -master %s -prefix %s -input %s -rmode Cu" % (parcfile, resamp_ppfile, orig_ppfile)
    script.append("# Resampling preproc data")
    script.append(resample_code2use)
    script.append("")

    # Step 2:  find the TR for the subject
    result = qc_data[qc_data['subid'] == int(sub)]
    site = result.site.tolist()[0]

    result = mri_params_data[mri_params_data['site'] == site]
    tr2use = result.tr_sec.tolist()[0]

    # Step 3:  compute turbulence
    parc_command = "ts = parcellate('%s', '%s')" % (resamp_ppfile, parcfile)
    turb_command = "[output] = turbulence('%s', ts, %s, '%s', '%s', '%s')" % (sub,tr2use, cogfile, labelsfile,outdir)
    matlab_command2run = "cd('%s'); maxNumCompThreads(%d); addpath %s; tic; %s; %s; toc; exit;" % (codedir, maxNumCompThreads, fsl_matlab_path, parc_command, turb_command)

    script.append("# Run turbulence pipeline in MATLAB")
    script.append('matlab -nodesktop -nosplash -r "%s"' % matlab_command2run)
    script.append("")

    # Step 4:  write the shell script to disk
    bash_script_name = "_03_postprocessing_%s_rest_%s.sh" % (sub,metric_type)
    out_log_name = "%s_postproc_%s_out.log" % (sub,metric_type)
    error_log_name = "%s_postproc_%s_error.log" % (sub,metric_type)

    script.append("# Delete the temporary resampled preproc data file")
    script.append("rm -Rf %s" % resamp_ppfile)

    script.append("")
    script.append("end=`date +%s`")
    script.append("runtime=$((end-start))")
    script.append("echo Time to complete = $runtime secs")
    script.append("")

    print("Saving %s postprocessing bash script for %s rest ..." %(metric_type,sub))
    fname2save = "%s/%s" % (outdir,bash_script_name)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()

    os.system("mkdir -p %s/_03_postprocessing_batches" % codedir)
    fname2save = "%s/_03_postprocessing_batches/%s" % (codedir,bash_script_name)
    fname = open(fname2save,"w")
    fname.write("\n".join(script)+"\n")
    fname.close()

    # if you're not just generating scripts and want to run the bash scripts
    # then do this below...
    if not script_only:

        if run_on_server:

            print("Inserting batch into the queue ...")

            # options for sbatch
            outlog = "%s/%s" % (outdir,out_log_name)
            errlog = "%s/%s" % (outdir,error_log_name)
            sbatch_opts = "-J %s --tasks-per-node=1 --cpus-per-task=1 --time=36:00:00 --no-requeue --mem=32g --output=%s --error=%s" % (sub,outlog,errlog)

            # run batch
            os.system("sbatch %s %s" % (sbatch_opts,fname2save))

        else:

            # run bash script
            os.system("bash %s" % fname2save)


