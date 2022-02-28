"""Usage: example_run.py -G VALUE -n VALUE -i FILE -o FILE -p VALUE -a VALUE -d VALUE [-l LENGTH] [-s SEED]

-G VALUE        global coupling scaling
-n VALUE        noise sigma
-i FILE         input file (npz)
-o FILE         output file (npz)
-p FILE         patient location or code
-a VALUE        distorsion parameter for SC weights matrix
-d VALUE        value of dt
-l LENGTH       simulation total length in ms [default: 1000]
-s SEED         random number generator seed [default: 42]
-h --help       show this
"""

from docopt import docopt
import os
import numpy as np
from src import data, analysis, simulation, viz  # Import analysis for fcd and clustering
from tvb.simulator.lab import *
import time
from scipy import signal 
from scipy import stats
import re

import matplotlib.pyplot as plt
from matplotlib import rcParams
rcParams['figure.figsize'] = 15,6 # Used for some figures, but check for what
from tvb.datatypes.time_series import * # Import TVB data
import socket # Import for Device namae
from datetime import datetime # Import for datetime
from scipy.stats import zscore


if __name__ == '__main__':
    args      = docopt(__doc__)
    G         = float(args["-G"])
    nsigma    = float(args["-n"])
    out_path  = args["-o"]
    in_path   = args["-i"]
    sim_len   = int(args["-l"])
    seed      = int(args["-s"])
    mysubj    = args["-p"]
    alpha1    = float(args["-a"])
    dt        = float(args["-d"])
    
    sys.stdout = open(f'{out_path}.log', 'w')
    # sys.stderr = sys.stdout

    t0        = time.time()

    try:

        jul_interface             = data.Julich_interface()

        Bold_input_archive,Bold_input_time,input_data = jul_interface.collect_dataset(in_path)

        try:
            Bold_input_data = Bold_input_archive[:,0,:,0]
        except:
            print('CORRUPTED BOLD ARCHIVE')
            Bold_input_data = []
            pass

        jul                                       = data.Julich()
        subjs                                     = jul.list_subjects()
        subj_age,gender,education,subj_ID,_,_,_,_ = jul.metadata()

        SUBJ_TARG     = [subj_loc for subj_loc in range(len(subj_ID)) if mysubj in subj_ID[subj_loc] ][0]
        myage         = subj_age[SUBJ_TARG]

        _, weights    = jul.load_subject_sc_100(mysubj)
        NHALF         = int(weights.shape[0]/2)

        print(mysubj,flush=True)
        print(subj_ID[SUBJ_TARG],flush=True)
        print(myage,flush=True)
        print(str(G),flush=True)

        flag_failure = 0

        if len(Bold_input_data)==0:
            print('Empty Input Data',flush=True)
            flag_failure = 1
        elif np.isnan(Bold_input_data).any():
            print('NaN in the Data',flush=True)
            flag_failure = 1
        elif np.sum(Bold_input_data)==0:
            print('Zero Data',flush=True)
            flag_failure = 1 
        
        if flag_failure:

            print('No Bold data in Input or NaN in Input',flush=True)
            print('Dt Recalibration',flush=True)
            
            print(weights,flush=True)

            MAT_WEIGHT    = weights
            magic_number  = 124538.470647693

            ## TRANSFORM SC MATRIX

            weights_symm  = MAT_WEIGHT/magic_number

            # SC AGEING

            mask_inter                              = np.zeros(weights.shape)
            mask_inter[0:NHALF,NHALF:NHALF*2]       = 1
            mask_inter[NHALF:NHALF*2,0:NHALF]       = 1

            weights_symm                            = weights_symm - alpha1*weights_symm*mask_inter

            print(dt,flush=True)
            
            conn               = connectivity.Connectivity(
                weights=weights_symm,
                region_labels=np.array(np.zeros(np.shape(weights_symm)[0]),dtype='<U128'),
                tract_lengths=np.zeros(np.shape(weights_symm)),
                areas =np.zeros(np.shape(weights_symm)[0]),
                speed=np.array(np.Inf,dtype=np.float),
                centres = np.zeros(np.shape(weights_symm)[0]))
            
            print('weight:',weights_symm,flush=True) 
            print('shape :',np.shape(conn.weights),flush=True)
            print('conn :', conn.weights,flush=True)

            mpr = models.MontbrioPazoRoxin(
                eta   = np.r_[-4.6],
                J     = np.r_[14.5],
                Delta = np.r_[0.7],
                tau   = np.r_[1.0],
            )
            
            sim = simulator.Simulator(
                model=mpr,
            connectivity=conn,
                coupling=coupling.Scaling(
                a=np.r_[G]
                ),
            conduction_speed=np.Inf,
                integrator=integrators.HeunStochastic(
                dt=dt,
                noise=noise.Additive(
                    nsig=np.r_[nsigma, nsigma*2],
                    noise_seed=seed
                )
            ),
            monitors=[
                monitors.TemporalAverage(period=0.1),
                ]
            )

            sim.configure()

            try:

                try:

                    (TemporalAverage_time, TemporalAverage_data), = simulation.run_nbMPR_backend(sim, simulation_length=sim_len)
                    TemporalAverage_time *= 10 # rescale time
                    
                    flag_failure = 0

                    if len(TemporalAverage_data)==0:
                        print('Empty Input Data',flush=True)
                        flag_failure = 1
                    elif np.isnan(TemporalAverage_data).any():
                        print('NaN in the Data',flush=True)
                        flag_failure = 1
                    elif np.sum(TemporalAverage_data)==0:
                        print('Zero Data',flush=True)
                        flag_failure = 1 

                    print('BOLD conversion')

                    if flag_failure==0:
                        print('Clean TAVG Data',flush=True)
                        Bold_time, Bold_data = simulation.tavg_to_bold(TemporalAverage_time, TemporalAverage_data, tavg_period=1., connectivity=sim.connectivity, svar=0, decimate=2000)
                    else:
                        print('Corrupted TAVG Data',flush=True)
                        Bold_time = []
                        Bold_data = []

                    print(Bold_data)

                    bold_sig       = Bold_data[:,0,:,0]

                except:

                    print('SIMULATION ERROR')

                    Bold_data      = np.zeros([int(sim_len/2000),weights.shape[0]])
                    Bold_time      = np.zeros([int(sim_len/2000),1])
                    bold_sig       = Bold_data

                    pass

                
                print('FCD Biomarker')

                win_FCD_vect  = [10e3,20e3,30e3,40e3,50e3]
                FCD_SUM_vect  = []
                FCD_MEAN_vect = []
                FCD_VAR_vect  = []
                FCD_OSC_vect  = []

                FCD_SUM_vect  = []
                FCD_MEAN_vect = []
                FCD_VAR_vect  = []
                FCD_OSC_vect  = []

                FCD_SUM_INTER_vect  = []
                FCD_MEAN_INTER_vect = []
                FCD_VAR_INTER_vect  = []
                FCD_OSC_INTER_vect  = []

                FCD_SUM_OV_vect  = []
                FCD_MEAN_OV_vect = []
                FCD_VAR_OV_vect  = []
                FCD_OSC_OV_vect  = []

                FCD_SUM_OV_INTER_vect  = []
                FCD_MEAN_OV_INTER_vect = []
                FCD_VAR_OV_INTER_vect  = []
                FCD_OSC_OV_INTER_vect  = []

                try:
                    
                    for win_idx in range(len(win_FCD_vect)):
                        win_FCD      = win_FCD_vect[win_idx]
                        transient    = int(5e3/2000)
                        FCD, fc_stack, speed_fcd      = analysis.compute_fcd(bold_sig[transient:,:], win_len=int(win_FCD/2000), win_sp=1)
                        fcd_inter, fc_stack_inter, _  = analysis.compute_fcd_filt(bold_sig[transient:,:],mask_inter,win_len=int(win_FCD/2000),win_sp=1)

                        FCD_OV, fc_stack_ov, _      = analysis.compute_fcd(bold_sig[transient:,:], win_len=int(win_FCD/2000), win_sp=int(win_FCD/2000))
                        FCD_OV_INTER, fc_stack_ov_inter, _      = analysis.compute_fcd_filt(bold_sig[transient:,:],mask_inter,win_len=int(win_FCD/2000),win_sp=int(win_FCD/2000))
                        print('FCD')

                        FCD_TRIU      = np.triu(FCD, k=1)
                        FCD_TRIU_OV   = np.triu(FCD, k=int(win_FCD/2000))

                        FCD_INTER_TRIU      = np.triu(fcd_inter, k=1)
                        FCD_INTER_TRIU_OV   = np.triu(fcd_inter, k=int(win_FCD/2000))

                        print('FCD_TRIU')

                        FCD_SUM       = sum(sum(FCD_TRIU))
                        FCD_MEAN      = np.mean(FCD_TRIU)
                        FCD_VAR       = np.var(FCD_TRIU)
                        FCD_OSC        = np.std(fc_stack)
                        FCD_OSC_INTER  = np.std(fc_stack_inter)

                        FCD_SUM_OV           = sum(sum(FCD_TRIU_OV))
                        FCD_MEAN_OV          = np.mean(FCD_TRIU_OV)
                        FCD_VAR_OV           = np.var(FCD_TRIU_OV)
                        FCD_OSC_OV           = np.std(fc_stack_ov)
                        FCD_OSC_OV_INTER     = np.std(fc_stack_ov_inter)

                        FCD_SUM_INTER        = sum(sum(FCD_INTER_TRIU))
                        FCD_MEAN_INTER       = np.mean(FCD_INTER_TRIU)
                        FCD_VAR_INTER        = np.var(FCD_INTER_TRIU)

                        FCD_SUM_OV_INTER     = sum(sum(FCD_INTER_TRIU_OV))
                        FCD_MEAN_OV_INTER    = np.mean(FCD_INTER_TRIU_OV)
                        FCD_VAR_OV_INTER     = np.var(FCD_INTER_TRIU_OV)

                        FCD_SUM_vect  +=[FCD_SUM]
                        FCD_MEAN_vect +=[FCD_MEAN]
                        FCD_VAR_vect  +=[FCD_VAR]
                        FCD_OSC_vect  +=[FCD_OSC]

                        FCD_SUM_INTER_vect  +=[FCD_SUM_INTER]
                        FCD_MEAN_INTER_vect +=[FCD_MEAN_INTER]
                        FCD_VAR_INTER_vect  +=[FCD_VAR_INTER]
                        FCD_OSC_INTER_vect  +=[FCD_OSC_INTER]

                        FCD_SUM_OV_vect  +=[FCD_SUM_OV]
                        FCD_MEAN_OV_vect +=[FCD_MEAN_OV]
                        FCD_VAR_OV_vect  +=[FCD_VAR_OV]
                        FCD_OSC_OV_vect  +=[FCD_OSC_OV]

                        FCD_SUM_OV_INTER_vect  +=[FCD_SUM_OV_INTER]
                        FCD_MEAN_OV_INTER_vect +=[FCD_MEAN_OV_INTER]
                        FCD_VAR_OV_INTER_vect  +=[FCD_VAR_OV_INTER]
                        FCD_OSC_OV_INTER_vect  +=[FCD_OSC_OV_INTER]

                        # FCD_list.append(FCD)
                except:
                    print('FCD ERROR')
                    pass

                ## FC,SC 

                print('FC - FS Biomarker')

                try:
                    bold_emp       = jul.load_subject_fc_100(mysubj)

                    rsFC_emp       = np.corrcoef(bold_emp.T)
                    rsFC_emp       = rsFC_emp * (rsFC_emp>0)
                    rsFC_emp       = rsFC_emp - np.diag(np.diagonal(rsFC_emp))

                    rsFC           = np.corrcoef(bold_sig.T)
                    rsFC           = rsFC * (rsFC>0)
                    rsFC           = rsFC - np.diag(np.diagonal(rsFC))

                    print(rsFC)

                    FC_CORR        = np.corrcoef([rsFC.ravel(),rsFC_emp.ravel()])[0,1]
                    FS_CORR        = np.corrcoef([rsFC.ravel(),weights_symm.ravel()])[0,1]

                    print(FC_CORR)

                except:

                    rsFC_emp = np.zeros([weights.shape[0],weights.shape[1]])
                    rsFC     = np.zeros([weights.shape[0],weights.shape[1]])
                    FC_CORR  = 0
                    FS_CORR  = 0

                    print('FC - FS ERROR')
                    pass

                np.savez(out_path, mysubj = mysubj, myage = myage, Bold_data = Bold_data, Bold_time = Bold_time,
                FCD_SUM_vect = FCD_SUM_vect, FCD_MEAN_vect = FCD_MEAN_vect,
                FCD_VAR_vect = FCD_VAR_vect, FCD_OSC_vect = FCD_OSC_vect,
                FCD_SUM_INTER_vect = FCD_SUM_INTER_vect, FCD_MEAN_INTER_vect = FCD_MEAN_INTER_vect,
                FCD_VAR_INTER_vect = FCD_VAR_INTER_vect, FCD_OSC_INTER_vect = FCD_OSC_INTER_vect,
                FCD_SUM_OV_vect = FCD_SUM_OV_vect, FCD_MEAN_OV_vect = FCD_MEAN_OV_vect,
                FCD_VAR_OV_vect = FCD_VAR_OV_vect, FCD_OSC_OV_vect = FCD_OSC_OV_vect,
                FCD_SUM_OV_INTER_vect = FCD_SUM_OV_INTER_vect, FCD_MEAN_OV_INTER_vect = FCD_MEAN_OV_INTER_vect,
                FCD_VAR_OV_INTER_vect = FCD_VAR_OV_INTER_vect, FCD_OSC_OV_INTER_vect = FCD_OSC_OV_INTER_vect,
                rsFC    = rsFC, rsFC_emp   = rsFC_emp, SC_emp = weights,
                FC_CORR = FC_CORR,FS_CORR = FS_CORR)


            except:

                Bold_data      = []
                Bold_time      = []
                rsFC_emp       = []
                rsFC           = []
                FC_CORR        = []
                FS_CORR        = []

                FCD_SUM_vect   = []
                FCD_MEAN_vect  = []
                FCD_VAR_vect   = []
                FCD_OSC_vect   = []

                FCD_SUM_INTER_vect  = []
                FCD_MEAN_INTER_vect = []
                FCD_VAR_INTER_vect  = []
                FCD_OSC_INTER_vect  = []

                FCD_SUM_OV_vect  = []
                FCD_MEAN_OV_vect = []
                FCD_VAR_OV_vect  = []
                FCD_OSC_OV_vect  = []

                FCD_SUM_OV_INTER_vect  = []
                FCD_MEAN_OV_INTER_vect = []
                FCD_VAR_OV_INTER_vect  = []
                FCD_OSC_OV_INTER_vect  = []


                weights        = []
                TemporalAverage_time = []
                TemporalAverage_data = []

                #SAVE DT
                np.savez(out_path, mysubj = mysubj, myage = myage, Bold_data = Bold_data, Bold_time = Bold_time,
                FCD_SUM_vect = FCD_SUM_vect, FCD_MEAN_vect = FCD_MEAN_vect,
                FCD_VAR_vect = FCD_VAR_vect, FCD_OSC_vect = FCD_OSC_vect,
                FCD_SUM_INTER_vect = FCD_SUM_INTER_vect, FCD_MEAN_INTER_vect = FCD_MEAN_INTER_vect,
                FCD_VAR_INTER_vect = FCD_VAR_INTER_vect, FCD_OSC_INTER_vect = FCD_OSC_INTER_vect,
                FCD_SUM_OV_vect = FCD_SUM_OV_vect, FCD_MEAN_OV_vect = FCD_MEAN_OV_vect,
                FCD_VAR_OV_vect = FCD_VAR_OV_vect, FCD_OSC_OV_vect = FCD_OSC_OV_vect,
                FCD_SUM_OV_INTER_vect = FCD_SUM_OV_INTER_vect, FCD_MEAN_OV_INTER_vect = FCD_MEAN_OV_INTER_vect,
                FCD_VAR_OV_INTER_vect = FCD_VAR_OV_INTER_vect, FCD_OSC_OV_INTER_vect = FCD_OSC_OV_INTER_vect,
                rsFC    = rsFC, rsFC_emp   = rsFC_emp, SC_emp = weights,
                FC_CORR = FC_CORR,FS_CORR = FS_CORR)

                print('CODE ERROR')
                pass
        else:

            print('Simulation data are present',flush=True)
            print(Bold_input_data)

            mysubj          = input_data['mysubj']
            myage           = input_data['myage']
            Bold_data       = input_data['Bold_data']
            Bold_time       = input_data['Bold_time']
            weights         = input_data['SC_emp']

            FCD_SUM_vect    = input_data['FCD_SUM_vect']
            FCD_MEAN_vect   = input_data['FCD_MEAN_vect']
            FCD_VAR_vect    = input_data['FCD_VAR_vect']
            FCD_OSC_vect    = input_data['FCD_OSC_vect']

            try:

                FCD_SUM_OV_vect    = input_data['FCD_SUM_OV_vect']
                FCD_MEAN_OV_vect   = input_data['FCD_MEAN_OV_vect']
                FCD_VAR_OV_vect    = input_data['FCD_VAR_OV_vect']
                FCD_OSC_OV_vect    = input_data['FCD_OSC_OV_vect']

                FCD_SUM_INTER_vect    = input_data['FCD_SUM_INTER_vect']
                FCD_MEAN_INTER_vect   = input_data['FCD_MEAN_INTER_vect']
                FCD_VAR_INTER_vect    = input_data['FCD_VAR_INTER_vect']
                FCD_OSC_INTER_vect    = input_data['FCD_OSC_INTER_vect']

                FCD_SUM_OV_INTER_vect    = input_data['FCD_SUM_OV_INTER_vect']
                FCD_MEAN_OV_INTER_vect   = input_data['FCD_MEAN_OV_INTER_vect']
                FCD_VAR_OV_INTER_vect    = input_data['FCD_VAR_OV_INTER_vect']
                FCD_OSC_OV_INTER_vect    = input_data['FCD_OSC_OV_INTER_vect']
            
            except:

                print('Re-compute INTER - OV FCD Biomarker')
                NHALF                                   = int(weights.shape[0]/2)
                mask_inter                              = np.zeros(weights.shape)
                mask_inter[0:NHALF,NHALF:NHALF*2]       = 1
                mask_inter[NHALF:NHALF*2,0:NHALF]       = 1

                bold_sig       =  Bold_data[:,0,:,0]

                win_FCD_vect   = [10e3,20e3,30e3,40e3,50e3]
                # win_FCD_vect  = [1e3,2e3]
                # FCD_list      = []

                FCD_SUM_INTER_vect  = []
                FCD_MEAN_INTER_vect = []
                FCD_VAR_INTER_vect  = []
                FCD_OSC_INTER_vect  = []

                FCD_SUM_OV_vect  = []
                FCD_MEAN_OV_vect = []
                FCD_VAR_OV_vect  = []
                FCD_OSC_OV_vect  = []

                FCD_SUM_OV_INTER_vect  = []
                FCD_MEAN_OV_INTER_vect = []
                FCD_VAR_OV_INTER_vect  = []
                FCD_OSC_OV_INTER_vect  = []

                try:
                    
                    for win_idx in range(len(win_FCD_vect)):
                        win_FCD      = win_FCD_vect[win_idx]
                        transient    = int(5e3/2000)
                        FCD, fc_stack, speed_fcd                = analysis.compute_fcd(bold_sig[transient:,:], win_len=int(win_FCD/2000), win_sp=1)
                        fcd_inter, fc_stack_inter, _            = analysis.compute_fcd_filt(bold_sig[transient:,:],mask_inter,win_len=int(win_FCD/2000),win_sp=1)

                        FCD_OV, fc_stack_ov, _                  = analysis.compute_fcd(bold_sig[transient:,:], win_len=int(win_FCD/2000), win_sp=int(win_FCD/2000))
                        FCD_OV_INTER, fc_stack_ov_inter, _      = analysis.compute_fcd_filt(bold_sig[transient:,:],mask_inter,win_len=int(win_FCD/2000),win_sp=int(win_FCD/2000))
                        print('FCD')

                        FCD_TRIU_OV          = np.triu(FCD, k=int(win_FCD/2000))

                        FCD_INTER_TRIU       = np.triu(fcd_inter, k=1)
                        FCD_INTER_TRIU_OV    = np.triu(fcd_inter, k=int(win_FCD/2000))

                        print('FCD_TRIU')

                        FCD_OSC_INTER        = np.std(fc_stack_inter)

                        FCD_SUM_OV           = sum(sum(FCD_TRIU_OV))
                        FCD_MEAN_OV          = np.mean(FCD_TRIU_OV)
                        FCD_VAR_OV           = np.var(FCD_TRIU_OV)
                        FCD_OSC_OV           = np.std(fc_stack_ov)
                        FCD_OSC_OV_INTER     = np.std(fc_stack_ov_inter)

                        FCD_SUM_INTER        = sum(sum(FCD_INTER_TRIU))
                        FCD_MEAN_INTER       = np.mean(FCD_INTER_TRIU)
                        FCD_VAR_INTER        = np.var(FCD_INTER_TRIU)

                        FCD_SUM_OV_INTER     = sum(sum(FCD_INTER_TRIU_OV))
                        FCD_MEAN_OV_INTER    = np.mean(FCD_INTER_TRIU_OV)
                        FCD_VAR_OV_INTER     = np.var(FCD_INTER_TRIU_OV)

                        FCD_SUM_INTER_vect  +=[FCD_SUM_INTER]
                        FCD_MEAN_INTER_vect +=[FCD_MEAN_INTER]
                        FCD_VAR_INTER_vect  +=[FCD_VAR_INTER]
                        FCD_OSC_INTER_vect  +=[FCD_OSC_INTER]

                        FCD_SUM_OV_vect  +=[FCD_SUM_OV]
                        FCD_MEAN_OV_vect +=[FCD_MEAN_OV]
                        FCD_VAR_OV_vect  +=[FCD_VAR_OV]
                        FCD_OSC_OV_vect  +=[FCD_OSC_OV]

                        FCD_SUM_OV_INTER_vect  +=[FCD_SUM_OV_INTER]
                        FCD_MEAN_OV_INTER_vect +=[FCD_MEAN_OV_INTER]
                        FCD_VAR_OV_INTER_vect  +=[FCD_VAR_OV_INTER]
                        FCD_OSC_OV_INTER_vect  +=[FCD_OSC_OV_INTER]

                except:
                    print('FCD ERROR')
                    pass

                pass

            rsFC                  = input_data['rsFC']
            rsFC_emp              = input_data['rsFC_emp']
            FC_CORR               = input_data['FC_CORR']
            FS_CORR               = input_data['FS_CORR']

            #SAVE DT
            np.savez(out_path, mysubj = mysubj, myage = myage, Bold_data = Bold_data, Bold_time = Bold_time,
            FCD_SUM_vect = FCD_SUM_vect, FCD_MEAN_vect = FCD_MEAN_vect,
            FCD_VAR_vect = FCD_VAR_vect, FCD_OSC_vect = FCD_OSC_vect,
            FCD_SUM_INTER_vect = FCD_SUM_INTER_vect, FCD_MEAN_INTER_vect = FCD_MEAN_INTER_vect,
            FCD_VAR_INTER_vect = FCD_VAR_INTER_vect, FCD_OSC_INTER_vect = FCD_OSC_INTER_vect,
            FCD_SUM_OV_vect = FCD_SUM_OV_vect, FCD_MEAN_OV_vect = FCD_MEAN_OV_vect,
            FCD_VAR_OV_vect = FCD_VAR_OV_vect, FCD_OSC_OV_vect = FCD_OSC_OV_vect,
            FCD_SUM_OV_INTER_vect = FCD_SUM_OV_INTER_vect, FCD_MEAN_OV_INTER_vect = FCD_MEAN_OV_INTER_vect,
            FCD_VAR_OV_INTER_vect = FCD_VAR_OV_INTER_vect, FCD_OSC_OV_INTER_vect = FCD_OSC_OV_INTER_vect,
            rsFC    = rsFC, rsFC_emp   = rsFC_emp, SC_emp = weights,
            FC_CORR = FC_CORR,FS_CORR = FS_CORR)

            print('Transfer to data output complete',flush=True)
        
    except:

        print('THIS WORKING WAS NOT COMPUTED',flush=True)

        Bold_data      = []
        Bold_time      = []
        rsFC_emp       = []
        rsFC           = []
        FC_CORR        = []
        FS_CORR        = []

        FCD_SUM_vect   = []
        FCD_MEAN_vect  = []
        FCD_VAR_vect   = []
        FCD_OSC_vect   = []

        FCD_SUM_INTER_vect  = []
        FCD_MEAN_INTER_vect = []
        FCD_VAR_INTER_vect  = []
        FCD_OSC_INTER_vect  = []

        FCD_SUM_OV_vect  = []
        FCD_MEAN_OV_vect = []
        FCD_VAR_OV_vect  = []
        FCD_OSC_OV_vect  = []

        FCD_SUM_OV_INTER_vect  = []
        FCD_MEAN_OV_INTER_vect = []
        FCD_VAR_OV_INTER_vect  = []
        FCD_OSC_OV_INTER_vect  = []


        weights        = []
        TemporalAverage_time = []
        TemporalAverage_data = []

        jul                                       = data.Julich()
        subjs                                     = jul.list_subjects()
        subj_age,gender,education,subj_ID,_,_,_,_ = jul.metadata()

        SUBJ_TARG     = [subj_loc for subj_loc in range(len(subj_ID)) if mysubj in subj_ID[subj_loc] ][0]
        myage         = subj_age[SUBJ_TARG]

        #SAVE DT
        np.savez(out_path, mysubj = mysubj, myage = myage, Bold_data = Bold_data, Bold_time = Bold_time,
        FCD_SUM_vect = FCD_SUM_vect, FCD_MEAN_vect = FCD_MEAN_vect,
        FCD_VAR_vect = FCD_VAR_vect, FCD_OSC_vect = FCD_OSC_vect,
        FCD_SUM_INTER_vect = FCD_SUM_INTER_vect, FCD_MEAN_INTER_vect = FCD_MEAN_INTER_vect,
        FCD_VAR_INTER_vect = FCD_VAR_INTER_vect, FCD_OSC_INTER_vect = FCD_OSC_INTER_vect,
        FCD_SUM_OV_vect = FCD_SUM_OV_vect, FCD_MEAN_OV_vect = FCD_MEAN_OV_vect,
        FCD_VAR_OV_vect = FCD_VAR_OV_vect, FCD_OSC_OV_vect = FCD_OSC_OV_vect,
        FCD_SUM_OV_INTER_vect = FCD_SUM_OV_INTER_vect, FCD_MEAN_OV_INTER_vect = FCD_MEAN_OV_INTER_vect,
        FCD_VAR_OV_INTER_vect = FCD_VAR_OV_INTER_vect, FCD_OSC_OV_INTER_vect = FCD_OSC_OV_INTER_vect,
        rsFC    = rsFC, rsFC_emp   = rsFC_emp, SC_emp = weights,
        FC_CORR = FC_CORR,FS_CORR = FS_CORR)

        print('ABSENT WORKING POINT',flush=True)
        pass

    CPU_TIME    = time.time() - t0
    print(f'CPU TIME {CPU_TIME}')
