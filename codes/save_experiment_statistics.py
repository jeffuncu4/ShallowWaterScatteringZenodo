'''
This file grabs the experiments which have been conducted and calculates the adjusted parameters 
and flux ratio and saves it to a csv file using the pandas library.
'''

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calculate_fluxratio import flux_ratio
from pathlib import Path
from AnalysisClass import Analysis
from scipy.optimize import curve_fit


def enstrophy(ana):
    return np.sum(ana.vorticity()**2)*ana.dx**2

def adjusted_params(Ro, Bu, Lr, cyclonic):
    ana = Analysis(Ro, Bu, Lr, 1000.0, cyclonic=cyclonic)
    print (ana.exp_name)
    exp_name  = ana.exp_name
    L = ana.L
    x = ana.x_axis
    f = ana.f
    u_cross = ana.u[0, 256, :]
    u_max = np.max(u_cross)
    umi = np.argmax(u_cross)
    R = np.abs(x[umi])
    print (f'scaled ratio {R/(L/np.pi)}')
    L_adj = np.pi*R
    lam_orig = L/Lr
    K_adj = L_adj/lam_orig
    # Bu_adj = Bu/L**2*L_adj**2
    Bu_adj = Bu*L**2/L_adj**2
    # Bu_adj = Bu*1
    Ro_bulk_adj = u_max/L_adj/f
    flux_rat = flux_ratio(Ro, Bu, Lr, cyclonic)
    ens_adj = enstrophy(ana)/L_adj**2/f**2
    print (f'k_adj: {K_adj}, Bu_adj: {Bu_adj}, Ro: {Ro}')
    local_Ro = np.max(np.abs(ana.vorticity()))/ana.f
    return exp_name, ens_adj, Ro_bulk_adj, local_Ro, Bu_adj, K_adj, flux_rat





def add_to_file(file_name, exp_name, Ro, Bu, Lr, cyclonic,
                ens_adj, Ro_bulk_adj, Bu_adj, K_adj, flux_rat, local_Ro):

    data = {'exp_name': [exp_name], 'Ro':Ro, 'Bu':Bu, 'Lr':Lr,\
            'cyclonic':cyclonic, 'ens_adj':ens_adj, 'Ro_bulk_adj':Ro_bulk_adj, 'Bu_adj':Bu_adj, \
            'K_adj':K_adj, 'flux_ratio': flux_rat, 'local_ro': local_Ro}
    
    
    data = pd.DataFrame(data)
    try:
        df = pd.read_csv(file_name)
        if exp_name not in df['exp_name'].values:
            df = pd.concat([data, df], ignore_index=True)
            df.to_csv(file_name, index=False)
#        print (df)
    except:
        data.to_csv(file_name, index=False)

def run_statistics(file_name, Ro, Bu, Lr, cyclonic):
    exp_name, ens_adj, Ro_bulk_adj, local_Ro, Bu_adj, K_adj, flux_rat = adjusted_params(Ro, Bu, Lr, cyclonic)

    add_to_file(file_name, exp_name, Ro, Bu, Lr, cyclonic,
                ens_adj, Ro_bulk_adj, Bu_adj, K_adj, flux_rat, local_Ro)





analytics_dir = Path('/media/jeff/Data/SWportableData/analytics/')
# name for file to save experiment statistics
pandas_filename = analytics_dir/ 'piR_adjusted_params_ens_local.csv'

ro_cyc_list = [0.01, 0.03, 0.04, 0.06]
ro_anti_list = [0.006, 0.01, 0.02, 0.03]
bu_list = [0.5, 0.9, 1.0, 1.1, 1.5]
lr_list= [1.,1.5,  2., 3., 4.]

def run_multiple(ro_list, bu_list, lr_list, cyclonic):
    for Ro in ro_list:
        print (Ro)
        for Bu in bu_list:
            for Lr in lr_list:
                run_statistics(pandas_filename, Ro, Bu, Lr, cyclonic)

run_multiple(ro_anti_list, bu_list, lr_list, cyclonic = False)
run_multiple(ro_cyc_list, bu_list, lr_list, cyclonic = True)

# def quick_plot(pandas_fn):
#     df = pd.read_csv(pandas_fn)
#     # Ros = df['bulk'].values
#     # Bus = df['Bu'].values
#     # Lrs = df['Lr'].values
#     # flux_rats = df['flux_ratio'].values
#     Ro_bulk_adj = df['Ro_bulk_adj'].values
#     Bu_adj = df['Bu_adj'].values
#     K_adj = df['K_adj'].values
#     flux_rats = df['flux_ratio'].values
    
#     plt.plot(Ro_bulk_ad**1.5/Bu_adj**0.75*K**2)
    # plt.show()
    # plt.plot(ens_adj**1/Bu_adj**1*K**2)
    
