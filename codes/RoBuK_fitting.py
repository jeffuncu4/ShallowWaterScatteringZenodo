import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from calculate_fluxratio import flux_ratio
from pathlib import Path
from AnalysisClass import Analysis
from scipy.optimize import curve_fit
from calculate_fluxratio import flux_ratio
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple
from matplotlib.legend import Legend



plt.rcParams['axes.titlesize'] = 'x-large'
plt.rcParams['axes.labelsize'] = 'x-large'
plt.rcParams['legend.fontsize'] = 'x-large'


analytics_dir = Path('/media/jeff/Data/SWportableData/analytics/')
# pandas_filename = analytics_dir/ 'piR_adjusted_params_ens.csv'
pandas_filename = analytics_dir/ 'piR_adjusted_params_ens_local.csv'



ro_cyc_list = [0.01, 0.03, 0.04, 0.06]
ro_anti_list = [0.006, 0.01, 0.02, 0.03]
bu_list = [0.5, 0.9, 1.0, 1.1, 1.5]
lr_list= [1.,1.5,  2., 3., 4.]


flux_max = 0.8 # no stricly neccesary but can be used to filter through experiments to do subfits
df = pd.read_csv(pandas_filename)
df = df.query('flux_ratio < @flux_max')
# df = df.query('Ro> 0.03 or Ro< 0.03')
# df = df.query('Ro> 0.021 or Ro< 0.021')
# df = df.query('Bu>0.5')
# df = df.query('Lr > 1')
# cyclonic=True
# print('cyclonic = ', cyclonic)
# df = df.query('cyclonic==@cyclonic')
Lr = df['Lr'].values
Bu = df['Bu'].values
Ro_bulk_adj = df['Ro_bulk_adj'].values
Bu_adj = df['Bu_adj'].values
K_adj = df['K_adj'].values
ens_adj = df['ens_adj'].values
flux_ratio = df['flux_ratio'].values
cyc_bool = df['cyclonic'].values
print (np.shape(ens_adj), type(ens_adj))

# just for local ro


local_ro = df['local_ro'].values

#L-adj  =K*adj*la_orig = K_adj*L/Lr
L = 25e3
L_adj = K_adj*L/Lr
# print (L_adj/L)
ens_orig = ens_adj*L_adj**2/L**2


# ens_adj = 1*ens_orig

Ro_bulk_orig = Ro_bulk_adj*L_adj/L

data_ens = np.stack([ens_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
data_ro = np.stack([Ro_bulk_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
data_zeta = np.stack([local_ro, Bu_adj, K_adj, flux_ratio], axis =-1)


print (np.shape(data_ro))


def fit_func(data, A, alpha, beta, gamma):
    ros = data[:,0]
    bus = data[:,1]
    ks = data[:,2]
    return 2/np.pi*np.arctan(A*ros**alpha*bus**beta*ks**gamma)
guess = (0.1, 1, -1, 2)

# THREE DIMENSIONAL FITTING
ens_p, ens_pcov = curve_fit(fit_func, data_ens[:,:3], data_ens[:,3], guess)
ro_p, ro_pcov = curve_fit(fit_func, data_ro[:,:3], data_ro[:,3], guess)
zeta_p, zeta_pcov = curve_fit(fit_func, data_zeta[:,:3], data_ro[:,3], guess)

#PRINT PARAMETERS
print ('enstrophy pararms',  '\n')
print ('A, alpha, beta, gamma')
for i in range(len(guess)):
    print (ens_p[i], ens_pcov[i,i]**0.5)
print (  '\n')
print ('ro_params', '\n')
for i in range(len(guess)):
    print (ro_p[i], ro_pcov[i,i]**0.5)
print (  '\n')
print ('zeta_params', '\n')
for i in range(len(guess)):
    print (zeta_p[i], zeta_pcov[i,i]**0.5)

#MAKE FITTED CURVES
ens_fit = fit_func(data_ens, *ens_p)
ro_fit = fit_func(data_ro, *ro_p)
zeta_fit = fit_func(data_zeta, *zeta_p)


logens = np.log(ens_fit )
logro = np.log(ro_fit )
logzeta = np.log(zeta_fit)
logflux = np.log(flux_ratio)

logens_p = np.polyfit(logens, logflux , 1)
logro_p = np.polyfit(logro, logflux , 1)
logzeta_p = np.polyfit(logzeta, logflux , 1)

log_ens_func = np.poly1d(logens_p)

# plt.plot( logens,log_ens_func(logens) )
# plt.show()

markersize = 3
# plt.scatter(ens_p[0]*ens_adj**ens_p[1]*Bu_adj**ens_p[2]*K_adj**ens_p[3], flux_ratio, c = 'r', label = 'ens')
# plt.scatter(ro_p[0]*Ro_bulk_adj**ro_p[1]*Bu_adj**ro_p[2]*K_adj**ro_p[3], flux_ratio, c = 'b',  label = 'ro')
# plt.legend()
# plt.show()

# plt.loglog(ens_p[0]*ens_adj**ens_p[1]*Bu_adj**ens_p[2]*K_adj**ens_p[3], flux_ratio, c = 'r', label = 'ens',\
#      marker='d', lw=0, markersize =markersize)
# plt.loglog(ro_p[0]*Ro_bulk_adj**ro_p[1]*Bu_adj**ro_p[2]*K_adj**ro_p[3], flux_ratio, c = 'b',  label = 'ro',\
#      marker='d', markersize =markersize, lw=0)
# plt.legend()
# plt.show()


Lr_dic = {1:'r',1.5:'g', 2:'b', 3:'k', 4:'y'}

# plt.scatter(ens_fit, flux_ratio, c = 'r', label = 'ens')
# plt.scatter(ro_fit, flux_ratio, c = 'b',  label = 'ro')
# plt.scatter(zeta_fit, flux_ratio, c = 'k',  label = 'zeta')


#LEGEND COLORS
Lr_dic = {1:'r',1.5:'g', 2:'b', 3:'k', 4:'y'}
bu_marker_dic  = {0.5:'s',0.9:'d', 1.0:'h', 1.1:'v', 1.5:'^'}
###better plotting
marker_factor = 20
# for i in range(len(df)):
#     rr = data_ens[i, 0]

#     bb = data_ens[i, 1]
#     ll = data_ens[i, 2]
#     fr = data_ens[i, 3]
#     lr = Lr[i]
#     cyc = cyc_bool[i]
#     x_axis = ens_p[0]*rr**ens_p[1]*bb**ens_p[2]#*ll**ens_p[3]
#     # x_axis = ens_p[0]*np.arctan(ens_p[1]*rr**ens_p[2]*bb**ens_p[3])#*ll**ens_p[3]
#     if cyc:
#         plt.plot(x_axis,fr, c =Lr_dic[lr], marker='d', \
#             markersize = rr*marker_factor+2 )
#     elif not cyc:
#         plt.plot(-x_axis,fr, c =Lr_dic[lr], marker='d', \
#             markerfacecolor='none',  markersize = rr*marker_factor+2 )
# plt.show()


fig, ax = plt.subplots( figsize=(8, 6))

marker_factor = 30
for index, row in df.iterrows():
    Lr = row['Lr']
    Bu = row['Bu']
    Ro_bulk_adj = np.array(row['Ro_bulk_adj'])
    Bu_adj = np.array(row['Bu_adj'])
    K_adj = np.array(row['K_adj'])
    ens_adj = np.array(row['ens_adj'])
    flux_ratio = np.array(row['flux_ratio'])
    
    cyc_bool = row['cyclonic']
    local_ro = row['local_ro']

    data_ens = np.stack([ens_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_ens = np.expand_dims(data_ens, axis = 0)
    ens_fit = fit_func(data_ens, *ens_p)

    data_zeta = np.stack([local_ro, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_zeta = np.expand_dims(data_zeta, axis = 0)
    zeta_fit = fit_func(data_zeta, *zeta_p)
    
    # data_ro = np.stack([Ro_bulk_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    # data_ro = np.expand_dims(data_zeta, axis = 0)
    # ro_fit = fit_func(data_zeta, *zeta_p)

    # rounded_fit = fit_func(data_ens, 0.651, 1, -1, 2, 0.5)
    # rounded_fit = Ro_bulk_adj**2.079/Bu_adj**0.975*K_adj**2*0.65*9.2
    rounded_fit = Ro_bulk_adj**2.0/Bu_adj**1*K_adj**2*5

    if cyc_bool:
        ax.plot(rounded_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
            markersize = Ro_bulk_adj*marker_factor+2 )
    elif not cyc_bool:
        ax.plot(-rounded_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
            markerfacecolor='none',  markersize = Ro_bulk_adj*marker_factor+2 )

    # plt.loglog(rounded_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #         markersize = Ro_bulk_adj*marker_factor + 2 )
    # plt.loglog(10*ens_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #         markersize = Ro_bulk_adj*marker_factor+2 )
  
print (index)

y_ext = flux_max
ax.set_title('Rounded Scaling')
# plt.loglog([-y_ext,0, y_ext], [y_ext,0, y_ext], label ='rounded_fit')
# plt.loglog([-y_ext*10,0, y_ext*10], [y_ext,0, y_ext], label ='optimized_fit')
ax.plot([-y_ext,0, y_ext], [y_ext,0, y_ext], 'k--', label ='rounded_fit')
# plt.plot([-y_ext,0, y_ext], [y_ext*1.5**2,0, y_ext*1.5**2], label ='rounded_fit')
# plt.plot([-y_ext,0, y_ext], [y_ext*2**2,0, y_ext*2**2], label ='rounded_fit')
# plt.plot([-y_ext,0, y_ext], [y_ext*4**2,0, y_ext*4**2], label ='rounded_fit')
# plt.plot([-y_ext,0, y_ext], [y_ext*3**2,0, y_ext*3**2], label ='rounded_fit')

ax.set_xlabel(r'$5Fr^2K^2$', fontsize = 'x-large')
ax.set_ylabel(r'$S$', fontsize = 'x-large')
ax.grid(True)
markersize=10
legend_elements = [Line2D([0], [0], color='r', marker='o', markersize =markersize, lw=0, label=r'$K$ = 1'),
                   Line2D([0], [0], color='g', marker='o', lw=0,markersize =markersize, label=r'$K$ = 1.5'),
                   Line2D([0], [0], color='b', marker='o',lw=0,markersize =markersize ,label=r'$K$ = 2'),
                   Line2D([0], [0], color='k', marker='o', lw=0,markersize =markersize, label=r'$K$ = 3'),
                   Line2D([0], [0], color='y', marker='o', lw=0,markersize =markersize, label=r'$K$ = 4')]

legend_elements2 = [Line2D([0], [0], color='k', marker='s', lw=0,markersize =markersize, label=r'$Bu$ = 0.5'),
                   Line2D([0], [0], color='k', marker='d', lw=0,markersize =markersize, label=r'$Bu$ = 0.9'),
                   Line2D([0], [0], color='k', marker='h',lw=0,markersize =markersize, label=r'$Bu$ = 1.0'),
                   Line2D([0], [0], color='k', marker='v', lw=0,markersize =markersize, label=r'$Bu$ = 1.1'),
                   Line2D([0], [0], color='k', marker='^', lw=0,markersize =markersize, label=r'$Bu$ = 1.5')]

first_legend = plt.legend(handles=legend_elements, loc='lower right', fontsize = 'medium')

plt.gca().add_artist(first_legend)
plt.legend(handles=legend_elements2, loc='center right', fontsize = 'medium')



bbox_to_anchor = (0.0, 0.6)
inset_ax = inset_axes(ax,
                        width='30%',  # Adjust the inset width as needed
                        height='30%',  # Adjust the inset height as needed
                        loc='upper center',
                        bbox_to_anchor=(0.-0.15,0.0,1,1),
                        bbox_transform=ax.transAxes)
    # Plot the zoomed-in data in the inset

for index, row in df.iterrows():
    Lr = row['Lr']
    Bu = row['Bu']
    Ro_bulk_adj = np.array(row['Ro_bulk_adj'])
    Bu_adj = np.array(row['Bu_adj'])
    K_adj = np.array(row['K_adj'])
    ens_adj = np.array(row['ens_adj'])
    flux_ratio = np.array(row['flux_ratio'])
    
    cyc_bool = row['cyclonic']
    local_ro = row['local_ro']

    data_ens = np.stack([ens_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_ens = np.expand_dims(data_ens, axis = 0)
    ens_fit = fit_func(data_ens, *ens_p)

    data_zeta = np.stack([local_ro, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_zeta = np.expand_dims(data_zeta, axis = 0)
    zeta_fit = fit_func(data_zeta, *zeta_p)
    
    # data_ro = np.stack([Ro_bulk_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    # data_ro = np.expand_dims(data_zeta, axis = 0)
    # ro_fit = fit_func(data_zeta, *zeta_p)

    # rounded_fit = fit_func(data_ens, 0.651, 1, -1, 2, 0.5)
    # rounded_fit = Ro_bulk_adj**2.079/Bu_adj**0.975*K_adj**2*0.65*9.2
    rounded_fit = Ro_bulk_adj**2.0/Bu_adj**1*K_adj**2*5

    if cyc_bool:
        inset_ax.plot(rounded_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
            markersize = Ro_bulk_adj*marker_factor+2 )
    elif not cyc_bool:
        inset_ax.plot(-rounded_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
            markerfacecolor='none',  markersize = Ro_bulk_adj*marker_factor+2 )
    inset_ax.set_xlim(-0.10, 0.10)  # Adjust the x-axis limits for zooming in
    inset_ax.set_ylim(0, 0.1)    # Adjust the y-axis limits for zooming in
    inset_ax.grid(True)
    inset_ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) 
    inset_ax.plot([-y_ext,0, y_ext], [y_ext,0, y_ext], 'k--', label ='rounded_fit')


plt.show()


marker_factor = 30
def plot_fit(flag, ax):

    
    for index, row in df.iterrows():
        Lr = row['Lr']
        Bu = row['Bu']
        Ro_bulk_adj = np.array(row['Ro_bulk_adj'])
        Bu_adj = np.array(row['Bu_adj'])
        K_adj = np.array(row['K_adj'])
        ens_adj = np.array(row['ens_adj'])
        flux_ratio = np.array(row['flux_ratio'])
        
        cyc_bool = row['cyclonic']
        local_ro = row['local_ro']

        data_ens = np.stack([ens_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
        data_ens = np.expand_dims(data_ens, axis = 0)
        ens_fit = fit_func(data_ens, *ens_p)

        data_zeta = np.stack([local_ro, Bu_adj, K_adj, flux_ratio], axis =-1)
        data_zeta = np.expand_dims(data_zeta, axis = 0)
        zeta_fit = fit_func(data_zeta, *zeta_p)
        
        data_ro = np.stack([Ro_bulk_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
        data_ro = np.expand_dims(data_ro, axis = 0)
        ro_fit = fit_func(data_ro, *ro_p)

        if flag==0:
            data_fit = ens_fit
        elif flag==1:
            data_fit = ro_fit
        elif flag==2:
            data_fit = zeta_fit


        if cyc_bool:
            ax.plot(data_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
                markersize = Ro_bulk_adj*marker_factor+2 )
        elif not cyc_bool:
            ax.plot(-data_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
                markerfacecolor='none',  markersize = Ro_bulk_adj*marker_factor+2 )
        # if cyc_bool:
        #     plt.plot(zeta_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
        #         markersize = Ro_bulk_adj*marker_factor+2 )
        # elif not cyc_bool:
        #     plt.plot(-zeta_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
        #         markerfacecolor='none',  markersize = Ro_bulk_adj*marker_factor+2 )
    print (index)

    y_ext = flux_max

    ax.plot([-y_ext,0, y_ext], [y_ext,0, y_ext], 'k--')



x = [1, 2, 3, 4]
y = [10, 15, 7, 12]

# Create a figure with three subplots
fig, axs = plt.subplots(1, 3, figsize=(11, 5))

axs[0].set_ylabel(r'$S$', fontsize='x-large')

xlabel_titles = [r'$S^{\theta}_{Z}(Z = \varepsilon)$',\
    r'$S^{\theta}_Z(Z=Ro_{b})$',\
        r'$S^{\theta}_Z(Z=Ro_{\zeta})$']

titles = ['(a) Enstrophy', '(b) Bulk Rossby Number', '(c) Local Rossby Number']

# Loop through each subplot
for i, ax in enumerate(axs):
    # Plot the main data in the subplot
    plot_fit(i, ax)
    ax.set_title(titles[i])
    ax.set_xlabel(xlabel_titles[i], fontsize='xx-large')
    ax.grid(True)
    # Create an inset axes in the current subplot
    inset_ax = inset_axes(ax,
                          width='33%',  # Adjust the inset width as needed
                          height='33%',  # Adjust the inset height as needed
                          loc='upper center')  # Position of the inset within the subplot

    # Plot the zoomed-in data in the inset
    plot_fit(i, inset_ax)
    inset_ax.set_xlim(-0.10, 0.10)  # Adjust the x-axis limits for zooming in
    inset_ax.set_ylim(0, 0.1)    # Adjust the y-axis limits for zooming in
    inset_ax.grid(True)
    inset_ax.tick_params(axis='y', which='both', left=False, right=False, labelleft=False) 
    
# Adjust layout and show the plot
plt.tight_layout()

plt.show()

fig, ax = plt.subplots(figsize=(8, 6))

marker_factor = 30
for index, row in df.iterrows():
    Lr = row['Lr']
    Bu = row['Bu']
    Ro_bulk_adj = np.array(row['Ro_bulk_adj'])
    Bu_adj = np.array(row['Bu_adj'])
    K_adj = np.array(row['K_adj'])
    ens_adj = np.array(row['ens_adj'])
    flux_ratio = np.array(row['flux_ratio'])
    
    cyc_bool = row['cyclonic']
    local_ro = row['local_ro']

    data_ens = np.stack([ens_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_ens = np.expand_dims(data_ens, axis = 0)
    ens_fit = fit_func(data_ens, *ens_p)

    data_zeta = np.stack([local_ro, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_zeta = np.expand_dims(data_zeta, axis = 0)
    zeta_fit = fit_func(data_zeta, *zeta_p)

    data_ro = np.stack([Ro_bulk_adj, Bu_adj, K_adj, flux_ratio], axis =-1)
    data_ro = np.expand_dims(data_ro, axis = 0)
    ro_fit = fit_func(data_ro, *ro_p)

    if cyc_bool:
        markerfacecolor = 'none'
    else:
        markerfacecolor=Lr_dic[Lr]
    # plt.loglog(ens_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #     markerfacecolor=markerfacecolor, markersize = Ro_bulk_adj*marker_factor+2 )
    

    # plt.loglog(ro_fit*10,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #     markerfacecolor=markerfacecolor,  markersize = Ro_bulk_adj*marker_factor+2 )
    
    # plt.loglog(zeta_fit*100,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #     markerfacecolor=markerfacecolor,  markersize = Ro_bulk_adj*marker_factor+2 )
    ax.loglog(ens_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
        markerfacecolor=markerfacecolor, markersize = 3 )
    

    ax.loglog(ro_fit*10,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
        markerfacecolor=markerfacecolor,  markersize = 3 )
    
    ax.loglog(zeta_fit*100,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
        markerfacecolor=markerfacecolor,  markersize = 3)



    # if cyc_bool:
    #     plt.plot(zeta_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #         markersize = Ro_bulk_adj*marker_factor+2 )
    # elif not cyc_bool:
    #     plt.plot(-zeta_fit,flux_ratio, c =Lr_dic[Lr], marker=bu_marker_dic[Bu], \
    #         markerfacecolor='none',  markersize = Ro_bulk_adj*marker_factor+2 )

ax.set_title('Log of Fits')
ax.loglog([0, 1], [0,1], 'k--', label = r'$\varepsilon$')
ax.loglog([0, 10], [0,1],'k--', label = r'$Ro_b$')
ax.loglog([0, 100], [0,1], 'k--', label = r'$Ro_{\zeta}$')
markersize=10


ax.set_ylabel(r'$S$')
ax.set_xlabel(r'$S^{\theta}_Z$')

legend_elements = [Line2D([0], [0], color='r', marker='o', markersize =markersize, lw=0, label=r'$K$ = 1'),
                   Line2D([0], [0], color='g', marker='o', lw=0,markersize =markersize, label=r'$K$ = 1.5'),
                   Line2D([0], [0], color='b', marker='o',lw=0,markersize =markersize ,label=r'$K$ = 2'),
                   Line2D([0], [0], color='k', marker='o', lw=0,markersize =markersize, label=r'$K$ = 3'),
                   Line2D([0], [0], color='y', marker='o', lw=0,markersize =markersize, label=r'$K$ = 4')]

legend_elements2 = [Line2D([0], [0], color='k', marker='s', lw=0,markersize =markersize, label=r'$Bu$ = 0.5'),
                   Line2D([0], [0], color='k', marker='d', lw=0,markersize =markersize, label=r'$Bu$ = 0.9'),
                   Line2D([0], [0], color='k', marker='h',lw=0,markersize =markersize, label=r'$Bu$ = 1.0'),
                   Line2D([0], [0], color='k', marker='v', lw=0,markersize =markersize, label=r'$Bu$ = 1.1'),
                   Line2D([0], [0], color='k', marker='^', lw=0,markersize =markersize, label=r'$Bu$ = 1.5')]

first_legend = plt.legend(handles=legend_elements, loc='lower right', fontsize = 'medium')

plt.gca().add_artist(first_legend)
plt.legend(handles=legend_elements2, loc='upper left', fontsize = 'medium')
angle = 60
fontsize_text = 18
plt.text(
    1e-3,             # x-coordinate
    1e-3*2,             # y-coordinate
    r'$\varepsilon$',  # text string
    fontsize=fontsize_text,   # font size
    color='black',   # text color
    rotation=angle)  # rotation angle

plt.text(
    1e-2,             # x-coordinate
    1e-3*2,             # y-coordinate
    r'$Ro_b$',  # text string
    fontsize=fontsize_text,   # font size
    color='black',   # text color
    rotation=angle)  # rotation angle

plt.text(
    1e-1,             # x-coordinate
    1e-3*2,             # y-coordinate
    r'$Ro_{\zeta}$',  # text string
    fontsize=fontsize_text,   # font size
    color='black',   # text color
    rotation=angle)  # rotation angle



plt.grid()
plt.show()
