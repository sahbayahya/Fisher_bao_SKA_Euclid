from math import *
from numpy import *
from scipy import *
import cosmolopy.constants as cc
import cosmolopy.distance as cd
import cosmolopy.perturbation as cp
import matplotlib.pyplot as plt
from scipy.integrate import quad
import matplotlib.ticker as tic
from scipy import special
#******************************************************

def D(z):
     h = 0.72
     return cd.angular_diameter_distance(z , **cosmo) * h


def H(z):
	return  cd.hubble_distance_z(z, **cosmo) 
     
     
(z0, Err_lnda0, Err_lnH0,R0,B0) = loadtxt('inputs/output_Fisher_bao_0n_rms_s30k.txt', unpack=True)#;z0 = z0[6:23] ;  Err_lnda0 = Err_lnda0[6:23]
#(z1_new, Err_lnda1_new, Err_lnH1_new,R1_new,B1_new) = loadtxt('inputs/output_Fisher_bao_Euclid_new_2.txt', unpack=True)#;z1 = z1[0:13] ;  Err_lnda1 = Err_lnda1[0:13]
(z1, Err_lnda1, Err_lnH1,R1,B1) = loadtxt('Output_2/output_Fisher_bao_Euclid_s15k_3.txt', unpack=True)#;z1 = z1[0:13] ;  Err_lnda1 = Err_lnda1[0:13]
(z3, Err_lnda3, Err_lnH3, R3,B3) = loadtxt('Output_2/output_Fisher_bao_Opt_3_Srms_corrected_nz_2.txt', unpack=True)#; z3 = z3[6:23] ;  Err_lnda3 = Err_lnda3[6:23]
(z7point3, Err_lnda7point3, Err_lnH7point3,R7point3,B7point3) = loadtxt('Output_2/output_Fisher_bao_Real_5_Srms_corrected_nz_2.txt', unpack=True)#; z7point3 = z7point3[6:20] ;  Err_lnda7point3 = Err_lnda7point3[6:20]
(z23, Err_lnda23, Err_lnH23,R23,B23) = loadtxt('Output_2/output_Fisher_bao_Pess_23_Srms_corrected_nz_2.txt', unpack=True) #;z23 = z23[6:14] ;  Err_lnda23 = Err_lnda23[6:14]
(z70, Err_lnda70, Err_lnH70,R70,B70) = loadtxt('Output_2/output_Fisher_bao_Opt_70_Srms_corrected_nz_2.txt', unpack=True) #; z70 = z70[4:9] ; Err_lnda70 = Err_lnda70[4:9]
(z100, Err_lnda100, Err_lnH100,R100,B100) = loadtxt('Output_2/output_Fisher_bao_Real_150_Srms_corrected_nz_2.txt', unpack=True)##; z100 = z100[4:8] ; Err_lnda100 =Err_lnda100[4:8] 
(z200, Err_lnda200, Err_lnH200,R200,B200) = loadtxt('Output_2/output_Fisher_bao_Pess_200_Srms_corrected_nz_2.txt', unpack=True)#; z200 = z200[4:7] ;  Err_lnda200 =Err_lnda200[4:7] 
#(z200, Err_lnda200_test, Err_lnH200,R200,B200) = loadtxt('output_Fisher_bao_SKA1MID_Pess_Srms_corrected_nz_test.txt', unpack=True)




#==================plot
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
ax.set_yscale('log')
fontsize = 20
# SKA1 
p1, = ax.plot(z70, Err_lnH70, color = '#204a87',linewidth=1.5, linestyle="-", label=r"${\rm SKA1}\ {\rm Opt}$")
p2, = ax.plot(z100, Err_lnH100,  color = '#3465a4',linewidth=1.5, linestyle="-", label=r"${\rm SKA1}$")
p3, = ax.plot(z200, Err_lnH200, color = '#729fcf',linewidth=1.5, linestyle="-", label=r"${\rm SKA1} \ {\rm Pess}$")
#SKA2
p11, = ax.plot(z3, Err_lnH3 , color = '#a40000',linewidth=1.5, linestyle="-", label=r"${\rm SKA2} \ {\rm Opt}$")
p21, = ax.plot(z7point3, Err_lnH7point3 , color = '#cc0000',linewidth=1.5, linestyle="-", label=r"${\rm SKA2}$" )
p31, = ax.plot(z23, Err_lnH23, color = '#ef2929',linewidth=1.5, linestyle="-", label=r"${\rm SKA2} \ {\rm Pess}$")
#Euclid
p4, = ax.plot(z1,  Err_lnH1 , color = '#edd400',linewidth=1.5, linestyle="-", label=r"${\rm Euclid}$")
ax.get_legend_handles_labels()
ax.legend(loc='upper right', ncol=2, frameon=False)
ax.set_title('')
ax.set_xlabel(r"$ { \rm redshift} (z)$", fontsize=fontsize)
ax.set_ylabel(r"$\sigma_{H}/H\%$", fontsize=fontsize)
ax.scatter(z1, Err_lnH1, s= 50, marker='o',   edgecolor = '#edd400', facecolor= '#edd400')
ax.scatter(z3, Err_lnH3,s=50, marker='o', edgecolor = '#a40000', facecolor= '#a40000')
ax.scatter(z7point3, Err_lnH7point3, s= 50, marker='o',  edgecolor = '#cc0000', facecolor='#cc0000')
ax.scatter(z23, Err_lnH23, s=50, marker='o',  edgecolor = '#ef2929',facecolor= '#ef2929')
ax.scatter(z70, Err_lnH70,  s=50, marker='o',  edgecolor = '#204a87',facecolor= '#204a87')
ax.scatter(z100, Err_lnH100,  s=50, marker='o',  edgecolor = '#3465a4',facecolor= '#3465a4')
ax.scatter(z200, Err_lnH200,  s=50, marker='o',  edgecolor = '#729fcf',facecolor= '#729fcf')
ax.get_yaxis().set_major_formatter(tic.ScalarFormatter())
ax.yaxis.set_major_formatter(tic.FormatStrFormatter('%f'))
ax.yaxis.set_major_formatter(tic.FuncFormatter(lambda x, pos: str(int(round(x)))))
#=======legend
l2 = plt.legend([p4,],[r"${\rm Euclid}$"], loc=4, frameon = False)
l1 = plt.legend([p1,p2, p3, p11, p21, p31,],[r"${\rm SKA1}\ {\rm Opt}$", r"${\rm SKA1}$"
    ,  r"${\rm SKA1} \ {\rm Pess}$", r"${\rm SKA2} \ {\rm Opt}$", r"${\rm SKA2}$"
    , r"${\rm SKA2} \ {\rm Pess}$" ] ,loc=1, ncol= 2, frameon= False)
#plt.gca().add_artist(l1)
plt.gca().add_artist(l2)


#======= y axis
yticks = [0.1, 0.2, 0.5, 1 ,2, 5, 10,20,25]
plt.yticks(yticks,[r'$0.1$', r'$0.2$', r'$0.5$',r'$1$',r'$2$',  r'$5$',r'$10$',r'$20$'])
plt.ylim(0.1,  25)
#======= x axis 
plt.xlim(0.,2.1)
xticks = arange(0, 3.0, 0.3)
#plt.xticks(xticks,[r'$0$', r'$0.3$',r'$0.6$', r'$0.9$',r'$1.2$',r'$1.5$',r'$1.8$', r'$2.1$',r'$2.4$',r'$2.7$', r'$3.0$'])
plt.savefig('output_lnH_mario_bias_corrected_nz.pdf')
plt.show()
