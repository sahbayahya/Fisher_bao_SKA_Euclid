'This Program to fit the Bais(z) for different rms sensitvities, produce the fitting values and plot  '
import numpy as np
import scipy.optimize as opt
from scipy import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm

def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    '''
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

# Interface to LineCollection:

def colorline(x, y, z=None, cmap=plt.get_cmap('hot'), norm=plt.Normalize(0.0, 1.0), linewidth=3, alpha=1.0):
    '''
    Plot a colored line with coordinates x and y
    Optionally specify colors in the array z
    Optionally specify a colormap, a norm function and a line width
    '''  
    # Default colors equally spaced on [0,1]:
    if z is None:
        z = np.linspace(0.0, 1.0, len(x))       
    # Special case if a single number:
    if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
        z = np.array([z])        
    z = np.asarray(z)
    segments = make_segments(x, y)
    lc = LineCollection(segments, array=z, cmap=cmap, norm=norm, linewidth=linewidth, alpha=alpha)
    ax = plt.gca()
    ax.add_collection(lc)
    return lc
            
def clear_frame(ax=None): 
    # Taken from a post by Tony S Yu
    if ax is None: 
        ax = plt.gca() 
    ax.xaxis.set_visible(False) 
    ax.yaxis.set_visible(False) 
    for spine in ax.spines.itervalues(): 
        spine.set_visible(False) 


''' '--------------------------------------------------------------------------------------------------------------------------'
This function is the function we think it will be easy
to fit its parameters to our data
'--------------------------------------------------------------------------------------------------------------------------'
'''
def func(p,x):
   w=p[0]*np.exp(p[1]*x)
   print w.size
   return w
''''--------------------------------------------------------------------------------------------------------------------------'
The purpose of this function is finding the 
difference between the theortical and the simulated data point 
at specific point x (or redshift).
'--------------------------------------------------------------------------------------------------------------------------'
'''
def residuals(p,x,y):
   w=p[0]* np.exp(p[1]*x)
   err=w-y
   err=err**2
   B=sum(err)
   return B
''''--------------------------------------------------------------------------------------------------------------------------'
This function call opt.fmin function to  minimize the error on
you paramters and give you the best fit parameters
'--------------------------------------------------------------------------------------------------------------------------'
'''
def find_parameters(p0, x, bias_rms , rms_value, xrange):
	plsqtot = opt.fmin(residuals, p0, args=(x,bias_rms), maxiter=10000, maxfun=10000)
	print  'b(z) = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1]
	ytot=func(plsqtot,xrange)
	return ytot
def read_input(v):
		z = [] ; bz = []
		for i in range(n):
			if Srms_[i] == v:
				z.append(redshift[i]); bz.append(bias[i])
		#print z, bz			
		return z, bz			
if __name__=='__main__':
	
	''' '--------------------------------------------------------------------------------------------------------------------------'
	Input data, note that the bais^2 in the input files
	'''
	(redshift, bias, Srms_) = np.loadtxt('inputs/bias_mill.out', unpack= True)
	n = len(redshift)
	y = ['y0', 'y1', 'y3', 'y5', 'y6', 'y7', 'y10', 'y23', 'y40', 'y70', 'y100', 'y150', 'y200']
	z = ['z0', 'z1', 'z3', 'z5', 'z6', 'z7', 'z10', 'z23', 'z40', 'z70', 'z100', 'z150','z200'] 
	b = ['b0', 'b1', 'b3', 'b5', 'b6', 'b7', 'b10', 'b23', 'b40', 'b70', 'b100', 'b150','b200'] 
	Srms= [0, 1, 3, 5, 6, 7.3, 10, 23,40, 70, 100,200] 
	#y = ['y0', 'y1', 'y3', 'y5', 'y10', 'y23', 'y100']
	#z = ['z0', 'z1', 'z3', 'z5', 'z10', 'z23', 'z100'] 
	#b =['b0', 'b1', 'b3', 'b5', 'b10', 'b23',  'b100'] 
	#Srms= [0, 1, 3, 5, 10, 23, 100, 150] 
	z_ = [] ; b_ = []
	print Srms[0], Srms[1] , Srms[2] , Srms[3] , Srms[4] , Srms[5] , Srms[6]  ,Srms[7]  
	for i in range (len(Srms)):
		z[i], b[i] = read_input(Srms[i])
		z_.append(z[i]) ; b_.append(b[i])
		print  len(z[i]), len(b[i])
	print len(z_),  z_[0], b_[0]
	''''--------------------------------------------------------------------------------------------------------------------------'
		The initial guess
	'''
	p0=[6.3,2.]    ;    p04=[5.74, 1.14]
	''' x range 
	'''
	xrange = linspace(0, 2.5, 200)
	''' '--------------------------------------------------------------------------------------------------------------------------'
	Fit the bias to a fucntion
	'''
	#y0 = find_parameters(p0, np.array(z_[0]), np.array(b_[0]) , 0, xrange)
	y_ = []
	for i in range(len(Srms)):
		y_.append(find_parameters(p0, np.array(z_[i]), np.array(b_[i]) , Srms[i], xrange) )
	#print y_	
	
	'''plot using cmap
	'''
	y = y_ ; z= b_
	tic = [ r'   $0 \mu$Jy',r'   $1 \mu$Jy', r'   $3 \mu$Jy', r'   $5 \mu$Jy',r'   $6 \mu$Jy',r'   $ 7.3 \mu$Jy' , \
	'   $10 \mu$Jy' ,  r'   $23\mu$Jy', r'   $40 \mu$Jy' , r'   $70\mu$Jy'  ,r'   $100\mu$Jy' , r'  $200\mu$Jy'] 
	p_tic =[0, 1, 3, 5, 6, 7.3,10, 23, 40, 70, 100,200]
	#tic = [ r'   $0 \mu$Jy',r'   $1 \mu$Jy', r'   $3 \mu$Jy', r'   $5 \mu$Jy', \
	#r'   $10 \mu$Jy' ,  r'   $23\mu$Jy' ,r'   $100\mu$Jy' ] 
	#p_tic =[0, 1, 3, 5,10, 23, 100]
	
	#=========plot
	fig, ax= plt.subplots()
	n = len(y_)
	for i in range(n):
	    color = i / float(n)	
	    heatmap=colorline(xrange, y_[i], color, cmap="RdBu" , linewidth=1.5)
    
	#============== sidebar    
	cbar = plt.colorbar(heatmap, ticks=())
	cbar.ax.set_yticklabels(tic)
	for j, lab in enumerate(tic):
	    cbar.ax.text(.6, (1* j) /float(n), lab, va='center')
	colors = ['#700f2c', '#9f1128' ,'#b6495b' , '#e0765d',  '#fbdccf', '#fcd7c2', '#fbf2ec','#edf3f7',  '#c7e0ee',  '#96c7df', '#6faed2', 'blue']
	#colors = ['#700f2c', '#9f1128' ,'#b6495b' , '#e0765d', '#fbf2ec','#edf3f7', '#6faed2']
	#========plot the data
	for i in range (len(z_)):
		   ax.scatter(z_[i], b_[i],s= 35, marker= 'o', edgecolor =colors[i], facecolor=colors[i])		
	''''--------------------------------------------------------------------------------------------------------------------------'
	y axis
	'''
	plt.ylabel(r"$b(z)$", fontsize=15)
	plt.ylim(0., 5)
	''''--------------------------------------------------------------------------------------------------------------------------'
	x axis
	''' 
	plt.xlabel(r"$ {\rm redshift} (z)$", fontsize=15)
	plt.xlim(0.01,2.2)
	''''--------------------------------------------------------------------------------------------------------------------------
	'''
	plt.savefig('outputs/fitted_bias.pdf')
	plt.show()
	print '---------------------- Program excuted successfully--------------------'
	print '----------------------------------The plot is saved in the outputs folder Thanks! -------------------------------------'
