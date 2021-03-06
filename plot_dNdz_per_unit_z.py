'''
This Program to fit dndz for different rms sensitvities,
produce the fitting values and plot
'''
import numpy as np
import scipy.optimize as opt
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection


def func_(p,x):
	''' This function is the function we think it will be easy
	to fit its parameters to our data
	'''
	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x)
	print w.size
	return w
	
	

def residuals_(p,x,y):
	'''The purpose of this function is finding the 
	difference between the theortical and the simulated data point 
	at specific point x (or redshift).
	'''
	w=10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x)
	err=w-y
	err=err**2
	B=sum(err)
	return B
def find_parameters_(p0, x, rms , rms_value, xrange):
	'''This function call opt.fmin function to  minimize the error on
	you paramters and give you the best fit parameters
	'''
	plsqtot = opt.fmin(residuals_, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value, 'p[0] = ' ,plsqtot[0],  'p[1] = ' , plsqtot[1],'p[2] = ', plsqtot[2]
	print '========================================'
	ytot=func_(plsqtot,xrange)
	return ytot, xrange
#==============================================================================
def func(p,x):
	''' This function is the function we think it will be easy
	to fit its parameters to our data
	'''
	w=10.**p[0]*  x**p[1]  * np.exp(-p[2]*x); w=np.log(w) 
	print w.size
	return w
	
	

def residuals(p,x,y):
	'''The purpose of this function is finding the 
	difference between the theortical and the simulated data point 
	at specific point x (or redshift).
	'''
	w=np.log(10.**p[0]* x**(p[1] ) *  np.exp(-p[2]*x))
	err=w-y
	err=err**2
	B=sum(err)
	return B

def find_parameters(p0, x, rms , rms_value, xrange):
	'''This function call opt.fmin function to  minimize the error on
	you paramters and give you the best fit parameters
	'''
	plsqtot = opt.fmin(residuals, p0, args=(x,rms), maxiter=10000, maxfun=10000)
	print  'rms = %d uJy' % rms_value,'p[0] =' ,plsqtot[0],  'p[1] =' , plsqtot[1], 'p[2]=', plsqtot[2]
	print '========================================'
	ytot=func(plsqtot,xrange)
	return ytot, xrange
def make_segments(x, y):
    '''
    Create list of line segments from x and y coordinates, in the correct format for LineCollection:
    an array of the form   numlines x (points per line) x 2 (x and y) array
    ''' 
    points = np.array([x, y]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)
    return segments

def colorline(x, y, z=None, cmap=plt.get_cmap('hot'), norm=plt.Normalize(0.0, 1.0), linewidth=2.0, alpha=1.0):
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
        
        
if __name__=='__main__':        
	'''
	Read your file where dN/dz [deg^-2 per unit z] are stored
	'''
	(x2, to, rm00muJy,rm01muJy, rm03muJy, rm05muJy, rm06muJy, rm073muJy, rm010muJy, rm023muJy, rm040muJy, rm070muJy, rm100muJy, rm150muJy, rm200muJy) = np.loadtxt('inputs/HIdndzb_modified_high.txt', unpack=True)
	(x1, to, rm0muJy,rm1muJy, rm3muJy, rm5muJy, rm6muJy, rm73muJy, rm10muJy, rm23muJy, rm40muJy, rm70muJy, rm0100muJy, rm0150muJy, rm0200muJy) = np.loadtxt('inputs/HIdndzb_modified.txt', unpack=True)
	(x3, tot, rs0muJy,rs1muJy, rs3muJy, rs5muJy, rs6muJy, rs73muJy, rs10muJy, rs23muJy, rs40muJy, rs70muJy, rs100muJy, rs150muJy, rs200muJy) = np.loadtxt('inputs/HIdndzb3.txt', unpack=True)
	(x, total, rms0muJy,rms1muJy, rms3muJy, rms5muJy, rms6muJy, rms73muJy, rms10muJy, rms23muJy, rms40muJy, rms70muJy, rms100muJy, rms150muJy, rms200muJy) = np.loadtxt('inputs/HIdndzb3_corrected.txt', unpack=True)
	#=========================================================================================
	
	''' p0 and p04 are the intial guess for your parameters 
	In this case its 3 parameters.
	'''
	p0=[5.52,  0.6, 4.6]
	p04=[5.74, 1.14, 3.95]

	'''Define x axis range (or redshift range)
	'''
	xrange = np.linspace(0, 3.0, 200)

	'''Call the functions to fit the data and get the best fit parameters
	'''
	print ' |   c1  | ',       '|         c2  |',        '|         c3  |'
	(ytot,xrange)= find_parameters_(p0,  x,total,0,xrange)
	(y0,xrange)=  find_parameters_(p0,  x,rms0muJy, 0,xrange)
	(y1,xrange) =  find_parameters_( p04, x1,rm1muJy, 1,xrange)
	(y3,xrange)= find_parameters_( p0, x1,rm3muJy, 3, xrange)
	(y5,xrange)= find_parameters_( p0, x1,rm5muJy ,5, xrange)
	(y6,xrange)= find_parameters_( p0, x1,rm6muJy, 6, xrange)
	(y7,xrange)= find_parameters_(p0, x1,rm73muJy, 7, xrange)
	(y10,xrange)= find_parameters_( p0, x1,rm10muJy, 10, xrange)
	(y23,xrange)= find_parameters_( p0, x1,rm23muJy, 23, xrange)
	(y40,xrange)= find_parameters_(p0, x1, rm40muJy, 40, xrange)
	(y70,xrange)= find_parameters_( p0, x1,rm70muJy, 70, xrange)
	(y100,xrange) = find_parameters_(p04, x2,rm100muJy,100, xrange)
	(y150,xrange) = find_parameters_(p04, x2,rm150muJy,150, xrange)
	(y200,xrange) =find_parameters_( p04, x2,rm200muJy,200, xrange)


# fit on the log space  fit log(dn/dz) to log(data)
	(yt,xrange)= find_parameters(p0,  x,np.log(total),0,xrange)
	(y0t,xrange)=  find_parameters(p0,  x,np.log(rms0muJy), 0,xrange)
	(y1t,xrange) =  find_parameters( p0, x,np.log(rms1muJy), 1,xrange)
	(y3t,xrange)= find_parameters( p0, x,np.log(rms3muJy), 3, xrange)
	(y5t,xrange)= find_parameters( p0, x,np.log(rms5muJy) ,5, xrange)
	(y6t,xrange)= find_parameters( p0, x,np.log(rms6muJy), 6, xrange)
	(y7t,xrange)= find_parameters(p0, x,np.log(rms73muJy), 7, xrange)
	(y10t,xrange)= find_parameters(p0, x,np.log(rms10muJy), 10, xrange)
	(y23t,xrange)= find_parameters_(p0, x1,rm23muJy, 23, xrange)
	(y40t,xrange)= find_parameters_(p0, x1, rm40muJy, 40, xrange)
	(y70t,xrange)= find_parameters_(p0, x1,rm70muJy, 70, xrange)
	(y100t,xrange) = find_parameters_(p04, x2,rm100muJy,100, xrange)
	(y150t,xrange) = find_parameters_(p04, x2,rm150muJy,150, xrange)
	(y200t,xrange) =find_parameters_(p04, x2,rm200muJy,200, xrange)
	print '============ Program excuted successfully ==========='
	'''
	Plot the results from this program 
	'''
	fig = plt.figure()
	ax = fig.add_subplot(1,1,1)
	ax.set_yscale('log')

	#y = log([y0, y1, y3, y5, y6, y7, y10, y23, y40, y70, y100, y150, y200])
	#yt = ([np.exp(y0t), np.exp(y1t), np.exp(y3t), np.exp(y5t), np.exp(y6t), np.exp(y7t), np.exp(y10t), y23t, y40t, y70t, y100t, y150t, y200t])
	yt = ([np.exp(y0t), np.exp(y1t), np.exp(y3t), np.exp(y5t), np.exp(y10t), y23t, y100t, y200t])
	#tic = [r'  $  \ 0 \mu$Jy',r'  $  \ 1 \mu$Jy', r'  $ \ 3 \mu$Jy', r'  $ \ 5 \mu$Jy',r'  $ \ 6 \mu$Jy',r'  $ \ 7.3 \mu$Jy' ,r'  $ \ 10 \mu$Jy' ,  r'  $ \ 23\mu$Jy',  r'  $  \ 40 \mu$Jy' , r'  $  \ 70\mu$Jy'  ,r'  $  \ 100\mu$Jy'  ,r'  $ \ 150\mu$Jy',  r'  $ \ 200\mu$Jy']

	tic = [r'  $  \ 0 \mu$Jy', r'  $ \ 1 \mu$Jy', r'  $ \ 3 \mu$Jy', r'  $ \ 5 \mu$Jy', r'  $  \ 10\mu$Jy' , r'  $ \ 23\mu$Jy' ,r'  $  \ 100\mu$Jy'  ,  r'  $ \ 200\mu$Jy']
	#=========plot
	fig, ax= plt.subplots()
	#ax.set_yscale('log')
	#n = len(y)
	#for i in range(n):
    	#	color = i / float(n)	
    	#	heatmap=colorline(xrange, y[i], color, cmap="RdBu")
    		
    	#fig, ax= plt.subplots()
	ax.set_yscale('log')
	nt = len(yt)
	for i in range(nt):
    		color = i / float(nt)	
    		heatmap2=colorline(xrange, yt[i], color, cmap="RdBu")
	#=========== 

	ax.scatter(x,(rms0muJy), s= 35, marker= 'o', edgecolor ='#71112e', facecolor='#71112e')
	ax.scatter(x,(rms1muJy), s= 35,marker ='o',edgecolor = '#c13d45', facecolor='#c13d45')
	ax.scatter(x,(rms3muJy) , s=35, marker = 'o', edgecolor = '#e68368', facecolor='#e68368')
	ax.scatter(x,(rms5muJy),  s= 35,marker ='o',edgecolor = '#fde3d5', facecolor='#fde3d5')
	#ax.scatter(x,(rms6muJy) ,  s= 35,marker ='o',edgecolor = '#fbdccf', facecolor='#fbdccf')
	#ax.scatter(x,(rms73muJy) ,  s= 35,marker ='o',edgecolor = '#fcd7c2', facecolor='#fcd7c2')
	ax.scatter(x,(rms10muJy),  s= 35,marker ='o',edgecolor = '#f9f9f9', facecolor='#f9f9f9')
	ax.scatter(x,(rms23muJy) ,  s= 35,marker ='o',edgecolor = '#c5e0ed', facecolor='#c5e0ed')
	#ax.scatter(x,(rms40muJy),  s= 35,marker ='o',edgecolor = '#c7e0ee', facecolor='#c7e0ee')
	#ax.scatter(x,(rms70muJy),  s= 35,marker ='o',edgecolor = '#96c7df', facecolor='#96c7df')
	ax.scatter(x,(rms100muJy) , s= 35,marker ='o',edgecolor = '#74b2d4', facecolor='#74b2d4')
	#ax.scatter(x,(rms150muJy),  s= 35,marker ='o',edgecolor = '#337eb8', facecolor='#337eb8')
	ax.scatter(x,(rms200muJy),  s= 35,marker ='o',edgecolor = '#2870b1', facecolor='#2870b1')    
	
	   
	#============== sidebar    
	cbar = plt.colorbar(heatmap2, ticks=())
	cbar.ax.set_yticklabels(tic)
	for j, lab in enumerate(tic):
	    cbar.ax.text(.6, (1* j) /float(nt), lab, va='center')
	
	plt.xlim(0.1,2.2 ,0.2)
	plt.ylim(1, 5e6)
	#plt.xlim(0.1,2.5)
	#plt.ylim(-5, 15)	
	
	
	#================ x axis 
	xticks = np.arange(min(x), max(x)+1, 0.3)
	plt.xlabel(r"${ \rm redshift} (z)$", fontsize=15)
	#============= y axis
	#yticks = [1, 10,100, 1e3, 1e4, 1e5, 1e6]
	#plt.yticks(yticks,[r'$1$', r'$10$',r'$10^2$', r'$10^3$',r'$10^4$',r'$10^5$',r'$10^6$'])
	plt.ylabel(r'$\frac{d{\rm N}}{dz}(z) \ [ {\rm deg}^{-2} \ {\rm per} \ {\rm unit} \ z ]$', fontsize= 15)
	#========= save fig 
	plt.savefig('outputs/fittingMario_dNOverdz_using_ObreschkowFunc.pdf')
	plt.show()
	print '==================The plot is save in the outputs folder Thanks! ======================'
