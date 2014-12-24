from pylab import *

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
	#for i in range(len(redshift)):
	#	print redshift[i] , bias[i] , Srms_[i]
	#y = ['y0', 'y1', 'y3', 'y5', 'y6', 'y7', 'y10', 'y23', 'y40', 'y70', 'y100']
	z = ['z0', 'z1', 'z3', 'z5', 'z6', 'z7', 'z10', 'z23', 'z40', 'z70', 'z100', 'z150', 'z200'] 
	b =['b0', 'b1', 'b3', 'b5', 'b6', 'b7', 'b10', 'b23', 'b40', 'b70', 'b100', 'b150', 'b200'] 
	Srms= [0, 1, 3, 5, 6, 7.3, 10, 23,40, 70, 100, 150, 200] 
	#y = ['y0', 'y1', 'y3', 'y5', 'y10', 'y23', 'y100', 'y200']
	#z = ['z0', 'z1', 'z3', 'z5', 'z10', 'z23', 'z100', 'z200'] 
	#b =['b0', 'b1', 'b3', 'b5', 'b10', 'b23',  'b100', 'b200'] 
	#Srms= [0, 1, 3, 5, 10, 23, 100, 200] 
	z_ = [] ; b_ = []
	#print Srms_[0], Srms_[1] , Srms_[2] , Srms_[3] , Srms_[4] , Srms_[5] , Srms_[6]  ,Srms_[7]  ,Srms_[8]  ,Srms_[9]  ,Srms_[10],Srms_[11]    
	for i in range (len(Srms)):
		z[i], b[i] = read_input(Srms[i])
		z_.append(z[i]) ; b_.append(b[i])
	print Srms[0], z[0] , b[0]

			

