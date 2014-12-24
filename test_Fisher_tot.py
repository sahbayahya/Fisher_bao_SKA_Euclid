''' This program to combine the  Fisher matrices (adding) for different bins and give the total fisher matrix.  
'''
from pylab import *
from numpy import linalg

(z, F11, F12, F21, F22) = loadtxt("Fisher_bao_Real_150_srms_corrected_nz_test.txt", unpack=True)

F  = []
for i in range(len(z)):
    	F.append(matrix([[F11[i], F12[i]], [F21[i], F22[i]] ]))
Ftot =  F[0] + F[1] + F[2] + F[3] + F[4]    	
print '_'*30   	
print  'Total F = '  , Ftot

Inv_Ftot = linalg.inv(Ftot)

print 'Err on DA = ' ,sqrt(Inv_Ftot[0,0])*100
print 'Err on H  = ' , sqrt(Inv_Ftot[1,1])*100

    	
