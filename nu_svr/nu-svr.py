from math import sqrt
import numpy as np
import matplotlib.pyplot as plt
from numpy import diag, matrix, inf
from openopt import QP
import math

def kernel_value(x,y):
	a=math.exp(-1*abs(x-y)**2)
	return a

def product(a,X,x):
	prod=0.0
	for i in range(len(a)):
		prod=prod+a[i]*kernel_value(X[i],x)
	return prod

#eps=0
#C=1.0
nu=0.5
C=10000.0

X=[]
Y=[]

#taking input from the file
tot_values=int(raw_input())
for i in range(tot_values):
	X.append(float(raw_input()))
	Y.append(float(raw_input()))

#H=kernel matrix
kernel=[[0.0 for i in range(2*tot_values)] for j in range(2*tot_values)]

for i in range(tot_values):
	for j in range(tot_values):
		kernel[i][j]=kernel_value(X[i],X[j])
		kernel[i+tot_values][j+tot_values]=kernel_value(X[i],X[j])
#negating the values for an cap
for i in range(tot_values):
	for j in range(tot_values):
		kernel[i+tot_values][j]=(-1.0)*kernel_value(X[i],X[j])
		kernel[i][j+tot_values]=(-1.0)*kernel_value(X[i],X[j])


#coeff of 2nd term to minimize
f=[0.0 for i in range(2*tot_values)]
for i in range(tot_values):
	f[i]=-float(Y[i])
for i in range(tot_values,2*tot_values):
	f[i]=float(Y[i-tot_values])


#constraints
lower_limit=[0.0 for i in range(2*tot_values)]
upper_limit=[float(C)/tot_values for i in range(2*tot_values)]
Aeq = [1.0 for i in range(2*tot_values)]
for i in range(tot_values,2*tot_values):
	Aeq[i]=-1.0
beq=0.0
A=[1.0 for i in range(2*tot_values)]
b=nu*float(C)





#coeff for 3rd constraint
#kernel=H
p = QP(np.asmatrix(kernel),np.asmatrix(f),lb=np.asmatrix(lower_limit),ub=np.asmatrix(upper_limit),Aeq=Aeq,beq=beq,A=A,b=b)
# or p = QP(H=diag([1,2,3]), f=[15,8,80], A= ...)
r = p._solve('cvxopt_qp', iprint = 0)
f_opt, x = r.ff, r.xf


support_vectors=[]
support_vectors_Y=[]
coeff=[]
b=0.0
#support vectors: points such that an-an' ! = 0
for i in range(tot_values):
	if not((x[i]-x[tot_values+i])==0):
		support_vectors.append( X[i] )
		support_vectors_Y.append(Y[i])
		coeff.append( x[i]-x[tot_values+i] )

#bias_term=tn-eps-(support vectors)*corresponding kernel
support_vector=[]
support_vector_Y=[]
low=min(abs(x))
for i in range(tot_values):
	if not(abs(x[i]-x[tot_values+i])<low+0.005):
		support_vector.append( X[i] )
		support_vector_Y.append(Y[i])

bias=0.0
for i in range(len(X)):
	bias=bias+float(Y[i]-product(coeff,support_vectors,X[i]))
bias=bias/len(X)


output_X=[]
output_Y=[]

output_X.append(0.0)

for i in range(350):
	output_X.append(output_X[-1]+float(10)/300)

for i in output_X:
	output_Y.append(product(coeff,support_vectors,i)+b)

plt.scatter(output_X,output_Y,marker='o')
plt.scatter(X,Y,marker='.')
print (x)
plt.scatter(support_vector,support_vector_Y,c="yellow",marker='x')

plt.show()

