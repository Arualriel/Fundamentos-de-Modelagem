#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 11 15:16:47 2020

@author: laura
"""


######## Bibliotecas #########

import numpy as np  # biblioteca numérica
import time # biblioteca para marcar o tempo
import matplotlib.pyplot as plt  # biblioteca gráfica

##############################


######## Constantes globais ########

global lamb1, deltah1,deltav1,gama1,mu1, media1, a1,f1,T1,Nh

global lamb2, deltah2,deltavi1,deltava1,gamai1,gamaa1,mu2, media2, a2,f2,T2

# #belford roxo sir+sis

# lamb1, deltah1,deltav1=5.42053196e+00,2.71887729e-02, 2.16177270e+03
# gama1,mu1= 6.55815765e+01,  1.88298626e+01
# media1, a1,f1,p1=1.46559090e+01,1.15112703e+01, 1.25813791e-01, 5.89995457e-01


# #belford roxo siar+sis

# lamb2, deltah2,deltavi1,deltava1=5.39772934e+00,2.91835124e-02, 2.17941095e+03, 2.16687948e+03
# gamai1,gamaa1,mu2= 6.53618105e+01, 6.50996308e+01, 1.87711807e+01
# media2, a2,f2,p2=1.46924679e+01,1.16615355e+01, 1.30569397e-01, 5.79489691e-01


#niteroi sir+sis

lamb1, deltah1,deltav1=3.06254644e+00,1.65112331e-02, 2.26643008e+03
gama1,mu1=  5.28295050e+01,1.58458478e+01
media1, a1,f1,p1= 1.62872421e+01,9.50314350e+00, 5.11436138e-02, 9.29835108e-01


#niteroi siar+sis

lamb2, deltah2,deltavi1,deltava1=3.06931477e+00,3.06621910e-02, 2.42215762e+03, 2.66414023e+03
gamai1,gamaa1,mu2= 5.61770018e+01, 5.76311087e+01, 1.67802625e+01
media2, a2,f2,p2= 1.41902316e+01,1.06372618e+01, 5.38594381e-02, 8.72882952e-01

global rho
rho=1.0/4.0

Nh=1.0

####################################


######## Funções ########

def f(T,n,t): # função a qual queremos achar o zero
    if(n%2==0):
        sh=T[0]
        ih=T[1]
        ah=T[2]
        rh=T[3]
        sv=T[4]
        iv=T[5]
        s1=deltah2*Nh-lamb2*iv*sh-deltah2*sh
        s2=rho*lamb2*iv*sh-deltah2*ih-gamai1*ih
        s3=(1.0-rho)*lamb2*iv*sh-deltah2*ah-gamaa1*ah
        s4=gamai1*ih+gamaa1*ah-deltah2*rh
        s5=(media2+a2*np.cos(2.0*np.pi*(t-f2)/p2))-deltavi1*sv*ih-deltava1*sv*ah-mu2*sv
        s6=deltavi1*ih*sv+deltava1*ah*sv-mu2*iv
        F=np.zeros(n)
        F[0],F[1],F[2],F[3],F[4],F[5]=s1,s2,s3,s4,s5,s6
    else:
        sh=T[0]
        ih=T[1]
        rh=T[2]
        sv=T[3]
        iv=T[4]
        s1=deltah1*Nh-lamb1*iv*sh-deltah1*sh
        s2=lamb1*iv*sh-deltah1*ih-gama1*ih
        s3=gama1*ih-deltah1*rh
        s4=media1+a1*np.cos(2.0*np.pi*((t-f1)/p1))-deltav1*sv*ih-mu1*sv
        s5=deltav1*ih*sv-mu1*iv
        F=np.zeros(n)
        print("aqui",s1,s2,s3,s4,s5,t)
        F[0],F[1],F[2],F[3],F[4]=s1,s2,s3,s4,s5
    return F

def jacobianaf(T,n): # jacobiana de f
    J=np.zeros((n,n))
    if(n%2==0):
        sh=T[0]
        ih=T[1]
        ah=T[2]
        rh=T[3]
        sv=T[4]
        iv=T[5]
        J[0,0],J[0,1],J[0,2],J[0,3],J[0,4],J[0,5]=-lamb2*iv-deltah2,0.0,0.0,0.0,0.0,-lamb2*sh
        J[1,0],J[1,1],J[1,2],J[1,3],J[1,4],J[1,5]=rho*lamb2*iv,-deltah2-gamai1,0.0,0.0,0.0,rho*lamb2*sh
        J[2,0],J[2,1],J[2,2],J[2,3],J[2,4],J[2,5]=(1.0-rho)*lamb2*iv,0.0,-deltah2-gamaa1,0.0,0.0,(1.0-rho)*lamb2*sh
        J[3,0],J[3,1],J[3,2],J[3,3],J[3,4],J[3,5]=0.0,gamai1,gamaa1,-deltah2,0.0,0.0
        J[4,0],J[4,1],J[4,2],J[4,3],J[4,4],J[4,5]=0.0,-deltavi1*sv,-deltava1*sv,0.0,-deltavi1*ih-deltava1*ah-mu2,0.0
        J[5,0],J[5,1],J[5,2],J[5,3],J[5,4],J[5,5]=0.0,deltavi1*sv,deltava1*sv,0.0,deltavi1*ih+deltava1*ah,-mu2
    else:
        sh=T[0]
        ih=T[1]
        rh=T[2]
        sv=T[3]
        iv=T[4]
        J[0,0],J[0,1],J[0,2],J[0,3],J[0,4]=-lamb1*iv-deltah1,0.0,0.0,0.0,-lamb1*sh
        J[1,0],J[1,1],J[1,2],J[1,3],J[1,4]=lamb1*iv,-deltah1-gama1,0.0,0.0,lamb1*sh
        J[2,0],J[2,1],J[2,2],J[2,3],J[2,4]=0.0,gama1,-deltah1,0.0,0.0
        J[3,0],J[3,1],J[3,2],J[3,3],J[3,4]=0.0,-deltav1*sv,0.0,-deltav1*ih-mu1,0.0
        J[4,0],J[4,1],J[4,2],J[4,3],J[4,4]=0.0,deltav1*sv,0.0,deltav1*ih,-mu1
        
    return J



def metododejacobi(J,b,n): # método iterativo de Jacobi
    l=1
    Nj=100
    t0=np.zeros(n)
    t=np.zeros(n)
    ti=t0
    TOL=10.0**(-8.0) # tolerância do método de Jacobi
    while(l<=Nj):
        for i in range(n):
            soma=0.0
            for j in range(n):
                if(j!=i):
                    soma=soma+J[i,j]*t0[j]
            t=-t0
            ti[i]=(1.0/J[i,i])*(-soma+b[i])
            t=t+ti
            print("Teta(",i+1,")=",ti)
            norma=0.0
            for k in range(n):
                norma=norma+t[k]**2.0
            norma=(norma)**(1.0/2.0)
            print("||Teta(",i+1,")-Teta(",i,")||=",norma)
        if(norma<TOL):
            print("JACOBI")
            print("Convergiu depois de ",l,"iterações.")
            return ti
        l=l+1
        t0=ti
    if(l>=Nj):
        print("Número máximo de iterações atingido.")
        return ti
 

def metododenewton(T0,tol,N,n): # método de Newton 
    l=1
    F=np.zeros(n)
    J=np.zeros((n,n))
    Ge=np.zeros(N)
    T=np.zeros(n)
    Y=np.zeros(n)
    for i in range(n):
        T[i]=T0[i]
    Y=T
    while(l<N):
        print("Teta(",l-1,") = ",T)
        F=f(T,n,l)
        for i in range(n):
            F[i]=-F[i]
        J=jacobianaf(T,n)
        print("determinante",np.linalg.det(J))
        print("Jf(T(",l-1,"))=",J)
        print("f(T(",l-1,"))=",F)
        Y=metododejacobi(J,F,n)
        T=T+Y
        print("Y(",l-1,") = ",Y)
        norma=0.0
        for i in range(n):
            norma=norma+Y[i]**2.0
        norma=(norma)**(1.0/2.0)
        Ge[l-1]=norma
        print("||Y(",l-1,")|| = ",norma)
        if(norma<tol):
            print("NEWTON")
            print("T(",l-1,")=",T,"É a aproximação de zero!!!")
            print("f(T)=",f(T,n,l))
            return T, l, Ge,J
        l=l+1
    if(l>=N):
        print("Número máximo de iterações excedido.")
        print("T(",l-1,")=",T,"É a aproximação de zero!!!")
        print("f(T)=",f(T,n,l))
        return T, l-1, Ge,J

tol=10.0**(-2.0) # tolerância

N=100# número máximo de iterações

n=6 # 5 6
T0=np.ones(n) # ponto inicial
T1=np.zeros(n)
J=np.zeros((n,n))
E1=np.zeros(N)

for i in range(n):
    T0[i]=T0[i]*0.05
l1=0

t=time.time()
T1,l1,E1,J=metododenewton(T0,tol,N,n)
print("tempo de execução",time.time()-t)
print("Convergiu depois de ",l1," iterações.")
print("jacobiana",J)
print("autovalores de J",np.linalg.eigvals(J))

Ge1=np.zeros(l1)

I1=np.zeros(l1)


for i in range(l1):
    Ge1[i]=np.log(E1[i])
    I1[i]=i

#########################################


######## Plotando os gráficos ########

pontos1=plt.scatter(I1,Ge1)
linha1,=plt.plot(I1,Ge1,label="ln(||X(l)-X(l-1)||)",color='purple',ls='-')
plt.xlabel("l")
plt.ylabel("ln(||X(l)-X(l-1)||)")
plt.legend(handles=[linha1,pontos1],loc='best')
plt.show() # plotando o log do erro para o método de Newton


# ######################################