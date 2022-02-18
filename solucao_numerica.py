#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  2 21:35:19 2020

@author: laura
"""

# bibliotecas

import numpy as np
import matplotlib.pyplot as plt
from random import *


global p,deltah,Nh,lamb,gama,deltav,mu,sh0,rh0,ih0,sv0,iv0

deltah=0.1*10.0**(-4.0)
deltav=0.6
Nh=1.0#515317.0 #513118.0
lamb=0.7
gama=0.08
mu=0.02
rh0=0.0
p=1.0/4.0
#ih0=7.0

iv0=0.3144
sv0=1.0-iv0


# funcao f

def F(X,phi):
    sh=X[0]
    ih=X[1]
    ah=X[2]
    rh=X[3] #rh=X[2]
    sv=X[4] #sv=X[3]
    iv=X[5] #iv=X[4]
    Sh = deltah*Nh - lamb*sh*iv - deltah*sh
    Ih = p*lamb*sh*iv - gama*ih - deltah*ih
    Ah = (1-p)*lamb*iv*sh - deltah*ah - gama*ah
    Rh = gama*ih - deltah*rh
    Sv = phi - deltav*sv*ih -mu*sv #+ mu*iv
    Iv = deltav*ih*sv - mu*iv
    
    saida=np.zeros(6)
    saida[0],saida[1],saida[2],saida[3],saida[4],saida[5]=Sh,Ih,Ah,Rh,Sv,Iv
    return saida


### solucao analitica
        


# arquivo de dados





# with open('dados_belford_roxo_casos.txt','r') as arquivo:
#     linhas1 = arquivo.read().splitlines()
    
# tamanho1=len(linhas1)
# tempo1=np.zeros(tamanho1, int)
# dados1=np.zeros(tamanho1)
# cont=0

# j=0
# for i in linhas1:
#     tempo1[j]=(j)*7.0
#     dados1[j]=float(i[0:2])
#     j=j+1
         
# (_, caps, _) = plt.errorbar(
#     tempo1, dados1,fmt='o-', markersize=5, capsize=10)

# plt.show()


with open('dados_niteroi_casos.txt','r') as arquivo:
    linhas2 = arquivo.read().splitlines()
    
tamanho2=len(linhas2)
tempo2=np.zeros(tamanho2, int)
dados2=np.zeros(tamanho2)
cont2=0

print(linhas2)
j=0
for i in linhas2:
    tempo2[j]=j*7.0
    dados2[j]=float(i[0:2])/515317.0
    j=j+1
         
# print("dados niteroi",tempo2,dados2)
# (_, caps, _) = plt.errorbar(
#     tempo2, dados2,fmt='o-', markersize=5, capsize=10)

# plt.show()


###################
tempo=len(tempo2)
num=tempo*7 #tempo total do experimento em semanas
X0=np.zeros(6)
ih0=dados2[0]
ah0=ih0*(1.0/(1-p))
sh0=Nh-ih0-ah0


X0[0],X0[1],X0[2],X0[3],X0[4],X0[5]=sh0,ih0,ah0,rh0,sv0,iv0

td=np.zeros(num)
for i in range(num):
    td[i]=i

# variaveis auxiliares

X1=np.zeros(6) #vetor do passo seguinte

#taxa de variacao do tempo
dt=10.0**(-5.0)
R=np.zeros((6,num)) #matriz de resultados
I=1000
X1=X0
R[0,0]=X0[0]
R[1,0]=X0[1]
R[2,0]=X0[2]
R[3,0]=X0[3]
R[4,0]=X0[4]
R[5,0]=X0[5]
phi=0.0
epsilon=10.0**(-15.0)
E=1.0
a=0.0
f=0.0
B=np.zeros((2,I))
erro=np.zeros(I)
cont=0
seed()
while((E>=epsilon)and(cont<I)):
    a=uniform(1,3)
    f=uniform(1,3)
    B[0,cont]=a
    B[1,cont]=f
    X0[0],X0[1],X0[2],X0[3],X0[4]=sh0,ih0,rh0,sv0,iv0
    X1=X0
    soma=0.0
    for k in range(num):
        phi=a*np.sin(2*np.pi*(k-f)/num)
        K1=F(X0,phi)
        K2=F(X0+dt*K1*0.5,phi)
        K3=F(X0+dt*K2*0.5,phi)
        K4=F(X0+dt*K3,phi)
        X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
       
        if(k<num-1):         
            R[0,k+1]=X1[0]
            R[1,k+1]=X1[1]
            R[2,k+1]=X1[2]
            R[3,k+1]=X1[3]
            R[4,k+1]=X1[4]
            R[5,k+1]=X1[5]
        X0=X1
    for k in range(len(tempo2)):
        soma=soma+(dados2[k]-R[1,tempo2[k]])**2.0
        
        
    erro[cont]=soma
    E=erro[cont]
    cont=cont+1
menor=erro[0]
print("erros e menor",erro,menor)
ind=0
for k in range(cont):
    if(erro[k]<=menor):
        menor=erro[k]
        ind=k
abest=B[0,ind]    
fbest=B[1,ind]
print("aaaaa",abest,fbest)
for k in range(num):
    phi=abest*np.sin(2*np.pi*(k-fbest)/num)
    K1=F(X0,phi)
    K2=F(X0+dt*K1*0.5,phi)
    K3=F(X0+dt*K2*0.5,phi)
    K4=F(X0+dt*K3,phi)
    X1=X0+(dt/6.0)*(K1+2*K2+2*K3+K4)
       
    if(k<num-1):         
        R[0,k+1]=X1[0]
        R[1,k+1]=X1[1]
        R[2,k+1]=X1[2]
        R[3,k+1]=X1[3]
        R[4,k+1]=X1[4]
        R[5,k+1]=X1[5]
    X0=X1




x=R[0,:]
y=R[1,:]
z=R[2,:]
w=R[3,:]
xv=R[4,:]
yv=R[5,:]   
print("R",R)
linha1,=plt.plot(td,x,label="Sh(t)",color='red',ls='-.')

linha2,=plt.plot(td,y,label="Ih(t)",color='blue',ls='-.')

linha3,=plt.plot(td,z,label="Ah(t)",color='pink',ls='-.')

linha4,=plt.plot(td,w,label="Rh(t)",color='green',ls='-.')

linha5,=plt.plot(td,xv,label="Sv(t)",color='orange',ls='-')

linha6,=plt.plot(td,yv,label="Iv(t)",color='purple',ls='-')



plt.xlabel("t")
plt.ylabel("Dinâmica(t)")
plt.legend(handles=[linha1,linha2,linha3,linha4,linha5,linha6],loc='best')

plt.show() 


linha2,=plt.plot(td,y,label="I(t)",color='blue',ls='-')

data=plt.scatter(tempo2,dados2,marker='o',color='purple',label="Dados")

plt.xlabel("t")
plt.ylabel("Infectados(t)")
plt.legend(handles=[linha2,data],loc='best')

plt.show() 
