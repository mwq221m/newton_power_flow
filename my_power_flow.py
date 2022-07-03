# -*- coding: utf-8 -*-
"""
Created on Wed Jul 14 10:20:54 2021

@author: ma
"""

from numpy import *
from numpy.linalg import *
import time
start=time.time()
def f_calculation(u,det,G,B,idx):
    length=shape(G)[0]
    temp1=0;temp2=0
    for j in range(length):
        temp1+=u[j]*(G[idx][j]*cos(det[idx]-det[j])+B[idx][j]*sin(det[idx]-det[j]))
    p=temp1*u[idx]
    for j in range(length):
        temp2+=u[j]*(G[idx][j]*sin(det[idx]-det[j])-B[idx][j]*cos(det[idx]-det[j]))
    q=temp2*u[idx]
    return (p,q)

def make_H(u,det,G,B,idx1,idx2):
    length=shape(G)[0]
    temp=0
    if idx1==idx2:
        for j in range(length):
            if j!=idx1:
                temp+=u[j]*(G[idx1][j]*sin(det[idx1]-det[j])-B[idx1][j]*cos(det[idx1]-det[j]))
        temp=temp*(-u[idx1])
    else:
        temp=u[idx1]*u[idx2]*(G[idx1][idx2]*sin(det[idx1]-det[idx2])-B[idx1][idx2]*cos(det[idx1]-det[idx2]))
    return temp

def make_N(u,det,G,B,idx1,idx2):
    length=shape(G)[0]
    temp=0
    if idx1==idx2:
        for j in range(length):
            if j!=idx1:
                temp+=u[j]*(G[idx1][j]*cos(det[idx1]-det[j])+B[idx1][j]*sin(det[idx1]-det[j]))
        temp=temp*u[idx1]+2*u[idx1]**2*G[idx1][idx1]
    else:
        temp=u[idx1]*u[idx2]*(G[idx1][idx2]*cos(det[idx1]-det[idx2])+B[idx1][idx2]*sin(det[idx1]-det[idx2]))
    return temp

def make_J(u,det,G,B,idx1,idx2):
    length=shape(G)[0]
    temp=0
    if idx1==idx2:
        for j in range(length):
            if j!=idx1:
                temp+=u[j]*(G[idx1][j]*cos(det[idx1]-det[j])+B[idx1][j]*sin(det[idx1]-det[j]))
        temp=temp*u[idx1]
    else:
        temp=-u[idx1]*u[idx2]*(G[idx1][idx2]*cos(det[idx1]-det[idx2])+B[idx1][idx2]*sin(det[idx1]-det[idx2]))
    return temp



def make_L(u,det,G,B,idx1,idx2):
    length=shape(G)[0]
    temp=0
    if idx1==idx2:
        for j in range(length):
            if j!=idx1:
                temp+=u[j]*(G[idx1][j]*sin(det[idx1]-det[j])-B[idx1][j]*cos(det[idx1]-det[j]))
        temp=temp*u[idx1]-2*u[idx1]**2*B[idx1][idx1]
    else:
        temp=u[idx1]*u[idx2]*(G[idx1][idx2]*sin(det[idx1]-det[idx2])-B[idx1][idx2]*cos(det[idx1]-det[idx2]))
    return temp

                
                

def make_Jac(u,det,G,B,information):
    length=len(information)
    Jac=zeros((length,length))
    #temp1,temp2=information[idx1]
    #temp3,temp4=information[idx2]
    for i in range(length):
        for j in range(length):
            temp1,temp2=information[i]
            temp3,temp4=information[j]
            if temp2==1 and temp4==1:
                Jac[i][j]=make_H(u,det,G,B,temp1,temp3)
            if temp2==1 and temp4==2:
                Jac[i][j]=make_N(u,det,G,B,temp1,temp3)
            if temp2==2 and temp4==1:
                Jac[i][j]=make_J(u,det,G,B,temp1,temp3)
            if temp2==2 and temp4==2:
                Jac[i][j]=make_L(u,det,G,B,temp1,temp3)
    return Jac

                
            
    
    
    
    
    
    
class MyPowerFlow(object):
    def __init__(self):
        self.u=[]
        self.det=[]
        self.iter_num=10
        self.Y=array([[]])
        self.Jac=array([[]])
        self.information=[]
        self.f=[]
        
    def add_node(self,node):
        #information=[]
        node_num=shape(node)[0]
        self.u=[1 for i in range(node_num)]
        self.det=[0 for i in range(node_num)]
        self.Y=zeros((node_num,node_num),dtype='complex_')
        for i in range(node_num):
            g=node[i][6];b=node[i][7]
            y=g+1j*b
            self.Y[i][i]+=y
            if node[i][0]==1:
                self.information.append((i,1))
                self.information.append((i,2))
                p=node[i][1]-node[i][3]
                q=node[i][2]-node[i][4]
                self.f.append(p)
                self.f.append(q)
                
            if node[i][0]==2:
                self.information.append((i,1))
                self.u[i]=node[i][5]
                p=node[i][1]-node[i][3]
                self.f.append(p)
            
            if node[i][0]==3:
                self.u[i]=node[i][5]
        
    def add_branch(self,branch):
        branch_num=shape(branch)[0]
        for i in range(branch_num):
            if branch[i][5]==1:
                #print('mark')
                idx1,idx2=branch[i][0],branch[i][1]
                #print(idx1);print(idx2)
                idx1=int(idx1);idx2=int(idx2)
                r,x=branch[i][2],branch[i][3]
                z=r+1j*x
                y=1/z
                b=branch[i][4]
                #print(b)
                #print(y)
                self.Y[idx1][idx1]+=y+b/2
                self.Y[idx2][idx2]+=y+b/2
                self.Y[idx1][idx2]+=-y
                self.Y[idx2][idx1]+=-y
                
    def power_flow(self):
        Jac_num=len(self.information)
        
        
        
        f_now=zeros((Jac_num,1))
        f=array(self.f).reshape(-1,1)
        
        for t in range(self.iter_num):
            self.Jac=make_Jac(self.u,self.det,self.Y.real,self.Y.imag,self.information)
            #vround=vectorize(round)
            #self.Jac=vround(self.Jac,5)
            for i in range(Jac_num):
                idx1,idx2=self.information[i]
                p,q=f_calculation(self.u,self.det,self.Y.real,self.Y.imag,idx1)
                #print(p);print(q)
                if idx2==1:
                    f_now[i][0]=p
                if idx2==2:
                    f_now[i][0]=q
            #print(f_now)
            #print(len(f_now))
            df=f-f_now
            #print(t)
            #print(df)
            dx=inv(self.Jac)@df
            for i in range(Jac_num):
                idx1,idx2=self.information[i]
                if idx2==1:
                    self.det[idx1]+=dx[i][0]
                    #if self.det[idx1]<det_min:
                        #self.det[idx1]=det_min
                    #if self.det[idx1]>det_max:
                        #self.det[idx1]=det_max
                if idx2==2:
                    self.u[idx1]+=dx[i][0]*self.u[idx1]
                    #if self.u[idx1]<u_min:
                        #self.u[idx1]=u_min
                    #if self.u[idx1]>u_max:
                        #self.u[idx1]=u_max
        for i in range(len(self.det)):
            self.det[i]*=180/pi
        
        
                
                
                
            
        
        
                
'''              
node=array([[3,0,0,0,0,1.05,0,0],\
            [2,0.2,0,0.5,0,1.03,0,0],\
                [1,0,0,0.6,0.25,1,0,0]])
branch=array([[0,1,0.08,0.24,0,1],[0,2,0.02,0.06,0,1],[1,2,0.06,0.18,0,1]])
'''



#node=array([],dtype='float_')#加一个小数点 使得矩阵以浮点数存储

'''
node=array([[3,0,0,0,0,1,0,0],\
            [1,0,0,100.0,60,1,0,0],\
                [1,0,0,90,40,1,0,0],\
                    [1,0,0,120,80,1,0,0],\
                        [1,0,0,60,30,1,0,0],\
                            [1,0,0,60,20,1,0,0],\
                                [1,0,0,200,100,1,0,0],\
                                    [1,0,0,200,100,1,0,0],\
                                        [1,0,0,60,20,1,0,0],\
                                            [1,0,0,60,20,1,0,0],\
                                                [1,0,0,45,30,1,0,0],\
                                                    [1,0,0,60,35,1,0,0],\
                                                        [1,0,0,60,35,1,0,0],\
                                                            [1,0,0,120,80,1,0,0],\
                                                                [1,0,0,60,10,1,0,0],\
                                                                    [1,0,0,60,20,1,0,0],\
                                                                        [1,0,0,60,20,1,0,0],\
                                                                            [1,0,0,90,40,1,0,0],\
                                                                               [1,0,0,90,40,1,0,0],\
                                                                                    [1,0,0,90,40,1,0,0],\
                                                                                        [1,0,0,90,40,1,0,0],\
                                                                                            [1,0,0,90,40,1,0,0],\
                                                                                                [1,0,0,90,50,1,0,0],\
                                                                                                    [1,0,0,420,200,1,0,0],\
                                                                                                        [1,0,0,420,200,1,0,0],\
                                                                                                            [1,0,0,60,25,1,0,0],\
                                                                                                                [1,0,0,60,25,1,0,0],\
                                                                                                                    [1,0,0,60,20,1,0,0],\
                                                                                                                        [1,0,0,120,70,1,0,0],\
                                                                                                                            [1,0,0,200,600,1,0,0],\
                                                                                                                                [1,0,0,150,70,1,0,0],\
                                                                                                                                    [1,0,0,210,100,1,0,0],\
                                                                                                                                        [1,0,0,60,40,1,0,0]])








    
#node[:][1:5]=node[:][1:5]/10000
#vfloat=vectorize(float)
#node[:,1:5]=vfloat(node[:,1:5])
node[:,1:5]=node[:,1:5]/10000

branch=array([[0,1,0.0922,0.047,0,1],\
              [1,2,0.493,0.2511,0,1],\
                  [2,3,0.366,0.1864,0,1],\
                      [3,4,0.3811,0.1941,0,1],\
                          [4,5,0.819,0.707,0,1],\
                              [5,6,0.1872,0.6188,0,1],\
                                  [6,7,0.7114,0.2351,0,1],\
                                      [7,8,1.03,0.74,0,1],\
                                          [8,9,1.044,0.74,0,1],\
                                              [9,10,0.1966,0.0650,0,1],\
                                                  [10,11,0.3744,0.1238,0,1],\
                                                      [11,12,1.468,1.155,0,1],\
                                                          [12,13,0.5416,0.7129,0,1],\
                                                              [13,14,0.591,0.526,0,1],\
                                                                  [14,15,0.7463,0.545,0,1],\
                                                                      [15,16,1.289,1.721,0,1],\
                                                                          [16,17,0.372,0.574,0,1],\
                                                                              [1,18,0.164,0.1565,0,1],\
                                                                                  [18,19,1.5042,1.3554,0,1],\
                                                                                      [19,20,0.4095,0.4784,0,1],\
                                                                                          [20,21,0.7089,0.9373,0,1],\
                                                                                              [2,22,0.4512,0.3083,0,1],\
                                                                                                  [22,23,0.898,0.7091,0,1],\
                                                                                                      [23,24,0.896,0.7011,0,1],\
                                                                                                          [5,25,0.203,0.1034,0,1],\
                                                                                                              [25,26,0.2842,0.1447,0,1],\
                                                                                                                  [26,27,1.059,0.9337,0,1],\
                                                                                                                      [27,28,0.8042,0.7006,0,1],\
                                                                                                                          [28,29,0.5075,0.2585,0,1],\
                                                                                                                              [29,30,0.9744,0.963,0,1],\
                                                                                                                                  [30,31,0.3105,0.3619,0,1],\
                                                                                                                                      [31,32,0.341,0.5362,0,1],\
                                                                                                                                          [7,20,2,2,0,0],\
                                                                                                                                              [8,14,2,2,0,0],\
                                                                                                                                                  [11,21,2,2,0,0],\
                                                                                                                                                      [17,32,0.5,0.5,0,0],\
                                                                                                                                                          [24,28,0.5,0.5,0,0]])



 
 

ub=12.66;sb=10
zb=ub**2/sb
branch[:,(2,3)]=branch[:,(2,3)]/zb
'''

node=array([[3,0,0,0,0,1,0,0],\
            [1,0,0,0.02,0.016,1,0,0],\
                [1,0,0,0.04,0.027,1,0,0],\
                    [1,0,0,0.01,0.009,1,0,0],\
                        [1,0,0.25e-2,0.02,0.008,1,0,0],\
                            [1,0,0.002,0.03,0.015,1,0,0],\
                                [1,0,0,0.015,0.012,1,0,0],\
                                    [1,0,0.012,0.05,0.03,1,0,0],\
                                        [1,0,0.037,0.045,0.02,1,0,0],\
                                            [1,0,0.006,0.006,0.001,1,0,0],\
                                                [1,0,0.018,0.01,0.009,1,0,0],\
                                                    [1,0,0,0.01,0.007,1,0,0],\
                                                        [1,0,0,0.01,0.009,1,0,0],\
                                                            [1,0,0.018,0.021,0.01,1,0,0]])
branch=array([[0,1,0.075,0.1,0,1],\
              [0,2,0.11,0.11,0,1],\
                  [0,3,0.11,0.11,0,1],\
                      [1,4,0.09,0.18,0,1],\
                          [1,5,0.08,0.11,0,1],\
                              [4,6,0.04,0.04,0,1],\
                                  [2,7,0.08,0.11,0,1],\
                                      [7,8,0.08,0.11,0,1],\
                                          [7,9,0.11,0.11,0,1],\
                                              [2,10,0.11,0.11,0,1],\
                                                  [3,11,0.09,0.12,0,1],\
                                                      [3,12,0.08,0.11,0,1],\
                                                          [12,13,0.04,0.04,0,1],\
                                                              [5,9,0.04,0.04,0,0],\
                                                                  [10,11,0.04,0.04,0,0],\
                                                                      [6,13,0.09,0.12,0,0]])


            
pf=MyPowerFlow()
pf.add_node(node)
pf.add_branch(branch)
#print(pf.Y)
#print(pf.f)
#print(pf.u)
#print(pf.det)
pf.iter_num=10

pf.power_flow()
print(pf.u)
print(pf.det)
end=time.time()
print('求解使用时间：%ss'%round(end-start,5))

'''
node2=array([[1,0,0,1.5,0.1,1,0,0],[1,0,0,2.5,0.18,1,0,0],[2,5,0,0,0,0.97,0,0],[1,0,0,7,0.2,1,0,0],[1,0,0,6.5,0.15,1,0,0],[3,0,0,0,0,1,0,0]])
branch2=array([[1,0,0.00071,0.00036,0,1],[2,1,0.00017,0.00009,0,1],[3,2,0.00153,0.00051,0,1],[0,3,0.0007,0.00023,0,0],[3,4,0.00035,0.00012,0,1],[0,5,0.00069,0.00035,0,0],[4,5,0.00133,0.00044,0,1]])
pf2=MyPowerFlow()
pf2.add_node(node2)
pf2.add_branch(branch2)
pf2.iter_num=10
pf2.power_flow()
print(pf2.u)
print(pf2.det)
'''







            
                    
                    
                    
        
        