import numpy as np
import math
import time
import baseToolbox as bt
from math import pi
import copy

InteractionsName = ['Tropism','ApicalTropism','Proprioception','ApicalRelative']
CollectiveInteractionsName = ['Apical','Global','Apical_Apical_Morse','Apical_Global_Morse','Global_Apical_Morse','Global_Global_Morse']
GrowthMode = ['Apical','Exponential']



class skeletonElements:
    def __init__(self,x,theta0,curvature0,s,ds,radius,psi0 = 0,psiC0 = 0,psiG0 = 0,dt = .01,L=3,growth = {'name':'no','intensity':1,'direction':0},T_G_noise=1,sig_G_noise=0):
        
        self.x = x
        self.theta = theta0
        self.curvature = curvature0
        self.ds = ds
        self.s = s
        self.radius = radius
        self.psiC = psiC0
        self.psiG = psiG0
        self.psi0 = psi0
        self.dt = dt
        self.growth = growth.copy()
        self.interactions = []
        self.thetaApical = theta0
        self.boxsize=L
        self.istoolong=0
        self.marker = 0
        self.sig_G_noise=sig_G_noise
        self.T_G_noise=T_G_noise
        self.G_noise=float(0)
        
    def updateCurvilinearAbscissa(self,s):
        self.s = s

    def updateApicalAngle(self,theta):
        self.thetaApical = theta

    def updateOrientation(self,theta):
        self.theta = theta  

    def updateSpatialPosition(self,x0,theta0,psi0,ds0):
        self.x[0]  = np.mod(x0[0] + ds0 * (np.sin(theta0) * np.cos(psi0)),self.boxsize)
        self.x[1]  = x0[1] + ds0 * (np.sin(theta0) * np.sin(psi0))
        self.x[2]  = x0[2] + ds0 * (np.cos(theta0))

    def istoolong_update(self,ds_max):
        if self.ds>=2*ds_max:
            self.istoolong=1
        else:
            self.istoolong=0
            
    def set_marker(self):
        self.marker = 1
        
    def addInteractions(self,name,intensity=0,direction = 0):
        
        self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
        
    def update(self,time_counter):
        t0 = time.time()

        delta = 0
        deltaParrallel = 0
        deltaPerpendicular = 0
        for interaction in self.interactions:
            #print('interaction : '+str(interaction['direction']))
            interactionName = interaction['name']
            if interaction['name'] == 'Proprioception':
                deltaParrallel +=interaction['intensity'] * self.curvature 
            else:
                if interaction['name'] == 'Tropism':
                    delta += interaction['intensity'] * np.sin(bt.angleDifference(self.theta,interaction['direction'])) #sin is added!!!!!!
                elif interaction['name'] == 'ApicalTropism':
                    delta += interaction['intensity'] * np.sin(bt.angleDifference(self.thetaApical,interaction['direction'])) #sin is added!!!!!!
                elif interaction['name'] == 'ApicalRelative':
                    delta += interaction['intensity'] * np.sin(interaction['direction']) #sin is added!!!!!!
        #print('angle : '+str(self.thetaApical))
        ###################################################         
        #DELTA MUST BE SMALLER THAN 1. MUST FIX INTESITIES#
        ####################################################
        delta=np.min((delta,1))
        deltaParrallel +=   delta * np.cos(self.psiG - self.psiC)
        deltaPerpendicular +=   delta * np.sin(self.psiG - self.psiC)
        
        
        self.curvature += deltaParrallel * self.growth['growthRate'] * self.dt /self.radius
        if np.abs(self.curvature) >0:
            self.psiC = (deltaPerpendicular * self.growth['growthRate'] * self.dt)/self.curvature
        if self.growth['name'] in GrowthMode:
            if self.growth['growthRate']:
#                if np.mod(time_counter,self.T_G_noise)==1:
#                    self.G_noise= float(np.random.uniform(-self.sig_G_noise,self.sig_G_noise,1))
#                self.ds += self.ds * (self.growth['growthRate'] +self.G_noise)* self.dt 
                self.ds += self.ds * (self.growth['growthRate'])* self.dt 


class Plant:
    def __init__(self,x0,theta0=0,curvature0=0,length0 = 1,radius=0.01,psi0 = 0,psiC0 = 0,psiG0 = 0,ds_max=0.01,dt=.01,L=15,growth = 'no',growthZone = 1,growthRate = 1,T_G_noise=1,sig_G_noise=0):
        
        self.x = []
        self.s = []
        self.theta = []
        self.curvature = [] 
        
        self.boxsize=L
        self.x0 = x0
        self.theta0 = theta0
        self.curvature0 = curvature0
        
        self.radius=radius
        self.length0 = length0
        self.length = length0
        self.ds = ds_max
        self.ds_max=ds_max
        self.psi0 = psi0
        self.psiG0 = psiG0
        self.psiC0 = psiC0
        self.dt = dt
        self.skeleton = []
        self.skeleton_copy=[]
        
        self.interactions = []
        x = self.x0
        theta = self.theta0
        curvature = self.curvature0

        self.growth = {'name':growth,'growthZone':growthZone,'growthRate':growthRate}
        self.T_G_noise=T_G_noise
        self.sig_G_noise=sig_G_noise
        self.G_noise=float(0)

        s = 0

        self.skeleton.append(skeletonElements(np.copy(x),theta,curvature,0,self.ds,self.radius,self.psi0,psiC0,psiG0,self.dt,self.boxsize,growth=self.growth,T_G_noise=self.T_G_noise,sig_G_noise=self.sig_G_noise))
        
        
        for k in range(1,int(np.floor(self.length/self.ds))):
            
            x[0]  = np.mod(x0[0] + self.ds * (np.sin(theta0) * np.cos(psi0)),self.boxsize)
            x[1]  = x[1] + self.ds * (np.sin(theta0) * np.sin(psi0))
            x[2]  = x[2] + self.ds * (np.cos(theta0))
            s += self.ds
            self.skeleton.append(skeletonElements(np.copy(x),theta,curvature,s,self.ds,self.radius,self.psi0,psiC0,psiG0,self.dt,self.boxsize,growth=self.growth,T_G_noise=self.T_G_noise,sig_G_noise=self.sig_G_noise))
        self.flatten()

    def set_markers(self,separation):
        for k in range(1,len(self.skeleton)):
            if np.mod(k,separation)==0:
                self.skeleton[k].set_marker()
                
        
    def updateSpatialPosition(self):
        s = 0
        for k in range(1,len(self.skeleton)):
            s += self.skeleton[k-1].ds
            self.skeleton[k].updateCurvilinearAbscissa(s)
            self.skeleton[k].updateOrientation(bt.angleDifference(self.skeleton[k-1].theta, - self.skeleton[k-1].curvature * self.skeleton[k-1].ds))
            self.skeleton[k].updateSpatialPosition(self.skeleton[k-1].x,self.skeleton[k-1].theta,self.skeleton[k-1].psiC,self.skeleton[k-1].ds)
        for skel in self.skeleton:
            skel.updateApicalAngle(self.skeleton[-1].theta)


    def addInteractions(self,name,intensity=0,direction = 0):
        if name in InteractionsName:
            self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
            for skel in self.skeleton:
                skel.addInteractions(name,intensity,direction)            
        else:
            print(' --- '+name+' is not part of the known interactions ')
            print(' --- please use one of the following interaction :')
            for names in InteractionsName:
                print(' --- --- '+str(names))

    def flatten(self):
        
        self.x=np.array([skel.x for skel in self.skeleton])
        self.s=np.array([skel.s for skel in self.skeleton])
        self.theta=np.array([skel.theta for skel in self.skeleton])
        self.curvature=np.array([skel.curvature for skel in self.skeleton])
        self.length = self.s[-1]

    def break_skeleton(self):
        self.skeleton_copy=[]        
        for skel in self.skeleton:
            skel.istoolong_update(self.ds_max)
            if skel.istoolong==0:
                self.skeleton_copy.append(skeletonElements(skel.x,skel.theta,skel.curvature,skel.s,skel.ds,skel.radius,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
                if skel.marker==1:
                    self.skeleton_copy[-1].set_marker()
            else:
                DS_temp=skel.ds/2
                self.skeleton_copy.append(skeletonElements(skel.x,skel.theta,skel.curvature,skel.s,DS_temp,skel.radius,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
                if skel.marker==1:
                    self.skeleton_copy[-1].set_marker()
                X_temp=np.copy(skel.x)
                X_temp[0]  = np.mod(X_temp[0] + DS_temp * (np.sin(skel.theta) * np.cos(skel.psi0)),skel.boxsize)
                X_temp[1]  = X_temp[1] + DS_temp * (np.sin(skel.theta) * np.sin(skel.psi0))
                X_temp[2]  = X_temp[2] + DS_temp * (np.cos(skel.theta))
                self.skeleton_copy.append(skeletonElements(np.copy(X_temp),skel.theta,skel.curvature,skel.s+DS_temp,DS_temp,skel.radius,skel.psi0,skel.psiC,skel.psiG,skel.dt,skel.boxsize,growth=skel.growth,T_G_noise=skel.T_G_noise,sig_G_noise=skel.sig_G_noise))
#                if skel.marker==1:
#                    self.skeleton_copy[-1].set_marker()
        self.skeleton=copy.deepcopy(self.skeleton_copy)
        self.skeleton_copy=[]
        self.updateSpatialPosition()
        self.flatten()

    def updateGrowth(self,time_counter):
        if np.mod(time_counter,self.T_G_noise)==1:
            self.G_noise= float(np.random.uniform(-self.sig_G_noise,self.sig_G_noise,1))

        if self.growth['name'] in GrowthMode:
            for skel in self.skeleton:
                
                if self.growth['name'] == 'Exponential':
                    skel.growth['growthRate'] = self.growth['growthRate']+self.G_noise 
                if self.growth['name'] == 'Apical':
                    
                    if (self.length-skel.s) < self.growth['growthZone']:

                        skel.growth['growthRate'] = self.growth['growthRate']+self.G_noise
                    else:
                        
                        skel.growth['growthRate'] = 0

    def update(self,time_counter):
#        self.break_skeleton()
        self.updateGrowth(time_counter)
        
        for skel in self.skeleton:
            skel.update(time_counter)
        self.updateSpatialPosition()

        self.flatten()

    def updateCollectiveInteraction(self):
        for skel in self.skeleton:
            skel.interactions = self.interactions

##############################################################################
##############################################################################
##############################################################################
##############################################################################
##############################################################################

class Roots:
    def __init__(self,radius,N,dx,T_P_noise,sig_P_noise,T_G_noise,sig_G_noise,theta0,ds_max=0.01,dt=0.01,growth = 'no',growthRate=1,growthZone = 1,elastic=0):
        self.time_counter=0
        self.radius=radius
        self.N = N
        self.dx = dx
        self.sig_P_noise=sig_P_noise
        self.T_P_noise=T_P_noise
        self.sig_G_noise=sig_G_noise
        self.T_G_noise=T_G_noise
        self.roots = []
        self.interactions = []
        self.collectiveInteractions = []
        self.collectiveInteractionsList = [] 
        self.P_Noise_vec=[]
        self.Noise_P_vec=[]
        self.growthZone=growthZone
        self.ds_max=ds_max
        self.elastic=elastic
#        self.G_Noise_vec=[]
        for k in range(0,N):
#            random_angle1=float(np.random.normal(0,20*pi/180,1))
            self.roots.append(Plant(x0 = [k*dx+0.5,0,0],theta0 = theta0,ds_max=ds_max,dt=dt,L=N*dx,radius=radius,growth = growth,growthRate=growthRate,growthZone = growthZone,sig_G_noise=self.sig_G_noise,T_G_noise=self.T_G_noise))

    def update(self):
        self.time_counter+=1
        if self.collectiveInteractions:
            for root in self.roots:
                root.break_skeleton()
            self.collectiveComputation()
            for k in range(0,self.N):
#                self.roots[k].break_skeleton()
                if type(self.intensityCollective)==type(np.ones(1)):
                    self.collectiveInteractionsList[k]['intensity'] = self.intensityCollective[k]
                    self.collectiveInteractionsList[k]['direction'] = self.directionCollective[k]
                    self.roots[k].updateCollectiveInteraction()
                elif type(self.intensityCollective)==type([1]):
                    qqq=0
                    for skel in self.roots[k].skeleton:
                        skel.interactions=[{'name': 'ApicalRelative', 'intensity': self.intensityCollective[k][qqq], 'direction': self.directionCollective[k][qqq]}]
                        qqq+=1
                self.roots[k].update(self.time_counter)
            if self.elastic:
                self.elastic_response()
        else:
            for root in self.roots:
                root.break_skeleton()
                root.update(self.time_counter)
            if self.elastic:
                self.elastic_response()
            
    def addInteractions(self,name,intensity=0,direction = 0):
        if name in InteractionsName:
            self.interactions.append({'name':name,'intensity':intensity,'direction':direction})
            for root in self.roots:
                root.addInteractions(name,intensity,direction)
        

        else:
            print(' --- '+name+' is not part of the known interactions ')
            print(' --- please use one of the following interaction :')
            for names in InteractionsName:
                print(' --- --- '+str(names))

    def addCollectiveInteraction(self,name,repulsionZone,attractionZone,repulsionIntensity,attractionIntensity):
        if name in CollectiveInteractionsName:
            self.collectiveInteractions.append({'name':name,'repulsionZone':repulsionZone,'attractionZone':attractionZone,'repulsionIntensity':repulsionIntensity,'attractionIntensity':attractionIntensity})
            for root in self.roots:
                root.addInteractions('ApicalRelative',0,0)
                self.collectiveInteractionsList.append(root.interactions[-1])
        else:
            print(' --- '+name+' is not part of the known collective interactions ')
            print(' --- please use one of the following collective interaction :')
            for names in CollectiveInteractionsName:
                print(' --- --- '+str(names))

    def collectiveComputation(self):
        self.flatten()
        for interaction in self.collectiveInteractions:
            if interaction['name'] == 'Apical':
                self.tipDistance()
                interactionTip =np.copy(self.distanceTip)
                interactionTip[(self.distanceTip>0) & (self.distanceTip<interaction['repulsionZone'])]=-1
                interactionTip[(self.distanceTip>interaction['repulsionZone']) & (self.distanceTip<interaction['attractionZone'])]=1
                interactionTip[self.distanceTip>interaction['attractionZone']]=0
                
                self.directionCollective = np.sum(interactionTip*self.alphaTip,0)/(self.N-1)
                
                self.intensityCollective = self.directionCollective*0
                self.intensityCollective[np.sum(np.abs(interactionTip),0)>0] =1.0
            
            if interaction['name'] == 'Apical_Apical_Morse':
                jnd=0
                self.directionCollective=np.zeros((self.N,))
                self.intensityCollective=np.zeros((self.N,))
                if np.mod(self.time_counter,self.T_P_noise)==1:
                    self.P_Noise_vec =np.random.normal(0,self.sig_P_noise, (np.size(self.directionCollective,0)))
                for root in self.roots:
                    R=[]
                    interactionTip=[]
                    X=root.skeleton[-1]
                    phasor=0.0
                    intensity=0.0
                    for root2 in self.roots:
                        if root2==root:
                            continue
                        R=bt.dist2d(X,root2.skeleton[-1],periodic=1,L=11)
                        angle=bt.angle2d(X,root2.skeleton[-1],periodic=1,L=11)
                        if interaction['attractionZone']==np.inf:
                            intensity=interaction['repulsionIntensity']-interaction['attractionIntensity']
                        else:
                            intensity=(interaction['repulsionIntensity']/interaction['repulsionZone'])*np.exp(-R/interaction['repulsionZone']) - (interaction['attractionIntensity']/interaction['attractionZone'])*np.exp(-R/interaction['attractionZone'])
                        phasor+=intensity*np.exp(1j*angle)
#                    print(np.angle(phasor))

                    self.directionCollective[jnd]=-bt.angleDifference(np.angle(phasor),X.theta)+self.P_Noise_vec[jnd]
                    self.intensityCollective[jnd]=np.abs(phasor)
                    jnd+=1
            
            if interaction['name'] == 'Apical_Global_Morse':
                jnd=0
                self.directionCollective=np.zeros((self.N,))
                self.intensityCollective=np.zeros((self.N,))
                if np.mod(self.time_counter,self.T_P_noise)==1:
                    self.P_Noise_vec =np.random.normal(0,self.sig_P_noise, (np.size(self.directionCollective,0)))
                for root in self.roots:
                    R=[]
                    interactionTip=[]
                    X=root.skeleton[-1]
                    phasor=0.0
                    intensity=0.0
                    for root2 in self.roots:
                        if root2==root:
                            continue
                        for skel in root.skeleton:
                            R=bt.dist2d(X,skel,periodic=1,L=11)
                            angle=bt.angle2d(X,skel,periodic=1,L=11)
                            if interaction['attractionZone']==np.inf:
                                intensity=interaction['repulsionIntensity']-interaction['attractionIntensity']
                            else:
                                intensity=(interaction['repulsionIntensity']/interaction['repulsionZone'])*np.exp(-R/interaction['repulsionZone']) - (interaction['attractionIntensity']/interaction['attractionZone'])*np.exp(-R/interaction['attractionZone'])
                            phasor+=intensity*np.exp(1j*angle)

                    self.directionCollective[jnd]=-bt.angleDifference(np.angle(phasor),X.theta)+self.P_Noise_vec[jnd]
                    self.intensityCollective[jnd]=np.abs(phasor)
                    jnd+=1

            if interaction['name'] == 'Global_Global_Morse':
                jnd=0
                self.directionCollective=[]
                self.intensityCollective=[]
                if np.mod(self.time_counter,self.T_P_noise)==1:
                    self.P_Noise_vec =np.random.normal(0,self.sig_P_noise, ((self.N,)))
                for root in self.roots:
                    self.directionCollective.append(np.zeros((np.size(root.skeleton),)))
                    self.intensityCollective.append(np.zeros((np.size(root.skeleton),)))                    
                    qq=0
                    for X in root.skeleton:
                        R=[]
                        interactionTip=[]
                        phasor=0.0
                        intensity=0.0
                        for root2 in self.roots:
                            if root2==root:
                                continue
                            for skel in root2.skeleton:
                                R=bt.dist2d(X,skel,periodic=1,L=11)
                                angle=bt.angle2d(X,skel,periodic=1,L=11)
                                if interaction['attractionZone']==np.inf:
                                    intensity=interaction['repulsionIntensity']-interaction['attractionIntensity']
                                else:
                                    intensity=(interaction['repulsionIntensity']/interaction['repulsionZone'])*np.exp(-R/interaction['repulsionZone']) - (interaction['attractionIntensity']/interaction['attractionZone'])*np.exp(-R/interaction['attractionZone'])
                                phasor+=intensity*np.exp(1j*angle)
                        self.directionCollective[jnd][qq]=-bt.angleDifference(np.angle(phasor),X.theta)+self.P_Noise_vec[jnd]
                        self.intensityCollective[jnd][qq]=np.abs(phasor)
                        qq+=1
                    jnd+=1

            if interaction['name'] == 'Global_Apical_Morse':
                jnd=0
                self.directionCollective=[]
                self.intensityCollective=[]
                if np.mod(self.time_counter,self.T_P_noise)==1:
                    self.P_Noise_vec =np.random.normal(0,self.sig_P_noise, ((self.N,)))
                for root in self.roots:
                    self.directionCollective.append(np.zeros((np.size(root.skeleton),)))
                    self.intensityCollective.append(np.zeros((np.size(root.skeleton),)))                    
                    qq=0
                    for X in root.skeleton:
                        R=[]
                        interactionTip=[]
                        phasor=0.0
                        intensity=0.0
                        for root2 in self.roots:
                            if root2==root:
                                continue
                            skel=root2.skeleton[-1]
                            R=bt.dist2d(X,skel,periodic=1,L=11)
                            angle=bt.angle2d(X,skel,periodic=1,L=11)
                            if interaction['attractionZone']==np.inf:
                                intensity=interaction['repulsionIntensity']-interaction['attractionIntensity']
                            else:
                                intensity=(interaction['repulsionIntensity']/interaction['repulsionZone'])*np.exp(-R/interaction['repulsionZone']) - (interaction['attractionIntensity']/interaction['attractionZone'])*np.exp(-R/interaction['attractionZone'])
                            phasor+=intensity*np.exp(1j*angle)
                        self.directionCollective[jnd][qq]=-bt.angleDifference(np.angle(phasor),X.theta)+self.P_Noise_vec[jnd]
                        self.intensityCollective[jnd][qq]=np.abs(phasor)
                        qq+=1
                    jnd+=1



    def flatten(self):
        
        self.xTip=np.array([root.x[-1] for root in self.roots])
        self.thetaTip=np.array([root.theta[-1] for root in self.roots])

    def tipDistance(self):
        self.distanceTip,self.alphaTip = bt.distPointToPoint(self.xTip,self.thetaTip,1,11)

    def elastic_response(self):
        FBR=0.5
        thresh = 5*pi/180 #angle threshold to avoid small forces
        extra_length=0.5*self.growthZone #extra length for bending "lignification"
        bool1=0
        qq=0
        qq2=0
        k=0
        while k<self.N:
            print(self.time_counter,k,qq2)
            for skel_gz in self.roots[k].skeleton:
                if skel_gz.s<=self.roots[k].length-self.growthZone:
                    continue
                X=copy.deepcopy(skel_gz)
                for k2 in range(0,self.N):
                    if 1:#k2!=k:
                        for skel in self.roots[k2].skeleton:
                            if k==k2:
                                if np.abs(skel.s-X.s)<0.2:
                                    continue
                            if bt.dist2d(X,skel,periodic=1,L=self.N*self.dx)<0.3:#self.ds_max*5:
                                (is_inter,angle_inter)=bt.intersection2D(X,skel,periodic=1,L=self.N*self.dx)
                                qq+=1
                                qq2+=is_inter
                                if qq2>100:
                                    break
                                
                                if is_inter>0:
                                    bool1=1
                                    if angle_inter>pi/2:
                                        angle_inter+=-pi
                                    if angle_inter<-pi/2:
                                        angle_inter+=pi
                                    for skel2 in self.roots[k].skeleton:
                                        if skel2.s<X.s and skel2.s > self.roots[k].length-self.growthZone-extra_length:
                                            if np.abs(angle_inter)>thresh:
                                                skel2.curvature+=-FBR*(skel2.s-(self.roots[k].length-self.growthZone-extra_length))*np.sin(angle_inter)
                                            else:
                                                skel2.curvature+=-FBR*(skel2.s-(self.roots[k].length-self.growthZone-extra_length))*np.sin(thresh)*np.sign(angle_inter)
                                    self.roots[k].updateSpatialPosition()
                                    self.roots[k].flatten()
                                    if skel.s>self.roots[k2].length-self.growthZone:
                                        for skel2 in self.roots[k2].skeleton:
                                            if skel2.s<skel.s and skel2.s > self.roots[k2].length-self.growthZone-extra_length:
                                                if np.abs(angle_inter)>thresh:
                                                    skel2.curvature+=FBR*(skel2.s-(self.roots[k2].length-self.growthZone-extra_length))*np.sin(angle_inter)
                                                else:
                                                    skel2.curvature+=FBR*(skel2.s-(self.roots[k2].length-self.growthZone-extra_length))*np.sin(thresh)*np.sign(angle_inter)
                                        self.roots[k2].updateSpatialPosition()
                                        self.roots[k2].flatten()                        
                                    break
                    if bool1:
                        break
                if bool1:
                    break
            if bool1:
                bool1=0
            else:
                k+=1
        for k in range(0,self.N):    
            self.roots[k].updateSpatialPosition()
            self.roots[k].flatten()                        
        self.flatten()