import numpy as np 
import math




def distPointToPoint(x,theta,periodic,L):
    distance = np.zeros((len(x[:,0]),len(x[:,0])))
    X = []
    for k in range(0,len(x[0,:])):
        X1 = np.tile(x[:,k],(len(x[:,k]),1))
        X2 = np.tile(x[:,k],(len(x[:,k]),1)).T # transpose of X1
        if k==0 and periodic==1:  # first column = x coordinates. PBC are only along x
            dis1=np.minimum(np.abs(X2-X1),np.abs(L-np.abs(X2-X1))) #PBC check along x only
            distance += dis1**2
            M1=np.sign(X2-X1)
            O1=np.ones((np.size(M1,0),np.size(M1,1)))
            O1[np.where(dis1!=np.abs(X2-X1))]=-1 # where PBC have been applied
            M1=np.multiply(M1,O1)
            dis2=np.minimum(np.abs(X2-X1),np.abs(L-np.abs(X2-X1)))*M1
            X.append(dis2)

        else:
            distance += (X2-X1)**2
            X.append((X2-X1))    
                
    X1 = X[0]*np.cos(-theta)+X[2]*np.sin(-theta)
    X2 = -X[0]*np.sin(-theta)+X[2]*np.cos(-theta)
    angle = np.arctan2(X1,X2)
    distance = np.sqrt(distance)

    np.fill_diagonal(distance,0)
    np.fill_diagonal(angle,0)

    return distance, angle



def angleDifference(A1, A2):
    A = A1 - A2
    A = (A + math.pi) % (2 * math.pi) - math.pi # % = Modulus - remainder of the division of left operand by the right
    return A



def intersection2D(A,B,periodic,L):
    Ax1=np.copy(A.x[0])
    Ax2=np.copy(A.x[0]+A.ds*np.sin(A.theta))
    Bx1=np.copy(B.x[0])
    Bx2=np.copy(B.x[0]+B.ds*np.sin(B.theta))
    if periodic:
        if (Ax1-Ax2)>(L-1):
           Ax1+=-L
        if (Ax2-Ax1)>(L-1):
           Ax2+=-L
        if (Bx1-Bx2)>(L-1):
           Bx1+=-L
        if (Bx2-Bx1)>(L-1):
           Bx2+=-L
    if np.abs(A.theta-B.theta)<1e-15:
        return 0,0.0
    elif np.abs(Ax1-Ax2)<1e-15:
        x=Ax1
        CB = B.x[2]-Bx1/np.tan(B.theta)
        z=CB+x/np.tan(B.theta)
        za2=A.x[2]+A.ds*np.cos(A.theta)
        if (x < np.max((np.min((Ax1,Ax2)),np.min((Bx1,Bx2)))) or x > np.min((np.max((Ax1,Ax2)),np.max((Bx1,Bx2)))) or (z<np.min((za2,A.x[2]))) or (z>np.max((za2,A.x[2])))):
            return 0,0.0
        else:
            return 1,angleDifference(A.theta,B.theta)
    elif np.abs(Bx1-Bx2)<1e-15:
        x=Bx1
        CA = A.x[2]-Ax1/np.tan(A.theta)
        z=CA+x/np.tan(A.theta)
        zb2=B.x[2]+B.ds*np.cos(B.theta)
        if (x < np.max((np.min((Ax1,Ax2)),np.min((Bx1,Bx2)))) or x > np.min((np.max((Ax1,Ax2)),np.max((Bx1,Bx2)))) or (z<np.min((zb2,B.x[2]))) or (z>np.max((zb2,B.x[2])))):
            return 0,0.0
        else:
            return 1,angleDifference(A.theta,B.theta)     
    else:
        CA = A.x[2]-Ax1/np.tan(A.theta)
        CB = B.x[2]-Bx1/np.tan(B.theta)
        x=(CB-CA)/(1/np.tan(A.theta) - 1/np.tan(B.theta))
        za=CA+x/np.tan(A.theta)
        zb=CB+x/np.tan(B.theta)
        if (x < np.max((np.min((Ax1,Ax2)),np.min((Bx1,Bx2)))) or x > np.min((np.max((Ax1,Ax2)),np.max((Bx1,Bx2)))) or np.abs(za-zb)>1e-14):
            return 0,0.0
        else:
            return 1,angleDifference(A.theta,B.theta)
        
        
def intersection2D_r(A,B,periodic,L,radius):
        D=dist2d(A,B,periodic,L)
        if D>2*radius:
            return 0,0.0
        else:
            return 1,angleDifference(A.theta,B.theta)


def dist2d(A,B,periodic,L):
    if periodic:
        dis1=np.abs(A.x[0]-B.x[0])
        dis2=np.minimum(dis1,np.abs(L-dis1))
        dis3=np.abs(A.x[1]-B.x[1])
        dis4=np.minimum(dis3,np.abs(L-dis3))
        return np.sqrt(dis2**2 + dis3**2 +(A.x[2]-B.x[2])**2)
    else:
        return np.sqrt((A.x[0]-B.x[0])**2 + (A.x[1]-B.x[1])**2 +(A.x[2]-B.x[2])**2)


def angle2d(A,B,periodic,L):
    if periodic:
        dis1=np.abs(A.x[0]-B.x[0])
        if dis1>np.abs(L-dis1):
            DX=-np.sign(B.x[0]-A.x[0])*np.abs(L-dis1)
        else:
            DX=B.x[0]-A.x[0]
    else:
        DX=B.x[0]-A.x[0]
    DZ=B.x[2]-A.x[2]
    return np.arctan2(DX,DZ)

def angle3d(A,B,periodic,L):
    if periodic:
        dis1=np.abs(A.x[0]-B.x[0])
        if dis1>np.abs(L-dis1):
            DX=-np.sign(B.x[0]-A.x[0])*np.abs(L-dis1)
        else:
            DX=B.x[0]-A.x[0]
        dis2=np.abs(A.x[1]-B.x[1])
        if dis2>np.abs(L-dis2):
            DY=-np.sign(B.x[1]-A.x[1])*np.abs(L-dis2)
        else:
            DY=B.x[1]-A.x[1]
    else:
        DX=B.x[0]-A.x[0]
        DY=B.x[1]-A.x[1]
    DZ=B.x[2]-A.x[2]
    THETA = np.arctan2(np.sqrt(DY**2 + DX**2),DZ)
    PHI = np.arctan2(DY,DX)
    return PHI,THETA
#np.arccos(np.cos(THETA)*np.cos(A.theta)+np.sin(THETA)*np.sin(A.theta)*np.cos(PHI-A.psiC))

        
def writeCsvRoots(X,name,path,writeMode = 0):
    if writeMode == 0:
        fd = open(path+name,'wb')
    elif writeMode == 1:
        fd = open(path+name,'ab')
    
    wri = X



    np.savetxt(fd,np.c_[wri],delimiter=',',fmt= '%5.5f')
    fd.close()
