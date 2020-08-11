#Packages imported
import numpy as np
import sympy as sp
import random as rnd

#random line pair generator
def rpoint():
    return [rnd.uniform(-10,10), rnd.uniform(-10,10), rnd.uniform(-10,10)]

def rseg():
    return [rpoint(), rpoint()]

def randpair():
    return[rseg(),rseg()]

print(randpair())

# triple product function
def trp (x,y,z):
    return np.dot(x,np.cross(y,z))

# denominator functions
def cyc (x,y,z):
    return np.dot(x,y)*np.linalg.norm(z)

#General AT function
def atan(a,b,alpha,d):

    if alpha == 0 or alpha == np.pi or d == 0:
        return 0

    else:
        t = a*b*np.sin(alpha) + d**2*(1/np.tan(alpha))
        m = d*np.sqrt(a**2 + b**2 - (2*a*b*np.cos(alpha))+d**2)

    return np.arctan(t/m)

#Construct invariants of two line segments from endpoints
class Segquad:

    def __init__(self,parameter):

        self.lonestart = np.array(parameter[0][0]).astype(float)
        self.loneend = np.array(parameter[0][1]).astype(float)
        self.ltwostart = np.array(parameter[1][0]).astype(float)
        self.ltwoend = np.array(parameter[1][1]).astype(float)

        #Parameterisation vctors
        self.lonevec = self.loneend - self.lonestart
        self.ltwovec = self.ltwoend - self.ltwostart

        #Parameterised line equations for segments
        self.loneeq = self.lonestart + t*self.lonevec
        self.ltwoeq = self.ltwostart + s*self.ltwovec

        #Length and cross product for segments
        self.cross = np.cross(self.lonevec, self.ltwovec)
        self.seg1len = np.linalg.norm( self.lonevec )
        self.seg2len = np.linalg.norm( self.ltwovec )

        #symbolic form of gauss integral function
        self.difference = (self.lonestart + t * self.lonevec) - (self.ltwostart + s * self.ltwovec)
        self.distance = sp.sqrt(self.difference[0]**2 + self.difference[1]**2 + self.difference[2]**2)
        self.gaussfunc = np.dot( np.cross( self.lonevec, self.ltwovec ), self.difference) / self.distance**3

        #Angle between segments
        self.angle = np.arccos( np.dot(self.lonevec,self.ltwovec) / (self.seg1len * self.seg2len) )
        #print( 'angle = ' + repr( self.angle ) )
        self.alpha = self.angle/2




    #Calculates the signed distance of a line pair
    def segdist(self):

        return trp(self.lonevec,self.ltwovec,self.startdist)/np.linalg.norm(self.cross)

    #Calculate the a1, a2 co-ordinates of a line pair
    def acoord(self):

        b = self.startdist/((np.sin(self.angle))**2)
        x = self.lonevec/self.seg1len
        y = self.ltwovec/self.seg2len

        a1 = np.dot(((y*np.cos(self.angle))-x),b)
        a2 = np.dot(((y - x*np.cos(self.angle))),b)

        return [a1,a2]

    #Calculates the lk funtion from four separate AT functions
    def segatan(self):

        a1 = self.acoord()[0]
        a2 = self.acoord()[1]
        b1 = self.acoord()[0] + self.seg1len
        b2 = self.acoord()[1] + self.seg2len
        alpha = self.angle
        d = self.segdist()

        return (atan(a1,b2,d,alpha) + atan(a2,b1,d,alpha) - atan(a1,a2,d,alpha) - atan(b1,b2,d,alpha))/(4*np.pi)


#Linking Number for the Hopf Link as two orthogonal triangles
hsimp_1 = [[[-1,0,1],[-1,0,-1]],[[-1,0,-1],[1,0,0]],[[1,0,0],[-1,0,1]]]
hsimp_2 = [[[0,0,0],[2,1,0]],[[2,1,0],[2,-1,0]],[[2,-1,0],[0,0,0]]]

simplink = 0

for i in hsimp_1:
    for j in hsimp_2:

        segp = Segquad([i,j])

        simplink += segp.segatan()
        print(simplink)

print(simplink)


#Calculation of linking number between a segment of the line x = 0, y = 0
#and a  1-periodic lattice of  orthogonal line segments in the plane z = 0
#between y = +/- 1
def per1segs(x):

    base_seg = [[0,0,-1],[0,0,1]]
    orth_segs = [[[1,-1,0],[1,1,0]],[[-1,1,0],[-1,-1,0]]]
    rightquad = Segquad([base_seg,orth_segs[0]])
    leftquad = Segquad([base_seg,orth_segs[1]])
    perlink = rightquad.segatan() + leftquad.segatan()

    chain = []

    for i in range (1,x):
        for j in orth_segs:
            shift = [4*i,0,0]
            pstart = [a + b for a,b in zip(j[0],shift)]
            pend = [c + d for c,d in zip (j[1],shift)]

            chain.append([pstart,pend])

            mstart = [a - b for a,b in zip(j[0],shift)]
            mend = [c - d for c,d in zip (j[1],shift)]

            chain.append([mstart,mend])


    for j in chain:
        segp = Segquad([base_seg,j])

        perlink += segp.segatan()

    return perlink
