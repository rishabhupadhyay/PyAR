import numpy as np
from math import *
class Molecule:
    def __init__(self, atomlist, coords):
        self.atomlist=atomlist
        self.coords=coords

    def find_centre(self):
        coords=self.coords
        nat = len(self.coords)
        xcentre = sum([coords[i][0] for i in range(nat)])/nat
        ycentre = sum([coords[i][1] for i in range(nat)])/nat
        zcentre = sum([coords[i][2] for i in range(nat)])/nat
        centre = xcentre, ycentre, zcentre
        return centre

    def find_centre_of_mass(self):
        coords=self.coords
        atomlist=self.atomlist
        ZAtom = {"H":2.0, "B": 10.0, "C":12.0, "O":16.0, "N":14.0, "S":32.0}
        nat = len(coords)
        mass=sum([ZAtom[atomlist[i]] for i in range(nat)])
        xcentre = sum([ZAtom[atomlist[i]]*coords[i][0] for i in range(nat)])/mass
        ycentre = sum([ZAtom[atomlist[i]]*coords[i][1] for i in range(nat)])/mass
        zcentre = sum([ZAtom[atomlist[i]]*coords[i][2] for i in range(nat)])/mass
        com = xcentre, ycentre, zcentre
        return com

    def translate(self, magnitude):
        coords=self.coords
        atomlist=self.atomlist
        nat=len(coords)
        translated = [(atomlist[i], coords[i][0] + magnitude[0], coords[i][1] + magnitude[1], coords[i][2] + magnitude[2])
                  for i in range(nat)]
        return translated

    
    def translate_to_centre(self):
        coords=self.coords
        centre = self.find_centre()
        centered = self.translate(centre)
        return centered
    
    def get_distance(self,a,b):
        dx = a[0] - b[0]
        dy = a[1] - b[1]
        dz = a[2] - b[2]
        r = sqrt(dx**2+dy**2+dz**2)
        return r

    def getAverageRadii(self):
        coords=self.coords
        natoms = len(coords)
        x, y, z = self.find_centre()
        return sum([self.get_distance([x,y,z],coords[i]) for i in range(natoms)])/natoms

    def standardDeviation(self):
        coords=self.coords
        natoms=len(coords)
        Av=self.getAverageRadii()
        x,y,z = self.find_centre()
        SD=0
        for i in range(natoms):
        	SD=SD+pow(self.get_distance([x,y,z],coords[i]),2)
        SD=sqrt(SD/natoms)
        return SD

    def tensor(self):
        coords=self.coords
        ZAtom = {"H":2.0, "B": 10.0, "C":12.0, "O":16.0, "N":14.0, "S":32.0}
        natoms=len(coords)
        moi=np.zeros((3,3)) 
        for i in range(natoms):
            Z = ZAtom[self.atomlist[i]]
            x = coords[i][0]
            y = coords[i][1]
            z = coords[i][2]

            moi[0][0]=moi[0][0]+Z*(y*y+z*z)/float(natoms)
            moi[0][1]=moi[0][1]-Z*(x*y)/float(natoms)
            moi[0][2]=moi[0][2]-Z*(x*z)/float(natoms)
    
            moi[1][0]=moi[1][0]-Z*(x*y)/float(natoms)           
            moi[1][1]=moi[1][1]+Z*(x*x+z*z)/float(natoms)
            moi[1][2]=moi[1][2]-Z*(y*z)/float(natoms)
    
            moi[2][0]=moi[2][0]-Z*(x*z)/float(natoms)
            moi[2][1]=moi[2][1]-Z*(y*z)/float(natoms)
            moi[2][2]=moi[2][2]+Z*(x*x+y*y)/float(natoms)

        return moi
'''
    def get_inertial_coords(self):
        eigenvalues, eigenvectors = la.eig(moi)
        geom = [(coords_com[i][1],coords_com[i][2],coords_com[i][3]) for i in range(len(coords_com))]
        geom = np.array(np.dot(geom, eigenvectors))
        coords_com=[(coords_com[i][0],geom[i][0],geom[i][1],geom[i][2])for i in range(len(geom))]
        # coords_com = [(geom[i][1],geom[i][2],geom[i][3]) for i in range(len(geom))]
        moi = tensor(coords_com)
        order = [0, 1, 2]
        for p in range(3):
            for q in range(p+1, 3):
                if (moi.item(p, p) < moi.item(q, q)):
                    temp = order[p]
                    order[p] = order[q]
                    order[q] = temp
        moveaxes = np.zeros((3, 3))
        for p in range(3):
            moveaxes[p][order[p]] = 1.0
        geom = np.dot(geom, moveaxes)
        coords_com=[(coords_com[i][0],geom[i][0],geom[i][1],geom[i][2])for i in range(len(geom))]
        return coords_com
'''
molecule1=Molecule(['O','H','H'],[[4,1,3],[1,2,3],[3,4,5]])
print molecule1.find_centre()
#print molecule1.translate([1,1,1])
#print molecule1.translate_to_centre()
#print molecule1.getAverageRadii()
#print molecule1.standardDeviation()
#print molecule1.tensor()
print molecule1.find_centre_of_mass()
