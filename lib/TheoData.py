import numpy as np

class TheoData(object):
    """Stores temporary data required during Theobald RMSD calculation.
    
    Notes:
    Storing temporary data allows us to avoid re-calculating the G-Values
    repeatedly. Also avoids re-centering the coordinates."""
    
    def __init__(self, XYZData, NumAtoms=None, G=None):
        """Create a container for intermediate values during RMSD Calculation.
        
        Notes:
        1.  We remove center of mass.
        2.  We pre-calculate matrix magnitudes (ConfG)"""
        
        if NumAtoms is None or G is None:
            NumConfs=len(XYZData)
            NumAtoms=XYZData.shape[1]
        
            self.centerConformations(XYZData)
        
            NumAtomsWithPadding=4+NumAtoms-NumAtoms%4
        
            # Load data and generators into aligned arrays
            XYZData2 = np.zeros((NumConfs, 3, NumAtomsWithPadding), dtype=np.float32)
            for i in range(NumConfs):
                XYZData2[i,0:3,0:NumAtoms] = XYZData[i].transpose()
            
            #Precalculate matrix magnitudes
            ConfG = np.zeros((NumConfs,),dtype=np.float32)
            for i in xrange(NumConfs):
                ConfG[i] = self.calcGvalue(XYZData[i,:,:])
            
            self.XYZData=XYZData2
            self.G=ConfG
            self.NumAtoms=NumAtoms
            self.NumAtomsWithPadding=NumAtomsWithPadding
            self.CheckCentered()
        else:
            self.XYZData = XYZData
            self.G = G
            self.NumAtoms = NumAtoms
            self.NumAtomsWithPadding = XYZData.shape[2]
    
    def CheckCentered(self, Epsilon=1E-5):
        """Raise an exception if XYZAtomMajor has nonnzero center of mass(CM)."""
        
        XYZ=self.XYZData.transpose(0,2,1)
        x=np.array([max(abs(XYZ[i].mean(0))) for i in xrange(len(XYZ))]).max()
        if x>Epsilon:
            raise Exception("The coordinate data does not appear to have been centered correctly.")
    
    @staticmethod
    def centerConformations(XYZList):
        """Remove the center of mass from conformations.  Inplace to minimize mem. use."""
        
        for ci in xrange(XYZList.shape[0]):
            X=XYZList[ci].astype('float64')#To improve the accuracy of RMSD, it can help to do certain calculations in double precision.  This is _not_ one of those operations IMHO.
            X-=X.mean(0)
            XYZList[ci]=X.astype('float32')
        return
    
    @staticmethod
    def calcGvalue(XYZ):
        """Calculate the sum of squares of the key matrix G.  A necessary component of Theobold RMSD algorithm."""
        
        conf=XYZ.astype('float64')#Doing this operation in double significantly improves numerical precision of RMSD
        G = 0
        G += np.dot(conf[:,0],conf[:,0])
        G += np.dot(conf[:,1],conf[:,1])
        G += np.dot(conf[:,2],conf[:,2])
        return G

    def __getitem__(self, key):
        # to keep the dimensions right, we make everything a slice
        if isinstance(key, int):
            key = slice(key, key+1)
        return TheoData(self.XYZData[key], NumAtoms=self.NumAtoms, G=self.G[key])
    
    def __setitem__(self, key, value):
        self.XYZData[key] = value.XYZData
        self.G[key] = value.Gdef

    def __len__(self):
        return len(self.XYZData)