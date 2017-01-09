
# coding: utf-8

# # Ligands Framework 

# In[ ]:

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit import DataStructs
from rdkit.Chem import MACCSkeys


# In[3]:

class MoleculeHelper(object):
    "Molecule helper class to create fingerprints and compute similarity"    
    
    __smiles = None
    __fingerprint = None
    
    def __init__(self, smiles):
        self.__smiles = smiles
    
    # Override this method to implement different fingerprint algorithms
    def fingerprint(self, molecule):
        "Returns Morgan's fingerprint"
        return AllChem.GetMorganFingerprintAsBitVect(molecule, 2)
    
    # Override this method to implement different similarity algorithms
    def similarityAlgorithm(self, otherFingerprint, metric=DataStructs.TanimotoSimilarity):
        "Returns similarity coefficient between two molecule fingerprints. Defaults to Tanimoto."        
        return DataStructs.FingerprintSimilarity(self.__getFingerprint(), otherFingerprint, metric)        
    
    def getSmiles(self):
        "Returns the SMILES representation of this molecule."
        return self.__smiles
    
    def __getFingerprint(self):
        if self.__fingerprint == None:
            try:
                m = Chem.MolFromSmiles(self.__smiles)
                if m != None:
                    self.__fingerprint = self.fingerprint(m)
            except Exception as ex:
                print("**** SPARK - fingerprint: %s" % self.__smiles)
                print(ex)
                self.__fingerprint = None
                #raise
        return self.__fingerprint
            
    def getFingerprint(self):
        "Returns the fingerprint represenatation of this molecule."
        return self.__getFingerprint()
    
    def similarity(self, otherMoleculeHelper, metric=DataStructs.TanimotoSimilarity):
        sim = -1 # error
        try:
            otherFingerprint = otherMoleculeHelper.getFingerprint()
            if self.__getFingerprint() != None and otherFingerprint != None:
                sim = self.similarityAlgorithm(otherFingerprint, metric)
        except Exception as ex:
            print("**** SPARK - similarity: %s" % self.__smiles)
            print(ex)
            #raise
        return sim
        
    
    def similarityFromSmiles(self, otherSmiles, metric=DataStructs.TanimotoSimilarity):
        return self.similarity(type(self)(otherSmiles), metric)


# In[5]:

# class that inherits MoleculeHelper and implements MACCS fingerprinting
class MoleculeMACCSHelper(MoleculeHelper):    
    
    def fingerprint(self, molecule):
        return MACCSkeys.GenMACCSKeys(molecule)    

