""" 
  BasicMol.py - base class, simple container class for a SMILES string and a 
  proprietary molecule object from a molecule class library (e.g. OEChem)
"""
# from openeye.oechem import OEGraphMol,OEParseSmiles,OEAddExplicitHydrogens,OECreateCanSmiString

class BasicMol(object):

  """
    @param name: a string for the name(i.e., title, id) for the molecule
    @param smiStrg: a Daylight SMILES string
    @param molObject: a molecule Object from a molecule class library (e.g. OEChem)
  """
  def __init__(self,name,smiStrg,molObject):
    if (name == None):
      self.name= ""
    else:
      self.name= name
      
    self.smiStrg= smiStrg
    self.molObject= molObject
  
  # ---------------------------------------------------------------------------
  """
    Returns a hash int calculated from the SMILES string, (must use canonical 
    SMILES for this to work as designed). This is overridden to enable membership
    tests in sets and dictionaries.
  """
  def __hash__(self):
    return(self.smiStrg.__hash__())
    
  # ---------------------------------------------------------------------------
  """
    Returns true if the SMILES strings of the 2 BasicMols are equal (must use 
    canonical SMILES). This is overridden to enable membership tests in sets and 
    dictionaries.

    @param other: another BasicMol
  """
  def __eq__(self,other):
    return(self.smiStrg.__eq__(other.smiStrg))

  # ---------------------------------------------------------------------------
  """ 
    Override str function to enable e.g. basic printout

    @returns: concatenation of the SMILES string, and name, tab delimited
  """
  def __str__(self):
    return(self.smiStrg + "\t" + self.name)

  # ---------------------------------------------------------------------------
  """ Simple single line printout using the __str__ method  """
  def printBrief(self,outFile):
    outFile.write(self.__str__() + "\n")

  # ---------------------------------------------------------------------------
  """
    Builds the proprietary molObject - this method must be changed to use a different
    mol software library.
    
    @param addExplicitHydrogens: boolean, if true, explicit hydrogens are added
      to the molObject constructed.
  """
  def buildMolObject(self,addExplicitHydrogens):

    mol= OEGraphMol()
    parsedOK= OEParseSmiles(mol,self.smiStrg)
    if (parsedOK == 0):
      raise StandardError("SMILES string parse error for: " + self.smiStrg + " in BasicMol.buildMolObject()")

    mol.SetTitle(self.name)  
    if (addExplicitHydrogens):
      OEAddExplicitHydrogens(mol)

    self.molObject= mol
    
  # ---------------------------------------------------------------------------
  """
    Deletes the proprietary molObject - e.g., to conserve memory
  """
  def clearMolObject(self):

    self.molObject.Clear()
    del(self.molObject)

  # ---------------------------------------------------------------------------
  """
    Returns the number of graphs in this mol's SMILES string
  """
  def numGraphs(self):
    if (self.smiStrg == None):
      return(0)
    else:
      nFound= 0
      for ch in self.smiString:
        if (ch == '.'):
          nFound= nFound + 1

      return(nFound+1)

  # ---------------------------------------------------------------------------
  """
    Replace's this mol's SMILES string with the substring corresponding to 
    the largest subgraph
  """
  def chooseLargestGraph(self,addExplicitHydrogens):
    if (self.smiStrg != None):
      tokens= self.smiStrg.split(".")
      maxHeavy= 0
      biggestSmiStrg= None
      mol= OEGraphMol()
      
      for token in tokens:
        try:
          parsedOK= OEParseSmiles(mol,token)
          if (parsedOK):
            numHeavy= 0
            for atom in mol.GetAtoms():
              atomicSymbol= OEGetAtomicSymbol(atom.GetAtomicNum())
              if (not atomicSymbol == 'H'):
                numHeavy= numHeavy + 1

            if (numHeavy > maxHeavy):
              maxHeavy= numHeavy
              biggestSmiStrg= token
        except:                       # Just ignore for now...
          pass
          
      mol.Clear()
      del(molObject)
        
      if (biggestSmiStrg != None):
        self.smiStrg= biggestSmiStrg
        self.assignCanonicalSmiles(addExplicitHydrogens)

  # ---------------------------------------------------------------------------
  """
    Replaces the existing SMILES string with one created by a canonicalization
    algorithm.
    
    @param addExplicitHydrogens: boolean, if true, explicit hydrogens are added
      to the molObject constructed.
  """
  def assignCanonicalSmiles(self,addExplicitHydrogens):
    self.buildMolObject(addExplicitHydrogens)
    self.smiStrg= OECreateCanSmiString(self.molObject)
    self.clearMolObject()
    
#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
   # Test the hashing function for set membership:
   molSet= set()
   molSet.add(BasicMol("Acetic Acid","CC(=O)O",None))
   molSet.add(BasicMol("Ethyl Acetate","CC(=O)OCC",None))
   molSet.add(BasicMol("Propyl Amine","CCCN",None))
   molSet.add(BasicMol("2-Methyl Pyridine","c1cc(C)ncc1",None))
   molSet.add(BasicMol("Methyl Ethyl Ketone","CC(=O)CC",None))
   molSet.add(BasicMol("Propanoyl Chloride","CCC(=O)Cl",None))

   bMol= BasicMol("Alias for Propyl Amine","CCCN",None)
   if (bMol in molSet):
     print("Already have bMol")
   else:
     print("bMol is new")       
   
   moreInput= True
   while (moreInput):
     inString= raw_input("Enter SMILES string: ")
     inString= inString.strip()
     if (len(inString) > 0):
       aMol= BasicMol("Some Mol",inString,None)
       if (aMol in molSet):
         print("Already have this mol")
       else:
         print("This mol is new")       
     else:
       moreInput= False

