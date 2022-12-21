"""
  ReactionTransformOE.py: OEChem edition of AbstractReactionTransform
"""

from AbstractReactionTransform import AbstractReactionTransform
# from openeye.oechem import *

from BasicMol import BasicMol
from MolCollection import MolCollection
from SynthesizedMol import SynthesizedMol

class ReactionTransformOE(AbstractReactionTransform):
  """
    @param smirksString: SMIRKS string for the whole generic reaction
    @param name: string - reaction identifier
  """  
  def __init__(self,smirksString,name=""):
    AbstractReactionTransform.__init__(self,smirksString,name)
        
    #self.reactantMatchers= list()
    #for smarts in self.rctSmarts:
    #  self.reactantMatchers.append(OESubSearch(smarts))

    #Setup a list of numReactants empty MolCollections
    self.rctRgts= [None] * len(self.rctSmarts) 
    for i in xrange(len(self.rctRgts)):
      self.rctRgts[i]= MolCollection()
    
  # ---------------------------------------------------------------------------
  """ 
    Returns the number of times the mol passed in matches the ith generic reactant template
    If it returns 0, the mol cannot serve as the ith reactant
    Else if it returns 1, the mol can unambiguously can serve as the ith reactant (only 1 susceptible group)
    Else there is more than one reactive group in this mol that matches the reactive group specified for
      the ith reactant
   
    @param basicMol: a BasicMol object
  """ 
  def numMatchesAsReactant(self,ith,basicMol):
    #matchIter= self.reactantMatchers[ith].Match(basicMol.molObject)
    # OESubSearch.Match() returns an iterator, of type OEMatchBaseIter, not a list
    
    # Some changes to try to get around memory issues
    matcher= OESubSearch(self.rctSmarts[ith])
    basicMol.buildMolObject(True)
    matchIter= matcher.Match(basicMol.molObject,True)
    basicMol.clearMolObject()
    
    if (matchIter == None):
      return(0)
    else:
      numMatches= 0
      for match in matchIter:
        numMatches= numMatches + 1
    
    del(matcher)
    del(matchIter)
    
    return(numMatches)
  
  
  
#==============================================================================
#                  BELOW IS ALL TESTING CODE
#==============================================================================
""" Test relying on keyboard input only """
def keyboardTest():
   smirksString= raw_input( "Enter a SMIRKS string: " )
   tfm= ReactionTransformOE(smirksString)
 
   print("This reaction has " + str(tfm.getNumReactants()) + " reactants")

   moreInput= True
   while (moreInput):
     smilesString= raw_input( "Enter a molecule as a SMILES string (Return to quit): " )
     smilesString= smilesString.strip()
     if (len(smilesString) > 0):
       mol= BasicMol("",smilesString)
       mol.buildMolObject(true)
   
       for i in xrange(tfm.getNumReactants()):
         num= tfm.numMatchesAsReactant(i,mol)
         print("Molecule matched reactant " + str(i) + ", " + str(num) + " times")
     else:
       moreInput= False

#==============================================================================
""" subroutine of fileTest(), below. """
def getLinesFile(fName,firstToken=False,doStrip=False):
  try:
    inFile= file(fName,'r')
  except Exception:
    print("Error opening file: " + fName)
    raise

  try:
    inLines= inFile.readlines()
    if (firstToken):
      for i in xrange(len(inLines)):
        tokens= inLines[i].split()
        inLines[i]= tokens[0]
    elif (doStrip):
      for i in xrange(len(inLines)):
        inLines[i]= inLines[i].strip()
      
  finally:
    inFile.close()
  
  return(inLines)
  
#==============================================================================
""" Test using file input of molecules and reactions, and file output of matching
    results """
def fileTest():
    fName= raw_input( "Enter infile with SMIRKS string: ")
    inLines= getLinesFile(fName,True)
    tfm= ReactionTransformOE(inLines[0])
 
    print("This reaction has " + str(tfm.getNumReactants()) + " reactants")

    fName= raw_input( "Enter infile with SMILES strings: ")
    inLines= getLinesFile(fName,True)

    outFileNames= [None] * tfm.getNumReactants()
    outFiles= [None] * tfm.getNumReactants()
    for i in xrange(len(outFileNames)):
      outFileNames[i]= "matched" + str(i) + ".smi"
      outFiles[i]= file(outFileNames[i],'w')
    noMatchFile= file("nomatch.smi",'w')

    for smilesString in inLines:
      mol= BasicMol("",smilesString,None)
      mol.buildMolObject(True)
   
      someMatch= False
      for i in xrange(tfm.getNumReactants()):
        num= tfm.numMatchesAsReactant(i,mol)
        if (num > 0):
          someMatch= True
          outFiles[i].write(smilesString + "  " + str(num) + "\n")
     
      if (not someMatch):
        noMatchFile.write(smilesString + "\n")



    print("Molecules not matching anything written to: nomatch.smi")
    print("Molecules matching I-th reactant written to matchedI.smi (e.g. matched0.smi)")

#==============================================================================
""" Test using file input of molecules and reactions, and file output of matching
    results.  This version also makes the products as a combi lib - not suitable
    for large input files.
"""
def fileTest2():
    fName= raw_input( "Enter infile with SMIRKS string: ")
    inLines= getLinesFile(fName,True)
    rxn= ReactionTransformOE(inLines[0])
 
    print("This reaction has " + str(rxn.getNumReactants()) + " reactants")

    fName= raw_input( "Enter infile with SMILES strings: ")
    inLines= getLinesFile(fName,True)

    outFileNames= [None] * rxn.getNumReactants()
    outFiles= [None] * rxn.getNumReactants()
    for i in xrange(len(outFileNames)):
      outFileNames[i]= "matched" + str(i) + ".smi"
      outFiles[i]= file(outFileNames[i],'w')
    noMatchFile= file("nomatch.smi",'w')
    prodFile= file("products.smi",'w')

    for smilesString in inLines:
      mol= BasicMol("",smilesString,None)
      mol.buildMolObject(True)
   
      someMatch= False
      for i in xrange(rxn.getNumReactants()):
        num= rxn.numMatchesAsReactant(i,mol)
        if (num > 0):
          someMatch= True
          outFiles[i].write(smilesString + "  " + str(num) + "\n")
     
      if (not someMatch):
        noMatchFile.write(smilesString + "\n")

    # Generate all possible products
    
    # First identify reagents
    for smilesString in inLines:
      mol= SynthesizedMol("",smilesString,None)
      mol.buildMolObject(True)
   
      for i in xrange(rxn.getNumReactants()):
        numMatches= rxn.numMatchesAsReactant(i,mol)
        if (numMatches == 1):                # Only use case of unambiguous match
          rxn.rctRgts[i].add(mol)            # For the rxn, save all instances of matches for all reactants

    # Setup the combi library
    libgen= OELibraryGen(rxn.smirksString)
    libgen.SetExplicitHydrogens(True)
    libgen.SetValenceCorrection(True)
    for i in xrange(rxn.getNumReactants()):
      for j in xrange(len(rxn.rctRgts[i])):
        basicMol= rxn.rctRgts[i][j]
        basicMol.buildMolObject(True)
        libgen.AddStartingMaterial(basicMol.molObject,i)

    # Make the products
    for product in libgen.GetProducts():
      smiStrg= OECreateCanSmiString(product)
      prodFile.write(smiStrg + "\n")


    print("Molecules not matching anything written to: nomatch.smi")   
    print("Molecules matching I-th reactant written to matchedI.smi (e.g. matched0.smi)")

#==============================================================================
# Test Driver code
if __name__ == "__main__":
   fileTest2()
