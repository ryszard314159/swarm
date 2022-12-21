import os.path
 
# from openeye.oechem import *

from BasicMol import BasicMol
from AbstractMolScoreFunc import AbstractMolScoreFunc

"""
  LingoScoreFunc.py: Simple 2D scoring function using the OpenEye OELingoSim class,
  which implements a 2D similarity metric.
  
  NOTE: Current done using the molObject, but there are options to just use SMILES
  (canonical smiles) to construct the OELingoSim object and obtain the similarity, 
  and this may be faster - algorithm seems to use SMILES directly.  Think that through
  before implementing though - which SMILES (isomeric, canonical...) matters, 
  
"""
class LingoScoreFunc(AbstractMolScoreFunc):
  # ---------------------------------------------------------------------------------------
  """
    @param templateSDF_FileName: filename of an SDF file that holds the
      3D coordinates of the ligand or substrate, etc that we want similarity to
  """
  def __init__(self, templateMolFName):
    self.templateMol= None
    self.comparator= None
    
    if (not (os.path.isfile(templateMolFName))):
      raise Exception("Error, template molecule file: " + templateMolFName + " not found in LingoScoreFunc.__init__()")
      
    inStream= oemolistream(templateMolFName)
    mol= OEGraphMol()    
    resultCode= OEReadMolecule(inStream,mol)
    if (resultCode == 1):
      smiStrg= OECreateCanSmiString(mol)
      self.templateMol= BasicMol("template",smiStrg,mol)
      #self.templateMol.assignCanonicalSmiles(True)         # This also adds hydrogens
      
      self.templateMol.buildMolObject(True)
    else:
      raise Exception("Error parsing template molecule in file: " + templateMolFName)

    # Based on canonical SMILES string
    #self.comparator= OELingoSim(smiStrg)
    

    self.comparator= OELingoSim(self.templateMol.molObject)
  # ---------------------------------------------------------------------------------------
  """
    Implementation of the standard value method, here corresponding to the Lingo score
      @param basicMol: a BasicMol to be evaluated, molObject is ignored
      
    Time scale here is (roughly) 0.2 milliseconds per call - fast!
  """
  def value(self,basicMol):
    if (basicMol is None or basicMol.smiStrg is None):
      return(0.0)
    else: 
      # Based on canonical SMILES string
      # return(self.comparator.Similarity(basicMol.smiStrg))
    
      # Based on molObject
      basicMol.buildMolObject(True)
      score= self.comparator.Similarity(basicMol.molObject)
      basicMol.clearMolObject()
      return(score)
  
