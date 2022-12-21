import sys

# from openeye.oechem import *
from AbstractMolFilter import AbstractMolFilter

"""
  BasicFilter.py: Simple filtering on number of carbons, rotors, heavy atoms, 
  membership of all heavy atoms in an allowed set, absence of a set of forbidden 
  functional groups.  Operates on BasicMol objects.
"""
class BasicFilter(AbstractMolFilter):

  """
    @param minCarbons: minimum number of carbons required
    @param minHeavyAtoms: minimum number of heavy atoms required
    @param maxHeavyAtoms: maximum number of heavy atoms allowed
    @param maxRotors: maximum number of bond rotors allowed
    @param maxChiralAtoms: maximum number of chiral centers allowed
    @param allowedAtoms: presence of any atoms outside this set causes failure, (set of Atomic symbol strings)
    @param forbiddenGroups: presence of any of the functional groups in this list causes failure (list of SMARTS strings)  
  """
  def __init__(self, minCarbons=0, minHeavyAtoms=1, maxHeavyAtoms=30, maxRotors=10, maxChiralAtoms=6, 
               allowedAtoms=set(['C','H','N','O','P','S','F','Cl','Br','I','Na','K','Ca']), forbiddenGroups=None):
    self.minCarbons= minCarbons
    self.minHeavyAtoms= minHeavyAtoms
    self.maxHeavyAtoms= maxHeavyAtoms
    self.maxRotors= maxHeavyAtoms
    self.maxChiralAtoms= maxChiralAtoms
    if (allowedAtoms == None):
      self.allowedAtoms= list()
    else:
      self.allowedAtoms= allowedAtoms
  
    if (forbiddenGroups == None):
      self.forbiddenGroups= list()
      self.forbiddenMatchers= list()
    else:
      self.forbiddenGroups= forbiddenGroups
      self.forbiddenMatchers= list()
      for smarts in self.forbiddenGroups:
        self.forbiddenMatchers.append(OESubSearch(smarts))

  # ---------------------------------------------------------------------------
  """ Returns a simple set of allowed atoms for an organic molecule/salt """
  def getDefaultAllowedAtoms():
    return(set(['C','H','N','O','P','S','F','Cl','Br','I','Na','K','Ca']))

  # ---------------------------------------------------------------------------
  """
    Returns true if the molecule passed in passes all encoded criteria
    @param basicMol: a BasicMol to be evaluated, molObject member must be an OEGraphMol
  """
  def passes(self, basicMol):
    basicMol.buildMolObject(True)
    mol= basicMol.molObject
  
    numCarbons= 0
    numForbiddenAtoms= 0
    numHeavyAtoms= 0
    numChiralAtoms= 0
    atomicSymbol= 'H'

    OEPerceiveChiral(mol) 
    for atom in mol.GetAtoms():
      atomicSymbol= OEGetAtomicSymbol(atom.GetAtomicNum())
      #print "Atomic Symbol: ", atomicSymbol
      if (atomicSymbol == 'C'):
        numCarbons= numCarbons + 1
      if (not atomicSymbol == 'H'):
        numHeavyAtoms= numHeavyAtoms + 1
      if (not atomicSymbol in self.allowedAtoms):
        numForbiddenAtoms= numForbiddenAtoms + 1
      if (atom.IsChiral()):
        numChiralAtoms= numChiralAtoms + 1
 
    numRotors = 0
    for bond in mol.GetBonds():
      if bond.IsRotor(): 
        numRotors = numRotors + 1

    noForbiddenGroups= True
    for matcher in self.forbiddenMatchers:
      if (matcher.SingleMatch(mol)):
        noForbiddenGroups= False
        break

    #print "numCarbons: ", numCarbons
    #print "numDisallowed: ", numDisallowed
    #print "numRotors: ", numRotors
    #print "numHeavyAtoms: ", numHeavyAtoms
    #print "numChiralAtoms: ", numChiralAtoms
    basicMol.clearMolObject()

    return (numCarbons >= 1 and numForbiddenAtoms == 0 and noForbiddenGroups and numRotors <= self.maxRotors and numHeavyAtoms >= self.minHeavyAtoms and numHeavyAtoms <= self.maxHeavyAtoms and numChiralAtoms < self.maxChiralAtoms)
      
 
#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
   filter= BasicFilter()
   print("Enter a SMILES string")
   smilesString= sys.stdin.readline() 
   smilesString= smilesString.strip()
 
   if (filter.passes(smilesString)):
     print("Molecule passes filter")
   else:
     print("Molecule fails filter")


'''
~kzth541/apps/openeye/bin/filter 
/usr/contrib/openeye/data/filter_lead.txt
'''
