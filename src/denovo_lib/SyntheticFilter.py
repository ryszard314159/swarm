from BasicFilter import BasicFilter

"""
  SyntheticFilter.py: Extends BasicFilter to test also the number of synthetic steps 
  used to make a SynthesizedMol - the type it operates on.  
"""
class SyntheticFilter(BasicFilter):

  # ---------------------------------------------------------------------------
  """
    @param minCarbons: minimum number of carbons required
    @param minHeavyAtoms: minimum number of heavy atoms required
    @param maxHeavyAtoms: maximum number of heavy atoms allowed
    @param maxRotors: maximum number of bond rotors allowed
    @param maxChiralAtoms: maximum number of chiral centers allowed
    @param allowedAtoms: presence of any atoms outside this set causes failure, (set of Atomic symbol strings)
    @param forbiddenGroups: presence of any of the functional groups in this list causes failure (list of SMARTS strings)  
    @param maxSyntheticSteps: maximum allowed synthetic steps leading to the molecule)  
  """
  def __init__(self, minCarbons=0, minHeavyAtoms=1, maxHeavyAtoms=30, maxRotors=10, maxChiralAtoms=6, 
               allowedAtoms=set(['C','H','N','O','P','S','F','Cl','Br','I','Na','K','Ca']), 
               forbiddenGroups=None, maxSyntheticSteps=20):
    BasicFilter.__init__(self, minCarbons, minHeavyAtoms, maxHeavyAtoms, maxRotors, maxChiralAtoms, 
               allowedAtoms,forbiddenGroups)
    self.maxSyntheticSteps= maxSyntheticSteps

  # ---------------------------------------------------------------------------
  """
    Returns true if the molecule passed in passes all encoded criteria
    @param synMol: a SynthesizedMol to be evaluated, molObject member must be an OEGraphMol
  """
  def passes(self, synMol):
    return(BasicFilter.passes(self,synMol) and synMol.numSyntheticSteps() < self.maxSyntheticSteps)
    