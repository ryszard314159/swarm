import sys

# from openeye.oechem import *
from AbstractMolFilter import AbstractMolFilter

"""
  ReactiveFilter.py: Filtering based on susceptibility to a set of reactions - compounds that
  are reactive, pass.
"""
class ReactiveFilter(AbstractMolFilter):

  """
    @param rxnList: list of AbstractReactionTransform objects
    @param allowMultiMatches: boolean, if True, then a reactant may match more than once for a given reaction,
      i.e., a compound with 3 Br's will be be matched for an SN1 reaction.  If False, then there must be exactly 
      one matching group in the molecule for a reaction term for it to match
  """
  def __init__(self, rxnList,allowMultiMatches):
    self.rxnList= rxnList
    self.allowMultiMatches= allowMultiMatches

  # ---------------------------------------------------------------------------
  """
    Returns true if the molecule passed in passes all encoded criteria
    
    @param basicMol: a BasicMol to be evaluated, molObject member must be an OEGraphMol
  """
  def passes(self, basicMol):
  
    for rxn in self.rxnList:
      if rxn.canBeAReactant(basicMol,self.allowMultiMatches):
        return(True)

    #print "UNREACTIVE: " + basicMol.smiStrg
    return(False)
    

#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
   amide_formation_smirks= "[H]O[C:1]([#6:4])=[O:2].[H][N+:3]([#6:5])=[H]>>[H][N:3]([#6:5])[C:1]([#6:4])=[O:2]"
   
   filter= ReactiveFilter()
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
