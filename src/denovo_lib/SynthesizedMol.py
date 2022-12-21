""" 
  SynthesizedMol.py - extends ScoredMol with a Synthetic Origin object to allow
  tracing back how the molecule was made in-silico.
"""

import sys

from ScoredMol import ScoredMol
from SyntheticOrigin import SyntheticOrigin

class SynthesizedMol(ScoredMol):
  """
    @param name: a string for the name(i.e., title, id) for the molecule
    @param smiStrg: a Daylight SMILES string
    @param molObject: a molecule Object from a molecule class library (e.g. OEChem)
    @param prodScore: a float reflecting the goodness of this molecule as a product,
                        normally in the range [0..1]
    @param rgtScore: a float reflecting the goodness of this molecule when used as
                        a reagent, normally in the range [0..1]  
    @param synOrigin: a SyntheticOrigin object describing how this molecule was 
                        "made" in-silico
  """
  def __init__(self,name,smiStrg,molObject,prodScore=0.0,rgtScore=0.5,synOrigin=None):
    ScoredMol.__init__(self,name,smiStrg,molObject,prodScore,rgtScore)
    self.synOrigin= synOrigin
    self.isTerminal= False         # flag to indicate if this molecule should be expanded on
                                   # termination would normally be due to complexity max-ing out
    
    self.myRxns= list()            # myRxns is actually a list of 2-elem lists.
                                   # Each sublist has elem[0] as the rxn, and elem[1] as the
                                   # integer index of the reactant position assigned to this mol


  # ---------------------------------------------------------------------------
  """
    Full printout of this SynthesizedMol, plus the entire synthetic pathway leading up to it.
    Transforms, reactants (and their origins) are written on separate lines and organized using
    indentation. Tree is traversed to enable printout using indirect recursion
  """
  def printFull(self,outFile,indentation=""):
    # First print the molecule rep
    outFile.write(indentation + self.__str__() + "\n")

    # Now print it's synthetic origin; note increase in indentation level
    if (self.synOrigin != None):
      self.synOrigin.printFull(outFile,indentation + "\t")

  # ---------------------------------------------------------------------------
  """ 
    Returns the total number of synthetic steps used to generate this molecule, 
    through recursion.
  """ 
  def numSyntheticSteps(self):
    numSteps= 0
    if (self.synOrigin != None):
      numSteps= 1
      for i in xrange(len(self.synOrigin.reactants)):
        numSteps= numSteps + self.synOrigin.reactants[i].numSyntheticSteps()
      
    return(numSteps)
      

#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
   # These SMIRKS reactions aren't validated, btw, just placeholders here
   bromo_ethane= SynthesizedMol("","CCBr",None)

   ethanol= SynthesizedMol("","CCO",None)
   ethanol.synOrigin= SyntheticOrigin([bromo_ethane],"[C:1H3]Br>>[C:1H3]O")

   acetic_acid= SynthesizedMol("","CC(=O)O",None)
   acetic_acid.synOrigin= SyntheticOrigin([ethanol],"[C:1H3][O:2H]>>[C:1H](=[O:2])O")

   ammonia= SynthesizedMol("","NH3",None)

   acetamide= SynthesizedMol("","CC(=O)N",None)
   acetamide.synOrigin= SyntheticOrigin([acetic_acid,ammonia],"[C:1H](=[O:2])O.[NH3]>>[C:1H](=[O:2])N")

   acetamide.printFull(sys.stdout)
   print("Acetamide was made in: " + str(acetamide.numSyntheticSteps()) + " steps")
   