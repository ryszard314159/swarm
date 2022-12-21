""" 
  SyntheticOrigin.py: This describes how a molecule was made - the reactants and the 
  reactive transform (generic chemical reaction) that was employed.
""" 
class SyntheticOrigin:
  """ 
    @param reactants: is a list of SynthesizedMol objects
    @param smirksString: is a Daylight SMIRKS string
    @param rxnName: string identifier for the reaction
  """ 
  def __init__(self,reactants,smirksString,rxnName=""):
    self.reactants= reactants
    self.smirksString= smirksString
    self.rxnName= rxnName

  # ---------------------------------------------------------------------------
  """ 
    Simple single line printout of the SMIRKS string
    @param file: is a Python file object
  """ 
  def printBrief(self,file):
    file.write(self.smirksString + " " + rxnName + "\n")

  # ---------------------------------------------------------------------------
  """ 
    Full printout of this SyntheticOrigin, plus that for each reactant serving as input.
    Transforms, reactants (and their origins) are written on separate lines and organized using
    indentation. Tree is traversed to enable printout using indirect recursion
    @param file: is a Python file object
    @param indentation: is a string composed of blanks that gets lengthened on recursion to
                        provide the correct level of indentation on an output line
  """ 
  def printFull(self,outFile,indentation=""):
    outFile.write(indentation + self.smirksString + " " + self.rxnName + "\n")

    for reactant in self.reactants:
      reactant.printFull(outFile,indentation + "\t")


