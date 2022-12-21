"""
  AbstractReactionTransform.py: provides facilities for a generic reaction transform,
  i.e., one in which reagent/reaction patterns are decribed, but specific reagents 
  are not. SMARTS patterns are used to describe the reagents that will match as reactants.
"""

class AbstractReactionTransform:
  """
    @param smirksString: SMIRKS string for the whole generic reaction
    @param name: string - reaction identifier
  """
  def __init__(self,smirksString,name=""):
    self.smirksString= smirksString  # 
    self.name= name                  # 
    self.rctSmarts= None             # List with 1 SMARTS string for each reactant
    self.rctRgts= None               # List of MolCollections, 1 for each reactant
                                     #   rctRgts[i][j] is the jth specific molecule that can play the role of the ith generic reactant

    pos= smirksString.find(">>")
    if (pos < 0):
      raise Exception("Invalid SMIRKS string in ReactionTransform.__init__(): " + smirksString)
    
    leftSide= smirksString[0:pos] 
    self.rctSmarts= leftSide.split('.')
  
  # ---------------------------------------------------------------------------
  """ Returns the number of reactants needed to fullfill this ReactionTransform """
  def getNumReactants(self):
    return(len(self.rctSmarts))
    
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
    raise NotImplementedError("numMatchesAtReactant must be defined in subclasses of AbstractReactionTransform!")

  # ---------------------------------------------------------------------------
  """
    Returns True if the molecule passed in can serve as one of the reactants for this reaction
    Multiple matches are allowed (i.e., multiple groups in a molecule that can be mapped to the
    susceptible group for one reaction component.
   
    @param basicMol: a BasicMol object
    @param allowMultiMatches: boolean, if True, then a reactant may match more than once for a given reaction,
      i.e., a compound with 3 Br's will be be matched for an SN1 reaction.  If False, then there must be exactly 
      one matching group in the molecule for a reaction term for it to match
  """
  def canBeAReactant(self,basicMol,allowMultiMatches):
    for i in xrange(self.getNumReactants()):
      numMatches= self.numMatchesAsReactant(i,basicMol)
      if (allowMultiMatches):
        if (numMatches > 0):
          return(True)
      elif (numMatches == 1):
        return(True)

    return(False)

#==============================================================================

# Test Driver code
if __name__ == "__main__":
   

  fName= raw_input("Enter SMIRKS reaction file name: ")
  inFile= open(fName,'r')
  smirksLines= inFile.readlines()
  for line in smirksLines:
    if (line[0] == '#'):
      continue
    line= line.strip()
    if (len(line) > 0):
      tokens= line.split()
      smirks= tokens[0]
      if (len(tokens) > 1):
        name= tokens[1]
      else:
        name= ""
      rxn= AbstractReactionTransform(smirks,name)
      
      print("For rxn: " + line)
      print("Components are: ")
      for smarts in rxn.rctSmarts:
        print("  " + smarts)

  inFile.close()
 