""" 
  ScoredMol.py - extends BasicMol with fitness scoring to facilitate use in 
  molecule ranking, selection, denovo design.
"""

from BasicMol import BasicMol


class ScoredMol(BasicMol):
  """ 
    Retain a class average reagent score - this is used as a substitution (dynamic default)
    for molecules which haven't been used yet - they are assumed to look like everything else.
    Reference these with class prefix, e.g., ScoredMol.avgRgtScore 
  """
  avgRgtScore= 0.5
  numRgtScores= 1


  """
    @param name: a string for the name(i.e., title, id) for the molecule
    @param smiStrg: a Daylight SMILES string
    @param molObject: a molecule Object from a molecule class library (e.g. OEChem)
    @param prodScore: a float reflecting the goodness of this molecule as a product,
                        normally in the range [0..1]
    @param rgtScore: a float reflecting the goodness of this molecule when used as
                        a reagent, normally in the range [0..1]  
  """
  def __init__(self,name,smiStrg,molObject,prodScore=0.0,rgtScore=0.0):
    BasicMol.__init__(self,name,smiStrg,molObject)
    self.prodScore= prodScore
    self.rgtScore= rgtScore
    self.numTimesRgt= 0                 # num of times it's been used as a reagent

  # ---------------------------------------------------------------------------
  """ 
    Use this method to update the reagent score - don't update the instance variable directly.
    @param derivedProdScore: a prodScore assigned to a product synthesized using this molecule
    
    @returns: nothing
  """
  def updateReagentScore(self,derivedProdScore):
    self.rgtScore= ((self.numTimesRgt * self.rgtScore) + derivedProdScore)/(self.numTimesRgt + 1)
    self.numTimesRgt= self.numTimesRgt + 1

    # Update class average vars
    ScoredMol.avgRgtScore= ((ScoredMol.numRgtScores * ScoredMol.avgRgtScore) + derivedProdScore)/(ScoredMol.numRgtScores + 1)
    ScoredMol.numRgtScores= ScoredMol.numRgtScores + 1

  # ---------------------------------------------------------------------------
  """ 
    Use this method to fetch the reagent scpre update the reagent score - don't reference
    the instance variable directly..

    @returns: the averaged reagent score for this molecule, or class average if it hasn't
              been used as a reagent yet.
  """
  def getReagentScore(self):
    if (self.numTimesRgt == 0):
      return(ScoredMol.avgRgtScore)
    else:
      return(self.rgtScore)

  # ---------------------------------------------------------------------------
  """ 
    Override str function to enable e.g. basic printout. 

    @returns: concatenation of the SMILES string, name, product score and reagent score,
              tab delimited
  """
  def __str__(self):
    return(self.smiStrg + "\t" + self.name + "\t" + "%5.4f" % self.prodScore + "\t" + "%5.4f" % self.rgtScore)


  # **************************************************************************
  # *                    ALL METHODS BELOW ARE STATIC METHODS                
  # **************************************************************************

  # ---------------------------------------------------------------------------
  """ Static compare function for sorting by product score """
  def compareOnProdScore(scoredMol1,scoredMol2):
    if (scoredMol1.prodScore > scoredMol2.prodScore):
      return(1)
    elif (scoredMol1.prodScore < scoredMol2.prodScore):
      return(-1)
    else:
      return(0)

  compareOnProdScore= staticmethod(compareOnProdScore)

  # ---------------------------------------------------------------------------
  """ Static compare function for sorting by reagent score """
  def compareOnRgtScore(scoredMol1,scoredMol2):
    if (scoredMol1.getRgtScore() > scoredMol2.getRgtScore()):
      return(1)
    elif (scoredMol1.getRgtScore() < scoredMol2.getRgtScore()):
      return(-1)
    else:
      return(0)

  compareOnRgtScore= staticmethod(compareOnRgtScore)

  # ---------------------------------------------------------------------------
  """ Static compare function for sorting by SMILES string """
  def compareOnSmiStrg(scoredMol1,scoredMol2):
    if (scoredMol1.smiStrg > scoredMol2.smiStrg):
      return(1)
    elif (scoredMol1.smiStrg < scoredMol2.smiStrg):
      return(-1)
    else:
      return(0)
    
  compareOnSmiStrg= staticmethod(compareOnSmiStrg)


