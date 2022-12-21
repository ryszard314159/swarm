"""
  IDSequence.py: simple class to generate a sequence of identifiers for
  molecules with a fixed string prefix and an increasing integer suffix.
  
  TODO: Add option for index fields with left padded zeroes out to some
  fixed width.
"""
class IDSequence(object):
  # ---------------------------------------------------------------------------
  """
    @param prefix: a string desired to be the prefix for this sequence
  """
  def __init__(self,prefix):
    self.prefix= prefix
    self.index= 1

  # ---------------------------------------------------------------------------
  """
    Allows you to start the integer suffix of the sequence at a value different
    than the default (1), allowing you to, e.g. continue an earlier sequence.

    @param index: a positive integer
  """
  def setStartingIndex(self,index):
    self.index= index
    
  # ---------------------------------------------------------------------------
  """ Generates the next ID string from this sequence """
  def getNextID(self):
    currValue= self.index
    self.index= self.index + 1
    return(self.prefix + str(currValue))
  
  #----------------------------------------------------------------------------
  """
    Assuming a list of BasicMols which have previously generated names from an
    IDSequence, this loops over all mols and returns the highest value of the 
    integer index found, as an int.
  
    This allows a previously used sequence to be continued. 
  """
  def getHighestIndex(self,prefix,molList):
    start= len(prefix)
    maxInd= -1
    for basicMol in molList:
      if (basicMol.name.find(prefix) == 0):
        ind= int(basicMol.name[start:len(basicMol.name)])
        if (ind > maxInd):
          maxInd= ind

    return(maxInd)
