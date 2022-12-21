import copy

"""
  CombiIter.py: A combinatorial iterator class, used for generating all index 
  combinations for a multidimensional combinatorial problem.  See the test code
  at bottom for hints on usage.
"""
class CombiIter:
  """
    @param dimsList: list of integer dimensions - each value is terminal values 
       of the indices at each position to be varied - they will each range from 
       0 to (terminus-1)
  """
  def __init__(self,dimsList):
    self.value= [0] * len(dimsList)
    self.dimsList= dimsList
    self.returnedMax= False
    
  # ---------------------------------------------------------------------------
  """
    For internal use, returns true if the current configuration is the maximum value
    allowed
  """
  def atMaxValue(self):
    for i in xrange(len(self.value)):
      if (self.value[i] < (self.dimsList[i] -1)):
        return(False)
    
    return(True)
      
  # ---------------------------------------------------------------------------
  """
    For internal use, increments the current configuration, recursive.
    
    @param pos: integer indicating the position being varied
  """
  def incr(self,pos):
    if (self.value[pos] < (self.dimsList[pos]) - 1):
      self.value[pos]= self.value[pos] + 1
    else:
      self.value[pos]= 0
      self.incr(pos+1)
      
  # ---------------------------------------------------------------------------
  """
    Returns the next configuration of indices - a list of integers, If no more,
    then returns None
  """
  def next(self):
     if (self.returnedMax):
       return(None)
     else:
       retValue= copy.deepcopy(self.value)
       if (self.atMaxValue()):
         self.returnedMax=  True
       else:
         self.incr(0)
       
     return(retValue)
     

#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
   #dimsList= [3,2,4]   
   dimsList= [2,2]   
   iter= CombiIter(dimsList)
   
   more= True
   numCombos= 0
   while (more):
     combo= iter.next()
     if (combo == None):
       more= False
     else:
       for coord in combo:
         print(str(coord) + " ")
       print("")
       numCombos= numCombos + 1

   print("")
   print("Total Combinations: " + str(numCombos))
  