"""
  Simple class to identify timeout errors, which might often have recovery strategies
"""
class TimeoutException(Exception):
  # ---------------------------------------------------------------------------
  def __init__(self,strg):
    if (strg == None):
      self.strg= ""
    else:
      self.strg= strg
  
  # ---------------------------------------------------------------------------
  def __str__(self):
    return(self.strg)

  # ---------------------------------------------------------------------------
  def __repr__(self):
    return("TimeoutException: " + self.strg)
  