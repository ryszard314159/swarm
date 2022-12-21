import time

"""
  Timer.py: Simple convenience class for measuring elapsed time between events
"""
class Timer(object):
  # ---------------------------------------------------------------------------
  """ 
    Construction initializes this timer's initialTime
  """
  def __init__(self):
    self.initialTime= time.time()
      
  # ---------------------------------------------------------------------------
  """ 
    Resets this timer's initialTime to the present
  """
  def reset(self):
    self.initialTime= time.time()

  # ---------------------------------------------------------------------------
  """ 
    Returns the time elapsed since either this timer's construction, or it's
    most recent reset, whichever came last, in seconds, as a float.  Returns
    a larger value with each call, since time passes.
  """
  def getElapsed(self):
    currentTime= time.time()
    return(currentTime - self.initialTime)
