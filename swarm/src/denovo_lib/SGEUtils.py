import os,sys,time
from subprocess import *

from Timer import Timer
from TimeoutException import TimeoutException

"""
  SGEUtils.py: Some static utility methods for working with the Sun Grid Engine.  
"""
class SGEUtils:

  # ---------------------------------------------------------------------------
  """
    Returns the number of processors available in a given queue

    @param queueName: name of the queue to be examined.
  """
  def getAvailableProcCount(self,queueName="blade.q"):
    cmdLine= "qstat -f -q " + queueName + " | grep 'BIP'"
    osProcess= Popen(cmdLine, shell=True, stdout=PIPE)
    stdoutFile= osProcess.stdout

    lineNum= 1
    inLine= ""
    foundCount= 0
    while (inLine is not None and lineNum < 500):
      inLine= stdoutFile.readline()
      if (inLine is None):
        break
      else:
        inLine= inLine.strip()
        if (len(inLine) == 0):
          break
      
      tokens= inLine.split()
      cpusToken= tokens[2]           # Result will look like 1/4, meaning 1 of 4 cores in use"
      cpusToken= cpusToken.split("/")
      filled= int(cpusToken[0])
      total= int(cpusToken[1])
      avail= total - filled
      foundCount= foundCount + avail

      lineNum= lineNum + 1
      
    stdoutFile.close()
    
    return(foundCount)
  

  # ---------------------------------------------------------------------------
  """
    Polls the Sun Grid Engine for the presence of any given jobs from a user,
    using the qstat utility.  Testing on my Linux box indicates a single call 
    takes roughly 1/10th second.
    
    @param userName: string - login name for the user to look for.  If None,
                     the current process username is used
    @param jobNames: set of strings corresponding to the job names to look for.
                     Can be omitted to select by userName only.
    @param queueName: name of the queue to be examined.

    @returns: an int - the number of such jobs found 
  """
  def pollQueue(self,userName=None,jobNames=None,queueName="blade.q"):
    if (userName is None):
      userName= str(os.getlogin())
      
    cmdLine= "qstat -u " + userName + " -q " + queueName 
    cmdTokens= cmdLine.split()
    osProcess= Popen(cmdLine, shell=True, stdout=PIPE)
    stdoutFile= osProcess.stdout
    
    lineNum= 1
    inLine= ""
    foundCount= 0
    while (inLine is not None and lineNum < 500):
      inLine= stdoutFile.readline()
      if (inLine is None):
        break
      else:
        inLine= inLine.strip()
        if (len(inLine) == 0):
          break
      
      if (lineNum >= 3):
        tokens= inLine.split()
        qJobName= tokens[2]
        qUserName= tokens[3]
        if (jobNames is None):
          if (qUserName == userName):
            foundCount= foundCount + 1
        else:
          if (qJobName in jobNames and qUserName == userName):
            foundCount= foundCount + 1

      lineNum= lineNum + 1
      
    stdoutFile.close()
    
    return(foundCount)
    
  # ---------------------------------------------------------------------------
  """
    Polls the Sun Grid Engine for the presence of any given jobs from a user,
    using the qstat utility, until there are no more remaining.
    
    @param userName: string - login name for the user to look for.  If None,
                     the current process username is used
    @param jobNames: set of strings corresponding to the job names to look for.
                     Can be omitted to select by userName only.
    @param queueName: name of the queue to be examined.
    @param maxWait: float - maximum time to wait for completion of the jobs, in seconds,
                    default is 14 hours
    @param sleepTime: float - sleep time between polling calls, in seconds
  """
  def pollQueueUntilDone(self,userName=None,jobNames=None,queueName="blade.q",maxWait=(14.0*3600.0),sleepTime=3.0):
    waitTimer= Timer()
    
    foundCount= 1
    sleepTime= 3.0  # seconds between polling
    print("Count of qstat jobs found is: ")
    while (waitTimer.getElapsed() < maxWait and foundCount > 0):
      foundCount= self.pollQueue(userName,jobNames,queueName)
      print(" " + str(foundCount)) # ,

      if (foundCount == 0):
        break
      else:
        time.sleep(sleepTime)
    
    if (foundCount > 0):
      raise TimeoutException("Exceeded wait time of " + str(maxWait) + " seconds for Sun Grid Engine in SGEUtils.pollQueueUntilDone()")
    else:
      print(" SGE jobs finished in " + str(waitTimer.getElapsed()) + " seconds")


#==============================================================================
# Test Driver code
if __name__ == "__main__":
  app= SGEUtils()
  queueName= "blade.q"

  avail= app.getAvailableProcCount(queueName)
  print("Found " + str(avail) + " processors available on the SGE in the " + queueName + " queue")

  app.pollQueueUntilDone(None,None,queueName,30.0)  
  
