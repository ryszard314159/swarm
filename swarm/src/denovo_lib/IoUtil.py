from SynthesizedMol import SynthesizedMol

import os,sys

"""
  IoUtil.py: Class with some static convenience methods for input and output.

  Methods including logging to process-global error and trace files.
  Because these are treated as static class variables, and the calls are all
  static, any line in any module can simply call, e.g., IoUtil.logErrorString("Whoops!")
  without any conditional testing.  Initialization of these files earlier
  in the process using the open methods determines if anything will actually
  happen when the call is made.
"""
class IoUtil(object):
  errorFile= None
  traceFile= None
  
  #----------------------------------------------------------------------------
  """
    Opens the process-global trace file with the desired name, for writing
    (static method).
  """
  def openTraceFile(fName):
    IoUtil.traceFile= file(fName,'w')
    
  openTraceFile= staticmethod(openTraceFile)

  #----------------------------------------------------------------------------
  """
    Opens the process-global error file with the desired name, for writing
    (static method).
  """
  def openErrorFile(fName):
    IoUtil.errorFile= file(fName,'w')
    
  openErrorFile= staticmethod(openErrorFile)

  #----------------------------------------------------------------------------
  """ 
    Writes the desired string to the process-global trace file - no newline
    is used, so append to argument string as needed (static method).
    
    @param echoToTerm: boolean, if True, the string is also echoed to the terminal
  """
  def logTraceString(aString,echoToTerm=False):
    if (IoUtil.traceFile != None):
      IoUtil.traceFile.write(aString)
      IoUtil.traceFile.flush()

    if (echoToTerm):
      print(aString)

  logTraceString= staticmethod(logTraceString)


  #----------------------------------------------------------------------------
  """ 
    Writes the desired string to the process-global error file - no newline
    is used, so append to argument string as needed (static method).

    @param echoToTerm: boolean, if True, the string is also echoed to the terminal
  """
  def logErrorString(aString,echoToTerm=False):
    if (IoUtil.errorFile != None):
      IoUtil.errorFile.write(aString)
      IoUtil.errorFile.flush()

    if (echoToTerm):
      print(aString)
      
  logErrorString= staticmethod(logErrorString)

  #----------------------------------------------------------------------------
  """ 
    Writes the desired molecule to the process-global trace file, in brief form.
    (static method)
    
    @param basicMol: a BasicMol
  """
  def logTraceMol(basicMol):
    if (IoUtil.traceFile != None):
      basicMol.printBrief(IoUtil.traceFile)
      IoUtil.traceFile.flush()

  logTraceMol= staticmethod(logTraceMol)

  #----------------------------------------------------------------------------
  """ 
    Writes the desired molecule to the process-global trace file, in full form
    if it's a SynthesizedMol, else in brief form.  (static method)

    @param basicMol: a BasicMol
  """
  def logErrorMol(basicMol):
    if (IoUtil.errorFile != None):
      if (isinstance(basicMol,SynthesizedMol)):
        basicMol.printFull(IoUtil.errorFile)
      else:
        basicMol.printBrief(IoUtil.errorFile)
      IoUtil.errorFile.flush()

  logErrorMol= staticmethod(logErrorMol)

  #----------------------------------------------------------------------------
  """ Closes the process-global trace file (static method). """
  def closeTraceFile():
    if (IoUtil.traceFile != None):
      IoUtil.traceFile.close()
    
  closeTraceFile= staticmethod(closeTraceFile)

  #----------------------------------------------------------------------------
  """ Closes the process-global error file (static method). """
  def closeErrorFile():
    if (IoUtil.errorFile != None):
      IoUtil.errorFile.close()
    
  closeErrorFile= staticmethod(closeErrorFile)

  #----------------------------------------------------------------------------
  """ 
    If the file name passed in doesn't have a path already, adds an explicit path
    to it - the current working directory. (static method)
    
    @returns: string - the processed filename
. """
  def prependCWD(fName):
    if (fName.find("/") == -1 and fName.find("\\") == -1):
      return(os.getcwd() + os.sep + fName)
    else:
      return(fName)
    
  prependCWD= staticmethod(prependCWD)

  #----------------------------------------------------------------------------
  """ 
    Looks for a file in the Python interpreter's current search path.  If found,
    returns the complete path to this file. (static method)
    
    @param fName: string - name of the desired file, without any path prefix

    @returns: if found, string - the file name with full path to it prepended,
      else None
. """
  def findFileInPyPath(fName):
    for dir in sys.path:
      pathFName = os.path.join(dir,fName)
      if (os.path.isfile(pathFName)):
        return pathFName
    
    return None

  findFileInPyPath= staticmethod(findFileInPyPath)
 