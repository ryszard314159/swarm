import os,sys,subprocess,time, os.path
from subprocess import *

from AbstractMolScoreFunc import AbstractMolScoreFunc
from IoUtil import IoUtil
from MolCollection import MolCollection
from SGEUtils import SGEUtils
from Timer import Timer
from TimeoutException import TimeoutException

"""
  ROCSScoreFunc.py: Mol scoring function leveraging the OpenEye ROCS shape
  matching application
"""
class ROCSScoreFunc(AbstractMolScoreFunc):
  """
    @param templateFName: filename of an SDF file that holds the
      3D coordinates of the ligand or substrate, etc that we want similarity to
  """
  def __init__(self,templateFName):
    self.templateFName= IoUtil.prependCWD(templateFName)
    
    # This external resource is found dynamically, but must be in the Python search path:
    scriptFileName= "rocs_score1.sh"
    self.scriptName= IoUtil.findFileInPyPath(scriptFileName)
    if (self.scriptName is None):
      raise Exception("ROCS driver shell script: " + scriptFileName + " not found")
    else:
      print("Located ROCS shell script: " + self.scriptName)
      
    self.queueName= "blade.q"          # SGE queue to use for parallel execution
    self.maxWaitSGE= 1 * 60.0 * 60.0   # Max time to wait for results to come back from SGE queue (secs)
    
  # ---------------------------------------------------------------------------
  """
    This method parses an OpenEye ROCS report file to retrieve the desired score metric from
    the first line of values - therefore the output option to rank by this desired metric
    should be used in the ROCS run.  
    
    NOTE: Both the query/template molecule and the molecule to be scored must have an identifier
    (molname) with a single token, when split by whitespace, so that column names may be matched
    to values in this report file.
   
    @param inFileName: filename of file to be parsed
    @param metricHdr: Header string that can be matched from the first line to identify the desired
        output column (case sensitive)
   
    @returns: the desired value from the top ranked alignment
  """
  def parseRocsReport1(self,inFileName,metricHdr):

    #print "Parsing ROCS report"
    inFile= open(inFileName,'r')

    # Read and process header line to locate desired column
    inLine= inFile.readline()
    hdrs= inLine.split('\t')

    i= 0
    col= -1
    for hdr in hdrs:
      if hdr == metricHdr:
        col= i
        break
      i= i+1

    # If desired header found, then read top values line and pick out the desired value
    if (col > -1):
      inLine= inFile.readline()
      strgValues= inLine.split('\t')
      strgValue= strgValues[col]
      inFile.close() 
    else:
      inFile.close() 
      raise LookupError("Expected column header: " + metricHdr + " not found in file: " + inFileName)

    value= float(strgValue)
    if (metricHdr == "ComboScore"):        # the ComboScore is out of 2
      value= value/2.0
    
    return(value)    


  # ---------------------------------------------------------------------------------------
  """
    Implementation of the standard value method, here corresponding to the ROCS score
      @param basicMol: a BasicMol to be evaluated, molObject is ignored
  """
  def value(self,basicMol):

    #print "Entering ROCSScoreFunc.value()"
    # Write out smiles string for the product to a file
    baseName= IoUtil.prependCWD("swarm_tmp")
    smilesFName= baseName + ".smi"
    reportFName= baseName + ".rpt"
    
    # Create the input file for scoring
    smilesFile= open(smilesFName, 'w')
    smilesFile.write(basicMol.smiStrg)
    smilesFile.close()

    # Invoke the shell script, with the new file specified as input, and the output filename as well
    cmdLine= self.scriptName + " " + smilesFName + " " + self.templateFName
    cmdTokens= cmdLine.split()  
    #print "Command Line Echo: ", cmdLine
    osProcess= Popen(cmdLine, shell=True)
    
    # osProcess= Popen(cmdTokens, shell=True)
    returnCode= osProcess.wait()
  
    # Cleanup - shell script should delete intermediate files
    # os.remove(smilesFName)

    value= -1.
    if (returnCode == 0):  # Read the final output file
      if (os.path.isfile(reportFName)):                # check if file exists first
        value= self.parseRocsReport(reportFName,"ComboScore")
      
        # Cleanup - shell script should delete intermediate files
        os.remove(reportFName)
      else:
        # This should really be an exception, but to keep long runs proceeding...
        print("Error: " + reportFName + " not found in ROCSScoreFunc.value... ignoring")
        value= 0.0
      
    else:
      raise OSError("OS Commandline failed: " + cmdLine)

    print("Leaving ROCSScoreFunc.value(), value= " + "%5.4f" % value)

    return(value)

    
  # ----------------------------------------------------------------------------
  """
    This is the multi-mol version of parseRocsReport1, reads through the entire
    ROCS report file, extracting ids and scores, and keeping the highest score
    for each ID found.
   
    @param inFileName: filename of file to be parsed
    @param metricHdr: Header string that can be matched from the first line to 
       identify the desired output column (case sensitive)
   
    @returns: the desired value from the top ranked alignment
  """
  def parseRocsReport2(self,inFileName,metricHdr):

    resultsDict= dict()
    
    #print("...parsing ROCS report")
    try:
      inFile= open(inFileName,'r')
    except IOError:
      print("Error: " + inFileName + "not found in ROCSScoreFunc.parseRocsReport2... ignoring")
      return(resultsDict)

    # Read and process header line to locate desired column
    inLine= inFile.readline()
    hdrs= inLine.split()
    #print("Found: " + str(len(hdrs)) + " tokens on header line")
    
    i= 0
    col= -1
    for hdr in hdrs:
      if (hdr == metricHdr):
        col= i
        break
      i= i+1

    # If desired header found, then read top values line and pick out the desired value
    if (col > -1):
      #print("Found metricHdr at column: " + str(col))
      more= True
      while (more):
        inLine= inFile.readline()
        inLine= inLine.strip()
        if (inLine is None or len(inLine) == 0):
          more= False
        else:
          tokens= inLine.split()
          #print "Found: " + str(len(tokens)) + " tokens on report line"

          molName= tokens[0]               # extract the name prefix from the conformer id suffix
          pos= molName.rfind("_")
          if (pos > 0):
            molName= molName[0:pos]
          
          strgValue= tokens[col]
          value= float(strgValue)
          if (metricHdr == "ComboScore"):        # the ComboScore is out of 2
            value= value/2.0
 
          # Can have multiple entries for different enantiomers - keep the highest
          existingValue= resultsDict.get(molName)
          if (existingValue != None):
            if (value > existingValue):
              resultsDict.update({molName:value})
          else:
            resultsDict.update({molName:value})
 
      inFile.close()
      
    else:
      inFile.close() 
      errStrg= "Expected ROCS report column header: " + metricHdr + " not found in file: " + inFileName + " in ROCSScoreFunc.parseRocsReport2()... ignoring"
      IoUtil.logErrorString(errStrg,True)
      # raise LookupError("Expected ROCS report column header: " + metricHdr + " not found in file: " + inFileName)

    return(resultsDict)    

  # ---------------------------------------------------------------------------------------
  """
    NOTE: This version of values does not use the new parallelism.  It does though write
    out a list of structures for scoring, and properly parse the multimolecule score report.
    
    Ensemble version of ROCS scoring value() method to allow groups to be evaluated at one time
    more efficiently through parallelization.

    @param basicMolList: a list of BasicMol objects to be evaluated
    @returns: a dict object, keyed by the molecule names passed in, with values equal
      to their corresponding scores
  """
  def valuesSingleThread(self,basicMolList):
  
    t1= time.time()
    print("Entering ROCSScoreFunc.values()")

    # Write out smiles string for the product to a file
    baseName= IoUtil.prependCWD("swarm_tmp")
    smilesFName= baseName + ".smi"
    reportFName= baseName + ".rpt"
    
    # Create the input file for scoring
    smilesFile= open(smilesFName, 'w')
    for basicMol in basicMolList:
      smilesFile.write(basicMol.smiStrg + "\t" + basicMol.name + "\n")
      
    smilesFile.close()

    # Invoke the shell script, with the new file specified as input, and the output filename as well
    cmdLine= self.scriptName + " " + smilesFName + " " + self.templateFName
    
    cmdTokens= cmdLine.split()  
    #print "Command Line Echo: ", cmdLine
    osProcess= Popen(cmdLine, shell=True)
    
    # osProcess= Popen(cmdTokens, shell=True)
    returnCode= osProcess.wait()
  
    # Cleanup - shell script should delete intermediate files
    # os.remove(smilesFName)

    if (returnCode == 0):  # Read the final output file
      if (os.path.isfile(reportFName)):                # check if file exists first
        resultsDict= self.parseRocsReport2(reportFName,"ComboScore")
      
        # Cleanup - shell script should delete intermediate files
        os.remove(reportFName)
      else:
        # This should really be an exception, but to keep long runs proceeding..
        errStrg= "Error: " + reportFName + " not found in ROCSScoreFunc.valuesSingleThread()... ignoring"
        IoUtil.logErrorString(errStrg,True)
        resultsDict= dict()   # empty dictionary

    else:
      raise OSError("OS Commandline failed: " + cmdLine)


    print("...leaving")

    t2= time.time()
    print(str(len(basicMolList)) + " mols scored in " + str(t2-t1) + " seconds")

    return(resultsDict)

  # ---------------------------------------------------------------------------------------
  """
    This is strictly for debugging - profile molecules with scores not found as to in which
    chunks, and in what positions in those chunks to error trace file to allow distribution
    of failures to be examined. Then stop run altogether to allow examination of the 
    intermediate files from parallel scoring which correspond to chunks with an error - 
    especially the SMILES files, which can be run through rocs_score1.sh manually to examine 
    if a failure occurs there and why.
    
    To use, comment out intermediate file deletion in values() and call this after the allResultsDict
    has been completely filled in.
  """
  def chunkChecker(self,chunkLists,allResultsDict):
    IoUtil.logErrorString("******Chunk Checker Output of Missing Scores******\n")
    for i in xrange(len(chunkLists)):
      IoUtil.logErrorString("Chunk " + str(i) + " ----------------------------------\n")
      IoUtil.logErrorString("  Chunk " + str(i) + " has " + str(len(chunkLists[i])) + " molecules\n")
      numMissing= 0
      for j in xrange(len(chunkLists[i])):
        score= allResultsDict.get(chunkLists[i][j].name)
        if (score is None):
          numMissing= numMissing + 1
          IoUtil.logErrorString("  Mol [" + str(i) +"][" + str(j) + "]  missing score\n")
      IoUtil.logErrorString("  Chunk " + str(i) + " missing " + str(numMissing) + " scores in all\n")
          
    IoUtil.logErrorString("******End Chunk Checker Output ******\n")
    
    IoUtil.closeErrorFile()
    if (chunkLists != None):
      sys.exit(0)

# ---------------------------------------------------------------------------------------
  """
    Ensemble version of ROCS scoring value() method to allow groups to be evaluated at one time
    more efficiently through parallelization.
    
    @param basicMolList: a list of BasicMol objects to be evaluated
    @returns: a dict object, keyed by the molecule names passed in, with values equal
      to their corresponding scores
  """
  def values(self,basicMolList):

    print("Entering ROCSScoreFunc.valuesParallel()")
    
    #  Include path in filenames, because qsub will run them from another machine
    baseName= IoUtil.prependCWD("swarm_tmp")
    sgeUtils= SGEUtils()
    
    # Decide on number of chunks:
    """
      Unfortunately, runtime per mol varies greately with:
        - num chiral centers (supralinear)
        - num rotors  (supralinear)
        - num atoms
      ... so runtime per mol can easily vary  by 100x... I've seen it from < 1 sec per mol to
      1 minute per mol

      Hmmm... there seems to be roughly a 15 second penalty just for using the cluster      
    """

    availProcessors= sgeUtils.getAvailableProcCount(self.queueName)
    print("Found " + str(availProcessors) + " processors available in queue: " + self.queueName)
    maxProcessorsToUse= int(0.80 * availProcessors)     # maxmimum number to run in parallel
    
    numMols= len(basicMolList)
    
    minChunkSize= 3    # minimum molecules in file per submission to the SGE queue
    maxNumChunks= maxProcessorsToUse
    if (numMols <= minChunkSize):
      numChunks= 1
    else:
      numChunks= (numMols/minChunkSize)
      if (numChunks > maxNumChunks):
        numChunks= maxNumChunks
    
    print("Using " + str(numChunks) + " chunks")
    
    # Reformat single list of mols to list of numChunks sublists
    chunkLists= self.breakIntoChunks(basicMolList,numChunks)
    # Write out numChunks SMILES files
    smiNames= self.writeChunks(chunkLists,baseName,".smi")
    
    # Generate names of the corresponding report files matching those generated by the ROCS script, 
    # which we'll read back in
    rptNames= list()
    for fName in smiNames:
      rptNames.append(fName.replace(".smi",".rpt"))
      
    print("Submitting " + str(len(smiNames)) + " simultaneous jobs to SGE through qsub")
    workingDir= os.getcwd()
    qsubBase= "qsub -r yes -p 0 -S /bin/bash -V -q " + self.queueName + " -e " + workingDir + " -o " + workingDir
    jobNames= set()
    ind= 0
    for smiName in smiNames:
      jobName= "swm_" + str(ind)
      ind= ind + 1
      jobNames.add(jobName)
      cmdLine= qsubBase + " -N " + jobName + " " + self.scriptName + " " + smiName + " " + self.templateFName
      cmdTokens= cmdLine.split()
      osProcess= Popen(cmdLine,shell=True)
      
    print("All submitted, now polling for SGE completion")
    rocsTimer= Timer()
    try:
      sgeUtils.pollQueueUntilDone(None,jobNames,self.queueName,self.maxWaitSGE)
    except TimeoutException(ex):
      IoUtil.logErrorString(ex.__repr__(),True)
    else:
      print("All SGE jobs completed in " +  str(rocsTimer.getElapsed()) + " seconds")

    #----
    # Okay, NOW, I need to parse all that output! 
    
    
    # First clean up all the temp SMILES files that were written    
    for smiName in smiNames:
      if (os.path.isfile(smiName)):        # check if file exists first
        os.remove(smiName)
        
    # Now clean up all the temp stderr, stdout files that were written by the SGE during execution
    # Using wildcard operator with workingDir/jobname base
    for jobName in jobNames:
      cmdLine= "rm -f " + workingDir + os.sep + jobName + ".*"
      cmdTokens= cmdLine.split()
      osProcess= Popen(cmdLine,shell=True)
      returnCode= osProcess.wait()
    
    print("Parsing output files")
    allResultsDict= dict()
    for i in xrange(len(rptNames)):
      reportFName= rptNames[i]
      if (os.path.isfile(reportFName)):                # check if file exists first
        resultsDict= self.parseRocsReport2(reportFName,"ComboScore")
      
        # Add chunk results to total results
        keyIter= resultsDict.iterkeys()
        for key in keyIter:
          value= resultsDict.get(key)
          allResultsDict.update({key:value})
        resultsDict.clear()
      
        # Cleanup temporary ROCS report files - script should have gotten the rest
        os.remove(reportFName)
      else:
        # This should really be an exception, but to keep long runs proceeding...
        errStrg= "Error: " + reportFName + " not found in ROCSScoreFunc.values()... ignoring"
        IoUtil.logErrorString(errStrg,True)
        
    #self.chunkChecker(chunkLists,allResultsDict)
      
    print("...leaving")

    return(allResultsDict)

  # ---------------------------------------------------------------------------------------
"""
    Simple test of parallel execution - 
      - does it run?
      - how long does it take?
      - what mols get values returned?
"""
def testParallelExec():

    templateFName= raw_input("Template molecule infile name (.mol/.pdb/.mol2): ")
    if (os.path.isfile(templateFName)):
      app= ROCSScoreFunc(templateFName)
    else:
      print("ERROR: file not found: " + templateFName)
      sys.exit(1)
    
    molsFName= raw_input( "SMILES test molecules infile name: " )
    if (os.path.isfile(molsFName)):
      inFile= file(molsFName,'r')
      molColl= MolCollection()
      molColl.loadFromFile(inFile)
    else:
      print("ERROR: file not found: " + molsFName)
      sys.exit(1)

    evalTimer= Timer()
    resultsDict= app.values(molColl.molList)
    print("Scored " +  str(len(molColl)) + " mols in parallel in " + str(evalTimer.getElapsed()) + " seconds")
    
    for mol in molColl:
      score= resultsDict.get(mol.name)
      if (score is None):
        mol.prodScore= -1
        print("No value returned for: " + mol.name)
      else:
        mol.prodScore= score

    outFName= raw_input("Output file name [None]: ")
    if (outFName != None):
      outFName= outFName.strip()
      if (len(outFName) > 0):
        outFile= file(outFName,'w')
        molColl.writeToFile(outFile,True)

  # ---------------------------------------------------------------------------------------
"""
    Simple test of ROCS report parsing, uses a pre-existing SMILES file and corresponding
    ROCS report file
"""
def testParse2():
    app= ROCSScoreFunc("dummy")

    molsFName= raw_input( "SMILES test molecules infile name: " )
    if (os.path.isfile(molsFName)):
      inFile= file(molsFName,'r')
      molColl= MolCollection()
      molColl.loadFromFile(inFile)
    else:
      print("ERROR: file not found: " + molsFName)
      sys.exit(1)

    rptFName= raw_input( "Matching report infile name: " )
    if (os.path.isfile(rptFName)):
      resultsDict= app.parseRocsReport2(rptFName,"ComboScore")
    else:
      print("ERROR: file not found: " + molsFName)
      sys.exit(1)

    for mol in molColl:
      score= resultsDict.get(mol.name)
      if (score is None):
        mol.prodScore= -1
        print("No value returned for: " + mol.name)
      else:
        mol.prodScore= score
    
    outFName= raw_input("Output file name [None]: ")
    if (outFName != None):
      outFName= outFName.strip()
      if (len(outFName) > 0):
        outFile= file(outFName,'w')
        molColl.writeToFile(outFile,True) 

  
#==============================================================================
# Test Driver code
if __name__ == "__main__":
  testParallelExec()
  #testParse2()
    
    