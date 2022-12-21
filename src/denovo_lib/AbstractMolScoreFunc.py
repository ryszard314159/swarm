
import math

""" 
  AbstractMolScoreFunc.py: Abstract baseclass for molecule scoring functions
  Constructors may be used to initialize the concrete scoring functions with options
""" 
class AbstractMolScoreFunc:
    
  # ---------------------------------------------------------------------------
  """ 
    The money method - must be implemented by all concrete subclasses
    @param basicMol: a BasicMol object to be evaluated
    @returns: a score in the range [0..1] -> 1 = best
  """ 
  def value(self,basicMol):
    raise NotImplementedError("value must be defined in subclasses of AbstractMolScoreFunc!")

  # ---------------------------------------------------------------------------
  """ 
    Ensemble version of value() method to allow groups to be evaluated at one time
    more efficiently through, e.g., parallelization.  This default implementation
    just calls the single mol method in a loop - override as desired to implement
    parallelization.

    @param basicMolList: a list of BasicMol objects to be evaluated
    @returns: a dict object, keyed by the molecule names passed in, with values equal
      to their corresponding scores
  """ 
  def values(self,basicMolList):
    resultsDict= dict()
    for basicMol in basicMolList:
      score= self.value(basicMol)
      resultsDict.update({basicMol.name:score})
      
    return(resultsDict)
       
  # ---------------------------------------------------------------------------
  """
    Returns the match of the rawValue passed in to a Gaussian filter function
    parameterized as specified.

    Note that this is not identical to the amplitude of the corresponding normal 
    distribution's amplitude, because the normalizing preterm of (1/sqrt(2pi*stdDevn))
    is omitted purposely here, so that the maximum output value will be 1.0

    @param rawValue: value to be processed
    @param mean: mean for the Gaussian function
    @param stdDevn: standard deviation for the Gaussian function
    @returns: exp(-sqr(x-xbar)/(2*sqr(std)))  
  """
  def gaussianMatch(self,rawValue,mean,gaussStdDevn):  
    ssq= (rawValue - mean) * (rawValue - mean)
    if (ssq == 0):
      return(1.0)
    else:
      return(math.exp(-(ssq/(2.0 * stdDevn * stdDevn))))

  # ---------------------------------------------------------------------------
  """
    Breaks up a list of molecules into a list of sublists (chunks).
    
    @param molList: a list of molecules, actually, this will for any list of objects ;-)
    @param nChunks: the desired number of chunks
    @returns: a list of nChunks lists... the list passed in is unchanged
  """
  def breakIntoChunks(self,molList,nChunks):
    if (len(molList) < nChunks):
      numChunks= len(molList)
    else:
      numChunks= nChunks
      
    outerList= [None] * numChunks         # Create list of numChunks object pointers 
    for i in xrange(numChunks):           # Now fill in each element as a sublist
      outerList[i]= list()
      
    molInd= 0
    listInd= 0
    while (molInd < len(molList)):
      outerList[listInd].append(molList[molInd])
      molInd= molInd + 1
      listInd= listInd + 1
      if (listInd == numChunks):          # Cycle this index
        listInd= 0
 
    return(outerList)

  # ---------------------------------------------------------------------------
  """
    Writes out a series of SMILES files corresponding to a chunk list returned 
    by breakIntoChunks(), above.  Only the SMILES string and name are written,
    to avoid confusing external programs.
    
    @param chunkLists: a list of lists of molecules
    @param baseName: the base name for the files to be written out
    @param extension: the extension for the filenames to be written out, e.g., ".smi"
    @returns: a list of filenames for the files written out
  """
  def writeChunks(self,chunkLists,baseName,extension=""):
    fNameList= list()
    ind= 0
    for subList in chunkLists:
      fName= baseName + str(ind) + extension
      ind= ind + 1
      fNameList.append(fName)
      smilesFile= open(fName, 'w')
      for basicMol in subList:
        smilesFile.write(basicMol.smiStrg + "\t" + basicMol.name + "\n")
      
      smilesFile.close()
      
    return(fNameList)

      