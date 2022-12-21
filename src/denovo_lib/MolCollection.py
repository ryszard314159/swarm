import math
import random
import gc

from SynthesizedMol import SynthesizedMol
from ScoredMol import ScoredMol
from BasicMol import BasicMol
from IoUtil import IoUtil

"""
  MolCollection.py: Implements a collection of ScoredMol objects, with iteration,
  efficient access to sorting (product score, reagent score, and SMILES string), 
  random, best, and random-best mols, plus fast membership testing.  Also
  implemented are I/O, subsetting through filters.
  
  This is done by using both a list and a set instance vars.  There is little memory
  overhead for doing this though, as only the pointers to the molecules are redundant.  
"""

class MolCollection:
  def __init__(self):
    self.molSet= set()      # Allows efficient membership testing by SMILES, using BasicMol.__hash__()
    self.molList= list()    # Allows facile iteration, ranking
    self.isPreScored= False    # Handle and flag the loading of pre-scored molecules with 
                               # their scores from a previous run with brief format output

  # ---------------------------------------------------------------------------
  """
    This allows use of the standard 'in' operator on MolCollection to test for set membership
    
    @param basicMol: a BasicMol object
    @returns: True, if a molecule already exists in this collection with the same SMILES string
  """
  def __contains__(self,basicMol):
    return(self.molSet.__contains__(basicMol))

  # ---------------------------------------------------------------------------
  """
    This plus the MolCollectionIter class, below, allows use of standard iteration
    syntax, e.g., 'for mol in coll' on MolCollection

    @returns: a MolCollectionIter
  """
  def __iter__(self):
    return(MolCollectionIter(self.molList))

  # ---------------------------------------------------------------------------
  """
    This allows use of standard [i] indexing operator for read access to elements of
    MolCollection, with integer indexing only.
    
    @param ind: An integer
    @returns: the [ind-th] molecule in this MolCollection
    
    See also:  __setitem__(), .__delitem__(), not yet implemented due to presence of
    mySet.

  """
  def __getitem__(self,ind):
    if (type(ind) is int):
      return(self.molList.__getitem__(ind))
    else:
      raise TypeError("Invalid indexing key for MolCollection.__getitem__(): " + str(ind))


  # ---------------------------------------------------------------------------
  """
    Allows use of the standard global function len() on MolCollection
 
    @returns: the number of molecules in this MolCollection
  """
  def __len__(self):
    return(self.molList.__len__())


  # ---------------------------------------------------------------------------
  """
    Add function, which checks the type being added (better be ScoredMol or subclass),
    which throws an exception, as well as for previous presence of one with same SMILES
    string, which simply returns with no change.
  """
  def add(self,scoredMol):
    
    if (isinstance(scoredMol,ScoredMol)):
      if (not scoredMol in self.molSet):
        self.molSet.add(scoredMol)
        self.molList.append(scoredMol)
    else:
      raise TypeError("Invalid object type passed to MolCollection.add(): " + str(type(scoredMol)))

  # ---------------------------------------------------------------------------
  """
    Clears all molecules from this MolCollection
  """
  def clear(self):
    del(self.molSet)
    del(self.molList)
    self.molSet= set()
    self.molList= list()
    gc.collect()

  # ---------------------------------------------------------------------------
  """
    Currently this is creating SynthesizedMol objects in the collections from
    the SMILES strings read in. 
    
    @param inFile: file handle for a file of SMILES strings
    
    This has been upgraded to read back product and reagent scores if present, but not
    yet full synthetic routes.
  """
  def loadFromFile(self,inFile):
    prodScore = 0
    errcount = 0
    lineNum= 0
    try:
      for line in inFile:
        lineNum= lineNum+1
        line= line.strip()      # drop leading whitespace, trailing eoln chars - cause SMILES parse error
        if (len(line) > 0):
          if (line[0] == "#"):
            continue
          
          tokens= line.split()
          smiles= tokens[0]

          if (len(tokens) > 1):
            name= tokens[1]
          else:
            name= ""

          # Read in any pre-existing scores: 
          if (len(tokens) > 2):
            try:
              prodScore= float(tokens[2])
              self.isPreScored= True
            except:
              errcount = errcount + 1
          else:
            prodScore= 0.0

          rgtScore= 0.5
          if (len(tokens) > 3):
            try:
              rgtScore= float(tokens[3])
            except:
              errcount = errcount + 1

          # Get rid of any 'extra' graphs by choosing largest
          if (smiles.find(".") > 0):
            tokens= smiles.split(".")
            maxLen= 0
            maxToken= None
            for token in tokens:
              if (len(token) > maxLen):
                maxToken= token
                maxLen= len(token) 
            smiles= maxToken

          try:          
            synMol= SynthesizedMol(name,smiles,None,prodScore,rgtScore)            
            synMol.assignCanonicalSmiles(True)
            self.add(synMol)
          except Exception:
            errmsg = "SMILES parsing error at line: " + str(lineNum) + " for SMILES= " + line
            print("SMILES parsing error at line: " + str(lineNum) + " for SMILES= " + line)
            #raise
            IoUtil.logErrorString(errmsg)
        
    finally:
      inFile.close()
    
    print("Read in " + str(self.size()) + " molecules")
    if errcount > 0: print("loadFromFile: ', errcount, ' errors while reading scores")

  # ---------------------------------------------------------------------------
  """
    Currently this is assuming SynthesizedMol objects in the collections. 
    
    @param outFile: file handle fro the output
    @param briefFormat: boolean, if True, the mols are written out in brief
      format - just the SMILES string, name and scores on 1 line; if False,
      their synthetic origin is written out as well
  """
  def writeToFile(self,outFile,briefFormat):

    try:
      if (briefFormat):
        outFile.write("# SMILES\tMOL_ID\tPRODUCT_SCORE\tREAGENT_SCORE\n")
        for mol in self.molList:
          mol.printBrief(outFile)
      else:
        outFile.write("# SMILES\tMOL_ID\tPRODUCT_SCORE\tREAGENT_SCORE\n")
        outFile.write("# Synthetic origin is indented under main product line \n")
        for mol in self.molList:
          mol.printFull(outFile,indentation="")

    finally:
      outFile.close()


  # ---------------------------------------------------------------------------
  """
    Appends a list of mols to the current contents.  See also clear()
    
    @param aMolList: a list of ScoredMols
  """
  def appendList(self,aMolList):
    for mol in aMolList:
      self.add(mol)
  
  # ---------------------------------------------------------------------------
  """ Sorts the molList into ascending order by their product score values """
  def sortByProductScore(self):
    self.molList.sort(ScoredMol.compareOnProdScore)
  
  # ---------------------------------------------------------------------------
  """ Sorts the molList into ascending order by their reagent score values """
  def sortByReagentScore(self):
    self.molList.sort(ScoredMol.compareOnRgtScore)

  # ---------------------------------------------------------------------------
  """ Sorts the molList into ascending order by their smiles string values """
  def sortBySmilesString(self):
    self.molList.sort(ScoredMol.compareOnSmiStrg)

  # ---------------------------------------------------------------------------
  """
    Useful post-sorting (which yields ascending order), to get descending order.    
  """
  def reverse(self):
    self.molList.reverse()

  # ---------------------------------------------------------------------------
  """
    @returns: integer giving the number of molecules held in this MolCollection
  """ 
  def size(self):
    return(self.__len__())    
  
  #----------------------------------------------------------------------------
  """
    @returns: the molecule wth the highest reagent score in this MolCollection
   """
  def getBestRgt(self):
    bestMol= self.molList[0]
    bestScore= bestMol.getReagentScore()
    
    for scoredMol in self.molList:
      if (scoredMol.getReagentScore() > bestScore):
        bestScore= scoredMol.getReagentScore()
        bestMol= scoredMol
      
    return(bestMol)

  #----------------------------------------------------------------------------
  """
    @returns: the molecule wth the highest product score in this MolCollection
   """
  def getBestProd(self):
    bestMol= self.molList[0]
    bestScore= bestMol.prodScore
    
    for scoredMol in self.molList:
      if (scoredMol.prodScore > bestScore):
        bestScore= scoredMol.prodScore
        bestMol= scoredMol
      
    return(bestMol)

  # ---------------------------------------------------------------------------
  """
    @returns: a list containing the num molecules with the highest product score 
      in this MolCollection
  """
  def getNBestReagents(self,num):
    if (num > self.size()):
      num= self.size()
      
    self.sortByReagentScore()

    last= self.size()-1
    first= last - (num-1)
    if (first < 0):
      first= 0

    return(self.molList[first:last])

  # ---------------------------------------------------------------------------
  """
    @returns: a list containing the num molecules with the highest product score 
      in this MolCollection
  """
  def getNBestProducts(self,num):
    if (num > self.size()):
      num= self.size()

    self.sortByProductScore()

    last= self.size()-1
    first= last - (num-1)
    if (first < 0):
      first= 0

    return(self.molList[first:last])

  # ---------------------------------------------------------------------------
  """
    @returns: a randomly chosen molecule held in this MolCollection
  """
  def getRandom(self):
    return(self.molList[random.randint(0,self.size()-1)])
    
  #----------------------------------------------------------------------------
  """
    Metropolis sampling based on product score - molecules with scores < best will be 
    returned with probability decreasing according to the Metropolis sampling function.
    temperature must be > 0, and increasing temperature increases the likelihood of 
    accepting poorer performing molecules

    @param temperature: a float from 0 to MAX_FLOAT
  """
  def getRandMetropolisOnProd(self,temperature):
    bestMol= self.getBestProd()
    
    found= False
    while (not found):
       aMol= self.getRandom()
       if (aMol.prodScore >= bestMol.prodScore):
         found= True
       else:
         delta= bestMol.prodScore - aMol.prodScore
         randFlt= random.random()  # uniform random in [0..1]
         if (randFlt < math.exp(-delta/temperature)):
           found= True
    
    return(aMol)


  #----------------------------------------------------------------------------
  """
    Metropolis sampling based on reagent score - molecules with scores < best will be 
    returned with probability decreasing according to the Metropolis sampling function.
    temperature must be > 0, and increasing temperature increases the likelihood of 
    accepting poorer performing molecules

    @param temperature: a float from 0 to MAX_FLOAT
  """
  def getRandMetropolisOnRgt(self,temperature):
    bestMol= self.getBestRgt()
    
    found= False
    while (not found):
       aMol= self.getRandom()
       if (aMol.getReagentScore() >= bestMol.getReagentScore()):
         found= True
       else:
         delta= bestMol.getReagentScore() - aMol.getReagentScore()
         randFlt= random.random()  # uniform random in [0..1]
         if (randFlt < math.exp(-delta/temperature)):
           found= True
    
    return(aMol)


  # ---------------------------------------------------------------------------
  """
    Replaces the underlying set, list of molecules with the subset which pass the filter passed in
    
    @param molFilter: an AbstractMolFilter
  """
  def filter(self,molFilter):
    passList= list()
    passSet= set()
    for mol in self.molList:
      if (molFilter.passes(mol)):
        passList.append(mol)
        passSet.add(mol)
        
    self.molList= passList
    self.molSet= passSet

  # ---------------------------------------------------------------------------
  """
    Returns 2 MolColections - 
      1. those which pass the filter passed in
      2. those which do not.
      
    This MolCollection is left unchanged  

    @param molFilter: an AbstractMolFilter
  """
  def split(self,molFilter):
    passList= list()
    passSet= set()
    failList= list()
    failSet= set()
    for mol in self.molList:
      if (molFilter.passes(mol)):
        passList.append(mol)
        passSet.add(mol)
      else:
        failList.append(mol)
        failSet.add(mol)
        
    passed= MolCollection()
    passed.molList= passList
    passed.molSet= passSet
    
    failed= MolCollection()
    failed.molList= failList
    failed.molSet= failSet
    
    return(passed,failed)

  # ---------------------------------------------------------------------------
  # Prints a list of ScoredMols
  def printList(self):
    print("SMILES_STRG PROD_SCORE REAGENT_SCORE")
    for mol in self.molList:
      print(mol.__str__())



# =============================================================================
"""
  This is just a helper class to allow iteration syntax to be used on MolCollection
  objects, i.e., "for mols in coll"
"""
class MolCollectionIter:
  def __init__(self,aList):
    self.aList= aList
    self.ind= 0              # need an instance var that is a counter

  def next(self):
    if (self.ind == len(self.aList)):
      raise StopIteration
    else:
      locInd= self.ind
      self.ind= self.ind + 1
      
      return(self.aList[locInd])


#==============================================================================
# Test Driver code
if __name__ == "__main__":
   
  coll= MolCollection()
  coll.add(ScoredMol("Acetic Acid","CC(=O)O",None,random.random()))
  coll.add(ScoredMol("Ethyl Acetate","CC(=O)OCC",None,random.random()))
  coll.add(ScoredMol("Propyl Amine","CCCN",None,random.random()))
  coll.add(ScoredMol("2-Methyl Pyridine","c1cc(C)ncc1",None,random.random()))
  coll.add(ScoredMol("Methyl Ethyl Ketone","CC(=O)CC",None,random.random()))
  coll.add(ScoredMol("Propanoyl Chloride","CCC(=O)Cl",None,random.random()))

  print("Initial collection:")
  coll.printList()
  print("")
   
  coll.sortByProductScore()
  print("Sorted by product score:")
  coll.printList()
  print("")
   
  coll.sortBySmilesString()
  print("Sorted by smiles string:")
  coll.printList()
  print("")
   
  print("10 random samples from original list")
  for i in xrange(1,10):
    mol= coll.getRandom()
    print(mol.__str__())
  print("")

  inString= raw_input("Enter temperature for Metropolis test(0-1): ")
  temp= float(inString)
   
  freqDict= dict()
  for mol in coll:                  # Also test the iteration syntax here
    freqDict.update({mol.smiStrg:int(0)})
   
  for i in xrange(1,1000):
    mol= coll.getRandMetropolisOnProd(temp)
    num= freqDict.get(mol.smiStrg)
    num= num + 1
    freqDict.update({mol.smiStrg:num})

  print("After 1000 Metropolis random samples from original list")
  coll.sortByProductScore()
  for mol in coll:
    num= freqDict.get(mol.smiStrg)
    print(mol.__str__() + " was sampled " + str(num) + " times")

  print("Set membership testing")
  more= True
  while (more):
    smiStrg= raw_input("Enter a SMILES string: ")
    smiStrg= smiStrg.strip()
    if (len(smiStrg) < 1):
      more= False
    else:
      basicMol= BasicMol("",smiStrg,None)
      if (basicMol in coll):
        print("MolCollection has this molecule already")
      else:
        print("MolCollection does not have this molecule yet")
     
