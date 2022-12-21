#!/usr/bin/env python

import sys
import random

# from openeye.oechem import *

from denovo_lib.MolCollection import MolCollection
from denovo_lib.ScoredMol import ScoredMol
from denovo_lib.SynthesizedMol import SynthesizedMol
from denovo_lib.ReactionTransformOE import ReactionTransformOE
from denovo_lib.RefiningDenovoOptimizer import RefiningDenovoOptimizer
from denovo_lib.ROCSScoreFunc import ROCSScoreFunc


"""
  Swarm.py: Main program module for SWARM (Stochastic Wandering of Actively Reacting Molecules)
"""
class Swarm:

  # ---------------------------------------------------------------------------
  """
    Loads the molecules from a file of Daylight SMILES strings, establishing the collection
    @returns: a MolCollection
  """
  def readMols(self,molsFName):
    try:
      inFile = file(molsFName,'r')
    except Exception:
      print("Error opening molecules file: " + molsFName)
      raise
      
    mols= MolCollection()
    mols.loadFromFile(inFile)

    return(mols)

  # ---------------------------------------------------------------------------
  """
    Loads the reactions from a file of Daylight SMIRKS strings
    @returns: a list of ReactionTransformOE objects
  """
  def readReaxns(self,rxnsFName):
    try:
      inFile = file(rxnsFName,'r')
    except Exception:
      print("Error opening reactions file: " + rxnsFName)
      raise

    try:
      smirksStrings= inFile.readlines()
    finally:
      inFile.close()
    
    rxns= list()
    for smirksString in smirksStrings:
      smirksString= smirksString.rstrip()
      if (len(smirksString) > 0 and smirksString[0] != '#'):
        tokens= smirksString.split()
        smirks= tokens[0]
        if (len(tokens) > 1):
          name= tokens[1]
        else:
          name= ""
          
        rxns.append(ReactionTransformOE(smirks,name))

    print("Read " + str(len(rxns)) + " reactions")

    return(rxns)

  # ---------------------------------------------------------------------------
  """
    Simple surrogate for scoring, for debugging purposes only
  """
  def fakeScoring(self,molList):
    resultsDict= dict()
    for mol in molList:
      score= random.random()
      resultsDict.update({mol.name:score})
      
    return(resultsDict)
    
  # ---------------------------------------------------------------------------
  """
    Rescores the top N molecules (as judged by the current scoring function),
    with another scoring function as desired.
    
    @param molCollection: the MolCollection to be evaluated
    @param scoringFunc: the AbstractMolScoreFunc to be used
    @param n: integer for the size of subset to be taken

    @returns: a list containing the subset of molecules with new product scores
  """
  def scoreTopN(self,molCollection,scoringFunc,n):
    subList= molCollection.getNBestProducts(n)
    
    resultsDict= scoringFunc.values(subList)
    #resultsDict= self.fakeScoring(subList)

    for mol in subList:
      scoreValue= resultsDict.get(mol.name)
      if (scoreValue is None):            # No score was returned for this molecule
        mol.prodScore= 0.0                
      else:
        mol.prodScore= scoreValue

    subCollection= MolCollection()
    subCollection.appendList(subList)
    return(subCollection)

  # ---------------------------------------------------------------------------
  """
    Main program execution
  """
  def run(self, opts=None):

    # Get the molecular building blocks
    if opts['reagents'] != None:
      molsFName = opts['reagents']
    else:
      molsFName= raw_input( "SMILES building block molecules infile name: " )
    molCollection= self.readMols(molsFName)
  
    # Get reactions
    if opts['transforms'] != None:
      rxnsFName = opts['transforms']
    else:
      rxnsFName= raw_input( "SMIRKS reactions infile name: " )
    rxnList= self.readReaxns(rxnsFName) 
    
    # Create an optimizer, set it up, and then set it to work
    optimizer= RefiningDenovoOptimizer(molCollection,rxnList)
    optimizer.setup(opts)    
    finalCollection= optimizer.optimize()

    # Sort them into descending order by product score
    finalCollection.sortByProductScore()
    finalCollection.reverse()
    
    # Now write them out    
    print("Done structural exploration, writing out returned scores")

    #finalCollection.writeToFile(file("swarm_brief_lingo.out",'w'),True)
    #finalCollection.writeToFile(file("swarm_full_lingo.out",'w'),False)
    
    finalCollection.writeToFile(file("swarm_brief.out",'w'),True)
    finalCollection.writeToFile(file("swarm_full.out",'w'),False)
      
    
    # EXTRA: After Lingo in optimization, rescore the best N with ROCS here and write that out also...
    """
    print "Calculating ROCS scores"
    scoringFunc= ROCSScoreFunc(optimizer.templateMolFName)
    numRescore= 2000    
    subCollection= self.scoreTopN(finalCollection,scoringFunc,numRescore)
    subCollection.sortByProductScore()
    subCollection.reverse() 
    
    subCollection.writeToFile(file("swarm_brief_rocs.out",'w'),True)
    subCollection.writeToFile(file("swarm_full_rocs.out",'w'),False)

    """

#------------------------------------------------------------------------------

def usage():
      """
      swarmp.py -h|--help              # print this help
                -r|--reactants smiles  # file with reactants; default=None
                -x|--transforms smirks # file with reactions; default=None
                -s|--score scoreName   # score name: Lingo|ROCS; default=Lingo
                --seed                 # seed for random number generator; default=None (i.e. initialized from the clock)
                -m|--maxProducts int   # maximum number of products to synthetize; default=1e99 
                -t|--maxRuntime float  # time limit in hours; default=14
                --template molfile     # template (query) molecule for ligand-similarity based scoring (absolute path to mol|pdb|mol2); default=None
                -f|--saveFreq float    # frequency of saving results in hours; default=1
                -p|--prefix prefix     # Mol id prefix for synthesized mols (default=ISS)
                -d|--debug
      """

# Driver code
if __name__ == "__main__":
   # parse a command line and create opts
   from getopt import getopt
   clopts, args = getopt(sys.argv[1:], 'hdf:r:x:p:s:m:t:', ['help', 'reagents=', 'transforms=', 'prefix=', 'score=', 'maxProducts=', 'maxRuntime=', 'template=',
                                                            'saveFreq=', 'seed='])
   opts = {'debug': 0, 'scoreName': 'Lingo', 'maxProducts': 1e99, 'maxRuntime': 14, 'saveFreq': 1, 'reagents': None, 'transforms': None, 'template': None,
           'prefix': 'ISS', 'seed': None}
   for o, a in clopts:
     if o in ['-h', '--help']:
       print(usage.__doc__)
       sys.exit(1)
     elif o in ['-d', '--debug']:
       opts['debug'] = 1
     elif o in ['-r', '--reagents']:
       opts['reagents'] = a
     elif o in ['-x', '--transforms']:
       opts['transforms'] = a
     elif o in ['--template']:
       opts['template'] = a
     elif o in ['-s', '--score']:
       opts['scoreName'] = a
     elif o in ['--seed']:
       opts['seed'] = eval(a)
     elif o in ['-p', '--prefix']:
       opts['prefix'] = a
     elif o in ['-m', '--maxProducts']:
       opts['maxProducts'] = eval(a)
     elif o in ['-t', '--maxRuntime']:
       opts['maxRuntime'] = eval(a)
     elif o in ['-f', '--saveFreq']:
       opts['saveFreq'] = eval(a)

   app = Swarm()
   app.run(opts)
