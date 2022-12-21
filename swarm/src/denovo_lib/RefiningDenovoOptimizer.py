import sys
import random

# from openeye.oechem import *

from AbstractDenovoOptimizer import AbstractDenovoOptimizer
from BasicMol import BasicMol
from BasicFilter import BasicFilter
from CombiIter import CombiIter
from IDSequence import IDSequence
from IoUtil import IoUtil
from LingoScoreFunc import LingoScoreFunc
from ReactiveFilter import ReactiveFilter
from ROCSScoreFunc import ROCSScoreFunc
from ScoredMol import ScoredMol
from SynthesizedMol import SynthesizedMol
from SyntheticOrigin import SyntheticOrigin
from SyntheticFilter import SyntheticFilter
from Timer import Timer

"""
 RefiningDenovoOptimizer.py: attempts to refine initially judged best molecules

 Basic Rules:
 1. All molecules are scored as products to have a chance to be a candidate
 2. Molecules are chosen for improvement based on their product score, which comes from
    it's own evaluation against the template
 3. Molecules are chosen as reagents based on their reagent score, which comes from the product
    scores where they served as a reagent, averaged together
"""
class RefiningDenovoOptimizer(AbstractDenovoOptimizer):
  """
    @param molCollection: MolCollection of molecules (SynthesizedMol) to work with
    @param rxnList: List of reactions (ReactionTransformOE) to work with
  """
  def __init__(self,molCollection,rxnList):
    AbstractDenovoOptimizer.__init__(self,molCollection,rxnList)
    self.templateMol= None  # This is the ligand/substrate we are making an analog to
    self.maxProducts= 100000
    
    # 2 Complexity limits
    #     - To score - avoid even scoring anything really big (runtime)
    #     - To extend - avoid extending something that is already big
    #   Note: to score limits should be > to extend limits because anything too big to score is
    #         certainly too big to extend
    self.maxSizeRatio2Score= 1.6
    self.maxHeavyAtoms2Score= 30                  # computed from template size and maxSizeRatio
    self.maxChiralAtoms2Score= 6                  # 2^5 = 32 enantiomers to evaluate
    self.maxRotors2Score= 12                      # rotational DOF for conformational searching

    self.maxSizeRatio2Extend= 1.30
    self.maxHeavyAtoms2Extend= 23                  # computed from template size and maxSizeRatio
    self.maxChiralAtoms2Extend= 5                  # 2^4 = 16 enantiomers to evaluate
    self.maxRotors2Extend= 8                       # rotational DOF for conformational searching

    self.maxSynSteps= 20                    # synthetic steps - leave this loose - shorter routes
                                            # may be found manually
    

    self.maxRgtsPerRct= 10           # Max reagents per reactant in a given product expansion
    self.templateMolFName= None
    self.traceOutFName= "swarm.tra"       # Optional optimization trace file
    self.errOutFName= "swarm.err"         # Optional error trace file
    self.maxRuntime= 3600.0 * 14.0   # default runtime(secs) is 14 hours - overnight
    
    self.minToScorePll= 300          # minimum to accumulate before parallel scoring
    
    self.molIDPrefix= "ISS"
    self.idSeq= None                 # mol ID sequence for synthesized mols
    self.scoreName= 'Lingo'          # score function name
    #self.scoreName= 'ROCS'
    
  # ---------------------------------------------------------------------------
  """
    Sets up parameters using prompted terminal interaction with user unless defined on the command line
  """
  def setup(self, opts):

    #Use defaults for parameters

    if opts['debug']:
      print('opts =', opts)
    if opts['template'] != None:
      self.templateMolFName = opts['template']
    else:
      self.templateMolFName= raw_input("Template molecule filename (.mol/.pdb/.mol2) with path: ")

    inStream= oemolistream(self.templateMolFName)
    mol= OEGraphMol()    
    resultCode= OEReadMolecule(inStream,mol)
    if (resultCode == 1):
      smiStrg= OECreateCanSmiString(mol)
      self.templateMol= BasicMol("template",smiStrg,mol)
    else:
      raise Exception("Error parsing template molecule in file: " + self.templateMolFName)  

    if opts['scoreName'] != None:
      self.scoreName = opts['scoreName']
    else:
      inString= raw_input("Scoring Function [" + str(self.scoreName) + "]: ")   
      inString= inString.strip()
      if (len(inString) > 0):
        self.scoreName = inString
    
    #Use defaults for parameters
    if opts['maxProducts'] != None:
      self.maxProducts = opts['maxProducts']
    else:
      inString= raw_input("Maximum molecules to synthesize [" + str(self.maxProducts) + "]: ")   
      inString= inString.strip()
      if (len(inString) > 0):
        self.maxProducts= int(inString)

    if opts['maxRuntime'] != None:
      self.maxRuntime = opts['maxRuntime'] * 3600
    else:
      inString= raw_input("Maximum runtime (hours) [" + "%5.4f" % (self.maxRuntime/3600.0) + "]: ")   
      inString= inString.strip()
      if (len(inString) > 0):
        self.maxRuntime= float(inString) * 3600.0

    #inString= raw_input("Maximum product size (as ratio of template size) [" + str(self.maxSizeRatio) + "]: ")   
    #inString= inString.strip()
    #if (len(inString) > 0):
    #  self.maxSizeRatio= float(inString)

    if opts['prefix'] != None:
      inString= opts['prefix']
    else:
      inString= raw_input("Mol id prefix for synthesized mols [" + self.molIDPrefix + "]: ")   
      inString= inString.strip()
    if (len(inString) > 0):
      self.molIDPrefix= inString
    self.idSeq= IDSequence(self.molIDPrefix)

  
  # ---------------------------------------------------------------------------
  """
    This is the money method - returns the input MolCollection extended with all 
    products made during the optimization.
  """
  def optimize(self):

    optimizationTimer= Timer()             # for evaluating overall optimization runtime

    if (self.traceOutFName != None):
      IoUtil.openTraceFile(self.traceOutFName)

    if (self.errOutFName != None):
      IoUtil.openErrorFile(self.errOutFName)

    templateSize= self.templateMol.molObject.NumAtoms()

    if self.scoreName == 'Lingo':
      # Fast Lingo based 2D method
      self.molScoringFunc= LingoScoreFunc(self.templateMolFName)  
    elif self.scoreName == 'ROCS':
      # Slow ROCS based 3D method
      self.molScoringFunc= ROCSScoreFunc(self.templateMolFName)
    else:
      raise Exception('unknown scoring function:' + self.scoreName)

    '''
    # TBD: maybe we can use eval? or something similar to avoid long if/then/elif sequence...
    try:
      self.molScoringFunc = eval(self.scoreName + 'ScoreFunc(self.templateMolFName)') 
    except:
      raise Exception, ('unknown scoring function:' + self.scoreName)
    '''

    # 1. Filter input mols
    #   A. Remove molecules that are nonreactive - and don't count multigroup matching as reactant group matches
    if (not self.molCollection.isPreScored):
      # For now, don't do this on molecules from a prior run...
      IoUtil.logTraceString("Removing nonreactive molecules from input reagents",True)
      reactiveFilter= ReactiveFilter(self.rxnList,False)
      preFilterSize= self.molCollection.size()
      self.molCollection.filter(reactiveFilter)    
      postFilterSize= self.molCollection.size()

      IoUtil.logTraceString("After filtering for reactivity constraints, " + str(postFilterSize) + " of " + str(preFilterSize) + " input mols remain",True)
      if (postFilterSize == 0):
        IoUtil.logErrorString("ERROR: No mols passed reactivity filtering, exiting...")
        sys.exit(1)

    #   B. Build complexity filters filter to avoid products too much larger than the template molecule
    minCarbons=0
    minHeavyAtoms= 1
    self.maxHeavyAtoms2Score= int(templateSize * self.maxSizeRatio2Score)
    self.maxHeavyAtoms2Extend= int(templateSize * self.maxSizeRatio2Extend)
    allowedAtoms= set(['C','H','N','O','P','S','F','Cl','Br','I','Na','K','Ca'])
    forbiddenGroups= None

    simpleEnough2Extend= SyntheticFilter(minCarbons,minHeavyAtoms,self.maxHeavyAtoms2Extend,self.maxRotors2Extend,self.maxChiralAtoms2Extend,
                                        allowedAtoms,forbiddenGroups,self.maxSynSteps)

    simpleEnough2Score= SyntheticFilter(minCarbons,minHeavyAtoms,self.maxHeavyAtoms2Score,self.maxRotors2Score,self.maxChiralAtoms2Score,
                                        allowedAtoms,forbiddenGroups,self.maxSynSteps)

    # 2. Iterate over all starting molecules, assess complexity, and assign product scores:     
    if (self.molCollection.isPreScored):
      highestInd= self.idSeq.getHighestIndex(self.molIDPrefix,self.molCollection.molList)
      self.idSeq.setStartingIndex(highestInd)
      IoUtil.logTraceString("Reagents were prescored, skipping scoring and starting ID sequence at: " + str(highestInd) + "\n",True)
    else:
      IoUtil.logTraceString("Scoring starting materials\n",True)
      molsToScore= self.preScoring(self.molCollection.molList,simpleEnough2Score)
      resultsDict= self.molScoringFunc.values(molsToScore)
      self.postScoring(molsToScore,resultsDict)

    # 3. First block complex molecules from reacting, and analyze the rest for reaction susceptibility
    IoUtil.logTraceString("Analyzing reactions versus reagents\n",True)
    self.tagTerminal(self.molCollection.molList,simpleEnough2Extend)
    for aMol in self.molCollection.molList:
      if (not aMol.isTerminal):
        self.analyzeRxns4Mol(aMol,self.rxnList)

    IoUtil.logTraceString("Reaction Analysis Of Input Mols:\n",True)
    for aMol in self.molCollection.molList:
      IoUtil.logTraceString(aMol.smiStrg + " " + str(len(aMol.myRxns)) + "\n")

    IoUtil.logTraceString("Reaction Analysis Of Input Rxns:\n",True)
    for rxn in self.rxnList:
      IoUtil.logTraceString("Rxn: " + rxn.smirksString + " " + rxn.name  + "\n")
      for i in xrange(len(rxn.rctRgts)):
        IoUtil.logTraceString("  Reactant " + str(i) + " " + str(len(rxn.rctRgts[i])) + "\n")
      IoUtil.logTraceString("\n")

    # Now eliminate any reactions that don't have at least one reagent match per reactant:
    IoUtil.logTraceString("Eliminating dead reactions:\n",True)
    workingRxnList= list()
    for rxn in self.rxnList:
      rxnSatisfied= True
      for i in xrange(rxn.getNumReactants()):
        if (len(rxn.rctRgts[i]) < 1):
          rxnSatisfied= False
      
      if (rxnSatisfied):
        workingRxnList.append(rxn)

    IoUtil.logTraceString("Kept " + str(len(workingRxnList)) + " out of " +  str(len(self.rxnList)) + "\n",True)
    self.rxnList= workingRxnList
    
    # 4. Setup main optimization loop    
    numFruitless= 0
    maxFruitless= 100         # Avoid endless loop if not making anything new
    numProductsTot= 0
    initialTemperature= 1.0
    finalTemperature= 0.5
    temperatureRange= initialTemperature - finalTemperature
    temperature= initialTemperature

    bestMol= self.molCollection.getBestProd()
    bestScore= bestMol.prodScore

    setupTime= elapsedTime= optimizationTimer.getElapsed()
    IoUtil.logTraceString("Elapsed time for setup was " + "%5.4f" % (setupTime/60.) + " minutes\n")

    IoUtil.logTraceString(str(self.molCollection.size()) + " original input molecules assessed, BEST_SCORE= " + "%5.4f" % bestScore + "\n")
    IoUtil.logTraceString("Optimization trace below...\n")
    IoUtil.logTraceString("----------------------------------------------\n")
    IoUtil.logTraceString("NUM_PRODUCTS BEST_SCORE\n")
      
    # 5 BEGIN MAIN OPTIMIZATION LOOP:
    
    bufferedmolsToScore= list()  # buffer molecules for parallel scoring
    
    IoUtil.logTraceString("Entering main optimization loop\n",True)
    while (numProductsTot < self.maxProducts and numFruitless < maxFruitless and elapsedTime < self.maxRuntime):
      # NOTE: This loop gets molecules for improvement using Metropolis sampling 
      # based on product score, and we test inside whether they are truly suitable 
      # for further improvement before exiting
      invalidSM= True
      while (invalidSM):
        # Choose molecule to expand using product score:
        molToImprove= self.molCollection.getRandMetropolisOnProd(temperature)
        invalidSM= (molToImprove.isTerminal or len(molToImprove.myRxns) == 0)
      
      # Okay, expand around this molecule
      products= self.getRandomProducts4Mol(molToImprove,self.maxRgtsPerRct,temperature)
      if (len(products) > 0):
        numFruitless= 0
      else:
        numFruitless= numFruitless + 1
        continue
      
      # Tag complex products as terminal, analyze for reaction susceptibility, add
      # all to molCollection
      self.tagTerminal(products,simpleEnough2Extend)
      for product in products:  
        self.molCollection.add(product)
        if (not product.isTerminal):          # Just flagged as terminal in preScoring, above
          self.analyzeRxns4Mol(product,self.rxnList)

      # Buffer those that are simple enough to score
      molsToScore= self.preScoring(products,simpleEnough2Score)
      for mol in molsToScore:
        bufferedmolsToScore.append(mol)
        
      # If buffer is full, score them, and feed scores back to reagents
      if (len(bufferedmolsToScore) >= self.minToScorePll):
        resultsDict= self.molScoringFunc.values(bufferedmolsToScore)
        self.postScoring(bufferedmolsToScore,resultsDict)
        for product in bufferedmolsToScore:  
          if (product.prodScore > bestScore):
            bestScore= product.prodScore
        
        IoUtil.logTraceString("%6d" % numProductsTot + "  " + "%5.4f" % bestScore + "\n")
        del(bufferedmolsToScore[0:len(bufferedmolsToScore)])   # flush this after use

      
      # Incremental cooling
      temperature= temperature - ((len(products) / float(self.maxProducts)) * temperatureRange)
      if (temperature < finalTemperature):
        temperature= finalTemperature

      numProductsTot= numProductsTot + len(products)
      
      elapsedTime= optimizationTimer.getElapsed()

    # END MAIN OPTIMIZATION LOOP:

    # Score the leftovers:
    IoUtil.logTraceString("Done main optimization loop, scoring leftovers...\n",True)
    if (len(bufferedmolsToScore) > 0):
      resultsDict= self.molScoringFunc.values(bufferedmolsToScore)
      self.postScoring(bufferedmolsToScore,resultsDict)
      for product in bufferedmolsToScore:  
        if (product.prodScore > bestScore):
          bestScore= product.prodScore
      IoUtil.logTraceString("%6d" % numProductsTot + "  " + "%5.4f" % bestScore + "\n")
      del(bufferedmolsToScore[0:len(bufferedmolsToScore)])   # flush this after use


    optimizationTime= optimizationTimer.getElapsed()
    IoUtil.logTraceString("Total Elapsed time for setup, optimization was " + "%5.4f" % (optimizationTime/3600.0) + " hours\n")

    IoUtil.closeTraceFile()
    IoUtil.closeErrorFile()

    return(self.molCollection)
    
  # ---------------------------------------------------------------------------
  """
    Loops over a list of mols, and tags as terminal those that fail the simpleEnough2Extend
    filter passed in.
    
    @param molsToEvaluate: list of mols to be evaluated
    @param simpleEnough2Extend: AbstractMolFilter that returns true if the
      molecule's size/complexity are low enough to permit extension in reactions.
      
  """
  def tagTerminal(self,molsToEvaluate,simpleEnough2Extend):
    for mol in molsToEvaluate:
      if (simpleEnough2Extend.passes(mol)):
        mol.isTerminal= False
      else:
        mol.isTerminal= True

  # ---------------------------------------------------------------------------
  """
    Loops over a list of mols that are candidates for scoring, selecting and returning
    those that are simple enough to score as a list.
    
    @param molsToEvaluate: list of mols to be evaluated
    @param simpleEnough2Score: AbstractMolFilter that returns true if the
      molecule's size/complexity are low enough to permit it to be scored
      
    @returns: a sublist of molecules which can be submitted for scoring
  """
  def preScoring(self,molsToEvaluate,simpleEnough2Score):
    molsToScore= list()

    for mol in molsToEvaluate:
      if (simpleEnough2Score.passes(mol)):
        molsToScore.append(mol)
      else:
        mol.prodScore= -5.0
        
    return(molsToScore)

  # ---------------------------------------------------------------------------
  """
    Loops over a list of mols that were submitted for scoring, and processes 
    their scores as returned from multiple-valued scoring.  Scores are fed
    back to reagents for synthesized molecules.
    
    @param molsSubmitted: list of mols just returned from scoring
    @param resultsDict: a dict of score values that was returned from the 
      scoring call
  """
  def postScoring(self,molsSubmitted,resultsDict):
    strg= str(len(molsSubmitted)) + " submitted  and " + str(len(resultsDict)) + " results returned in postScoring"
    IoUtil.logErrorString(strg)
    
    for mol in molsSubmitted:
      scoreValue= resultsDict.get(mol.name)
      if (scoreValue is None):         # No score was returned for this molecule
        if (mol.prodScore >= -0.5):
          mol.prodScore= -10.0
        mol.isTerminal= True        
        IoUtil.logErrorString("ERROR - no score was returned for mol:\n")
        IoUtil.logErrorMol(mol)
      else:
        mol.prodScore= scoreValue
        if (mol.synOrigin != None):
          for reactant in mol.synOrigin.reactants:
            reactant.updateReagentScore(mol.prodScore)
        
  # ---------------------------------------------------------------------------
  """
    This does 2 things:
      1. Adds an instance variable for a SynthesizedMol - the list of reactions that it
         is susceptible to.
      2. Updates the list of reagents that each matching reaction can work with
  """
  def analyzeRxns4Mol(self,synMol,rxnList):
    if (synMol.isTerminal):                        # Don't register reactivity for terminal mol
      return

    for rxn in rxnList:
      foundMatch4Rxn= False
      for i in xrange(rxn.getNumReactants()):
        numMatches= rxn.numMatchesAsReactant(i,synMol)
        if (numMatches == 1):                   # Only use case of unambiguous match
          rxn.rctRgts[i].add(synMol)            # For the rxn, save all instances of matches for all reactants
          if (not foundMatch4Rxn):
            synMol.myRxns.append([rxn,i])       # For the mol, save it's rxn match, and position, for the 1st match
            foundMatch4Rxn= True  
    
  # ---------------------------------------------------------------------------
  """
    Returns a list of reaction products for the SynthesizedMol passed in, which must have 
    the myRxns list appended.  The reaction products are returned as SynthesizedMol,
    with the name and SyntheticOrigin assigned
      1. A single reaction that synMol can work in is chosen at random
      2. synMol plays a fixed role at the position determined in synMol.myRxns[rxnInd][1]
      3. Other reactants are subsituted by reagents chosen randomly with bias towards
         good prior reagent scores
      4. Products are checked internally for novelty against self.MolCollection
  """
  def getRandomProducts4Mol(self,synMol,maxRgtsPerRct,temperature):
    noProducts= True
    numFruitless= 0
    maxFruitless= 100         # Avoid endless loop if not making anything new
    rxnInd= 0
    rxn= None
    selected= set()
    products= list()    
    
    # To provide traceablility back to the reactants used for each product, first 
    # assemble lists of chosen reagents to use for each reactant, and then make a 
    # bunch of single product libraries using each reagent combination
    while (noProducts and numFruitless < maxFruitless):      # Loop because attempt may be fruitless
      rxnInd= random.randint(0,len(synMol.myRxns)-1)         # Pick a random reaction
      rxn= synMol.myRxns[rxnInd][0]                      
      
      # Now set the reaction up
      libgen= None
      libgen= OELibraryGen(rxn.smirksString)
      libgen.SetExplicitHydrogens(True)
      libgen.SetValenceCorrection(True)
     
      # First assign synMol to the reactant it was originally found as
      synMolRctInd= synMol.myRxns[rxnInd][1]
      
      # rgtsToUse[i][j] will be the jth reagent to use for the ith reactant 
      rgtsToUse= [None] * rxn.getNumReactants() 
      for i in xrange(len(rgtsToUse)):
        rgtsToUse[i]= list()

      rgtsToUse[synMolRctInd].append(synMol)
     
      # Now choose/assign other reactants
      foundAllRgts= True
      for i in xrange(rxn.getNumReactants()):
        if (not i == synMolRctInd):          
          if (len(rxn.rctRgts[i]) < 1):                  # Failure, no reagents for this reactant
            foundAllRgts= False
            break
          elif (len(rxn.rctRgts[i]) <= maxRgtsPerRct):   # Can use all for this reactant
            for rgt in rxn.rctRgts[i]:
              rgtsToUse[i].append(rgt)
          else:                                          # Too many - choose by reagent score
            selected.clear()
            while (len(selected) < maxRgtsPerRct):
              rgt= rxn.rctRgts[i].getRandMetropolisOnRgt(temperature)
              if (not rgt in selected):
                selected.add(rgt)
            
            while (len(selected) > 0):
              rgt= selected.pop()
              rgtsToUse[i].append(rgt)
      
      # Now generate the products as a combinatorial number of 1-product libraries, so that we can
      # assign the synthetic origin of each.
      if (foundAllRgts):
        dimsList= list()
        for i in xrange(len(rgtsToUse)):
          dimsList.append(len(rgtsToUse[i]))
         
        iter= CombiIter(dimsList)
        more= True
        while (more):
          # returns a different index configuration each call, specifying a reagent combination:
          indices= iter.next()
          if (indices == None):
            more= False
            break
          else:
            reactants= list()         # Keep track of reactants for synthetic origin
            for i in xrange(len(rgtsToUse)):
              rgt= rgtsToUse[i][indices[i]]
              reactants.append(rgt)
              rgt.buildMolObject(True)     
              libgen.ClearStartingMaterial(i)
              libgen.SetStartingMaterial(rgt.molObject,i)
              rgt.clearMolObject()              
            
            # There should be only 1 product in each pass
            # NOTE: The outputs of libgen.GetProducts() are NOT OEGraphMol objects
            for product in libgen.GetProducts():
              smiStrg= OECreateCanSmiString(product)   # Check if: product's already been made
              product.Clear()
              del(product)
              
            newSynMol= SynthesizedMol(None,smiStrg,None)
              
            
            if (not (newSynMol in self.molCollection)):
              try:
                newSynMol.assignCanonicalSmiles(True)
                newSynMol.name= self.idSeq.getNextID()
                newSynMol.synOrigin= SyntheticOrigin(reactants,rxn.smirksString)
                products.append(newSynMol)
                noProducts= False
              except StandardError:
                print("Invalid smiles string for product: " + smiStrg + " ignoring...")
        
        if (noProducts):
          numFruitless= numFruitless + 1          
      else:
        numFruitless= numFruitless + 1
             
      
    return(products)          

    
