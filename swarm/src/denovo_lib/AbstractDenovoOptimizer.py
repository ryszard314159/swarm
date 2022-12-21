"""
  AbstractDenovoOptimizer.py: Base class for optimization strategies for de novo design
  Enables pluggable optimization strategies
"""
class AbstractDenovoOptimizer:
  """
    @param molCollection: MolCollection of molecules (SynthesizedMol) to work with
    @param rxnList: List of reactions (AbstractReactionTransform) to work with
  """
  def __init__(self,molCollection,rxnList):
    self.molCollection= molCollection        # 
    self.rxnList= rxnList                    # 
    self.molScoringFunc= None                # Molecule scoring function (AbstractMolScoreFunc)
  
  # ---------------------------------------------------------------------------
  """ 
    This may be used to do any setup work, including getting info from the user
    thru the cmd line
  """ 
  def setup(self):
    raise NotImplementedError("setup must be defined in subclasses of AbstractDenovoOptimizer!")

  # ---------------------------------------------------------------------------
  """ 
    This is the money method - returns a MolCollection of the best products (or 
    all, since they are scored)
  """ 
  def optimize(self):
    raise NotImplementedError("optimize must be defined in subclasses of AbstractDenovoOptimizer!")
