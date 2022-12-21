"""
  AbstractMolFilter.py: Abstract base class for molecule filtering functions

  Constructors may be used to initialize the concrete filtering functions with options
"""
class AbstractMolFilter:
  """
    This method must be implemented by all concrete subclasses
    @param basicMol: a BasicMol object to be evaluated
    @returns: True if the molecule passes the filter
  """
  def passes(self,basicMol):
    raise NotImplementedError("value must be defined in subclasses of AbstractMolFilter!")
