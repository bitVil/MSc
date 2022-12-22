import numpy as np
from scipy.optimize import root
import pandas as pd


# # Masker helper class
#
# push to separate library once the functionality has settled down.
# The intention of this class is to separate the logic of projecting
# down to a subspace or to lift out of a subspace.

class Masker:
    """
    Helper class to quickly project on a subspace given 
    by a set of indices
    """
    def __init__(self, x0, indices):
        """
        `indices` can be list of label names or list of integer numbers.
        In the case of label names, we need `x0` to be a named tuple, so 
        that we can map to actual indices.  Internally we always store as 
        np.arrays for efficiency.
        """

        if len(indices) > 0 and isinstance(indices[0], str):
            indices = list( map(x0._fields.index, indices) )
            self.namedtuple = x0.__class__
        
        self.x0 = np.array(x0)
        
        self.indices = np.array(indices,dtype=np.int)
        self.x0_mask = self.x0[self.indices]

    def mask(self, x):
        """
        project to subspace
        """
        return x[self.indices]
    
    def _unmask(self, x):
        """
        internal unmask, do not convert to 
        namedtuple for efficiency
        """
        x0 = np.copy(self.x0)
        x0[self.indices] = x
        return x0
    def unmask(self, x):
        """
        lift from subspace
        """
        x0 = self._unmask(x)
        if self.namedtuple:
            x0 = self.namedtuple(*x0)
        return x0
    
    def f_unmask(self, f):
        """
        unmask first argument then apply f
        """
        return lambda x: f( self.unmask(x) )
        


# # Root finding wrapper

def root_indices(f,x0,indices):
    """
    wrapper for scipy.optimize.root to work with indices
    """    
    mask = Masker(x0, indices)    
    myroot = root( mask.f_unmask(f), mask.x0_mask )
    myroot["x"] = mask.unmask(myroot.x)
    
    return myroot 


# # Path following algorithm
#
# use `yield` instead of `return` to feed into `for` loops and `pandas`

# +
def follow_path(f, x0, t, epsilon, maxN = 100):  
    
    xg = x0 + t/np.linalg.norm(t)*epsilon
    for i in range(maxN):
        def tracer(x):
            return np.append(f(x), np.linalg.norm(x-x0)-epsilon)

        sol = root(tracer, xg)
        if not sol.success:
            print(f"did not converge after {i} steps at {xg}")
            break
         
        yield sol.x
        xg = 2*sol.x - x0
        x0 = sol.x

def follow_path_indices(f, x0, indices, *args, **kwargs ):
    mask = Masker(x0, indices)
    for i in follow_path(mask.f_unmask(f), mask.x0_mask, *args, **kwargs):
        yield mask.unmask(i)

def follow_path_indices_up_down(f, x0 ,indices, eps, Nmax):
    """
    follow f(x)==0 with x0 as starting guess and projected to indices.
    
    """
    x0 = root_indices(f, x0, indices[1:]).x
    t = np.zeros(len(indices))
    t[0] = 1.0

    result = pd.DataFrame(
        follow_path_indices(f, x0, indices, -t, eps, Nmax)
    )
    result = result[::-1].drop(0)
    result.index = -result.index
    
    result = result.append( 
        pd.DataFrame(
            follow_path_indices(f, x0, indices, t, eps, Nmax)
        ) 
    )
    return result



# -
