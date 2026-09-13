"""Memory-bounded nonlinear alternative using training-selected sequence features."""
import numpy as np

def dense_float32(x):
    return x.toarray().astype(np.float32,copy=False)
