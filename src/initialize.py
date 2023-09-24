import numpy as np
from constants import EPS, NX, NY

# Random initial bed state
def initialize():
    h0 = EPS * np.random.randn(NY, NX)
    return h0