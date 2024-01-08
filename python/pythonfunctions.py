import numpy as np

def gaussian(x, a, b, c):
    z = (x - b) / c;
    return a * np.exp(-0.5 * z * z);

def sro_model(T, a0, a1, a2,):
    return (np.exp(-a0*(T**(-1))) - 1)*(a1 + a2*(T**(-1)))
