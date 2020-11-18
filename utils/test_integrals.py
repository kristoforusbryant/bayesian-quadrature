import numpy as np 
### GENZ FUNCTIONS 
# 1. CONTINUOUS INTEGRAND FAMILY (cont) 
def genz1(a: np.ndarray, u: np.ndarray, x: np.ndarray)->float:
    return np.exp(- np.sum(a * np.abs(x-u))) 

# 2. CORNER PEAK INTEGRAND FAMILY (copeak)
def genz2(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.power(1 + np.sum(a*x), -len(a) - 1) 

# 3. DISCONTINUOUS INTEGRAND FAMILY (disc) 
def genz3(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    if all((x - u) > 0): return 0 
    else: 
        return np.exp(- np.sum(a * x))

# 4. GAUSSIAN PEAK INTEGRAND FAMILY (gaussian)
def genz4(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.exp(-np.sum(np.power(a,2) * np.power(x-u,2)))

# 5. OSCILLATORY INTEGRAND FAMILY (oscil)
def genz5(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.cos(2 * np.pi * u[0] + np.sum(a * x))

# 6. PRODUCT PEAK INTEGRAND FAMILY (prpeak)
def genz6(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.prod(1/(1/np.power(a, 2) + np.power(x-u, 2)))

# 7. STEP FUNCTION
def step(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float:
    return 1 if (1 >= x[0] and x[0] > .5) else 0 

genz = {1: genz1,
        2: genz2,
        3: genz3,
        4: genz4,
        5: genz5,
        6: genz6}