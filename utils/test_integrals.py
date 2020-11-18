import numpy as np 
from scipy.stats import norm
from itertools import combinations
from math import factorial

### GENZ FUNCTIONS 
# 1. CONTINUOUS INTEGRAND FAMILY (cont) 
def genz1(a: np.ndarray, u: np.ndarray, x: np.ndarray)->float:
    return np.exp(- np.sum(a * np.abs(x-u))) 

def genz1_int(a: np.ndarray, u: np.ndarray)->float:
    prod = 1
    for i in range(len(a)): 
        prod *= (1/a[i])*(2-np.exp(a[i] * u[i]) - np.exp(a[i]*(u[i] - 1)))
    return prod 

# 2. CORNER PEAK INTEGRAND FAMILY (copeak)
def genz2(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.power(1 + np.sum(a*x), -len(a) - 1) 

def genz2_int(a: np.ndarray, u: np.ndarray)->float:
    d = len(a)
    acc = 0 
    for k in range(d):
        for I in combinations(list(range(d)),k+1): 
            acc += np.power(-1, k+d) / (1 + np.sum(a) - np.sum(a[list(I)]))
    return acc / (factorial(d) * np.prod(a))

# 3. DISCONTINUOUS INTEGRAND FAMILY (disc) 
def genz3(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    if all((x - u) > 0): return 0 
    else: 
        return np.exp(- np.sum(a * x))
    
def genz3_int(a: np.ndarray, u: np.ndarray)->float:
    if len(a) == 1: 
        return (1/a[0]) * (np.exp(a[0] * u[0]) -1)
    elif len(a) > 1: 
        prod = 1
        for i in range(len(a)): 
            prod *= (1/a[i]) * (np.exp( a[i] * np.min((1, u[i])) ) - 1)
        return prod 
    
# 4. GAUSSIAN PEAK INTEGRAND FAMILY (gaussian)
def genz4(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.exp(-np.sum(np.power(a,2) * np.power(x-u,2)))

def genz4_int(a: np.ndarray, u: np.ndarray)->float:
    prod = np.power(np.pi, len(a)/2)
    for i in range(len(a)): 
        prod *= (1/a[i]) * (norm.cdf(np.sqrt(2) * a[i] * (1-u[i])) - norm.cdf(-np.sqrt(2) * a[i] * u[i]))
    return prod 

# 5. OSCILLATORY INTEGRAND FAMILY (oscil)
def genz5(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.cos(2 * np.pi * u[0] + np.sum(a * x))

h = {1: np.sin, 2: lambda x: -np.cos(x), 3: lambda x: -np.sin(x), 0: np.cos}
def genz5_int(a: np.ndarray, u: np.ndarray)->float:
    d = len(a)
    acc = 0 
    for k in range(d):
        for I in combinations(list(range(d)),k+1): 
            acc += np.power(-1, k) * h[np.mod(d, 4)](2*np.pi*u[0] + np.sum(a) - np.sum(a[list(I)]))
    return acc / np.prod(a)    
    
# 6. PRODUCT PEAK INTEGRAND FAMILY (prpeak)
def genz6(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float: 
    return np.prod(1/(1/np.power(a, 2) + np.power(x-u, 2)))

def genz6_int(a: np.ndarray, u: np.ndarray)->float:
    prod = 1
    for i in range(len(a)): 
        prod *= (a[i]) * (np.arctan(a[i]*(1-u[i])) - np.arctan(-a[i] * u[i]))
    return prod

# 7. STEP FUNCTION
def step(a: np.ndarray, u: np.ndarray, x:np.ndarray)->float:
    return 1 if (1 >= x[0] and x[0] > .5) else 0 

def step_int(a: np.ndarray, u: np.ndarray)->float:
    return .5

f_list = {1: genz1,
        2: genz2,
        3: genz3,
        4: genz4,
        5: genz5,
        6: genz6, 
        7: step}

int_list = {1: genz1_int,
            2: genz2_int,
            3: genz3_int,
            4: genz4_int,
            5: genz5_int,
            6: genz6_int,
            7: step_int}

default_a = {1: lambda d: 150/np.power(d,3),
             2: lambda d: 600/np.power(d,3), 
             3: lambda d: 10/np.power(d,3), 
             4: lambda d: 100/np.power(d,2), 
             5: lambda d: 110/np.power(d, 5/2), 
             6: lambda d: 600/np.power(d,3), 
             7: lambda d: 0}

class IntegralFunctions: 
    def __init__(self, d: int, u=None, a=None):
        self.d = d 
        self.u = np.array([.5 for _ in range(d)]) if not u else u  
        self.a = {i: np.array([default_a[i](d) for _ in range(d)]) for i in range(1,8)} if not a else a
        self.f_list = f_list 
        self.integrals = [int_list[i](self.a[i], self.u) for i in int_list]
        
    def __call__(self, x: np.ndarray, i: int): 
        return self.f_list[i](self.a[i],self.u,x)
   