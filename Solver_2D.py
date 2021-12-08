'''
Created on Sep 29, 2016

@author: Keisuke
'''
from numpy import *
from pylab import *

class Solver_2D:
    def __init__(self,X):
        self.X=X

def solver_2D(a,b,c):
    

    X0=(-b+(b**2-4*a*c)**(1/2))/(2*a)
    
    X=real(X0)
    X1=Solver_2D(X)
    return(X1)