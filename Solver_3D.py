'''
Created on Sep 29, 2016

@author: Keisuke
'''
from numpy import *
from pylab import *


class Solver_3D:
    def __init__(self,X):
        self.X=X

def solver_3D(a,b,c,d):
    

    X0=-b/(3*a)\
        -((2**(1/3)*(-b**2+3*a*c))\
        /(3*a*(-2*b**3+9*a*b*c-27*a**2*d+(4*(-b**2+3*a*c)**3+(-2*b**3+9*a*b*c-27*a**2*d)**2)**(1/2))**(1/3)))\
        +((-2*b**3+9*a*b*c-27*a**2*d+(4*(-b**2+3*a*c)**3+(-2*b**3+9*a*b*c-27*a**2*d)**2)**(1/2))**(1/3)\
        /(3*2**(1/3)*a))
    
    X=real(X0)
    X1=Solver_3D(X)
    return(X1)