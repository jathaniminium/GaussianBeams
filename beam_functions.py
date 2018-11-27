from numpy import *
import pickle
from matplotlib import *
from constants import *
from scipy.special import jn
from scipy.misc import factorial as factorial
from scipy.integrate import quad as integrate

#Function to calculate the value of the Laguerre polynomial L_nm(u) of degree n at m=0 for input function u(r,z).
#u is a function of beam radius, r, and distance along optical axis, z.
#p is is the polynomial degree.
def GL_mode(u,p):
	r = len(u[:,0])
	z = len(u[0,:])
	L = zeros((r,z))
	for i in range(p+1):
		L += float(factorial(p))*(-u)**i/( float(factorial(p - i))*(float(factorial(i)))**2. )
	return L

def GL_mode_1dim(u,p):
	L = zeros(1)
	for i in range(int(p+1)):
		L += float(factorial(p))*(-u)**i/( float(factorial(p - i))*(float(factorial(i)))**2. )
	return L

#Calculate integrand for inner product integral.
def integrand(x,s,t,p):
	Lp = GL_mode_1dim(x,p)
	result = Lp * exp(-x/2.) * jn(0,s*t*sqrt(x))
	return result.real

#Function to calculate coupling coefficient for Gauss-Laguerre mode p.  For example, with w_0 = 0.644a, c0 = 0.99.
#r is the beam radius.
#w_0 is the beam waist radius.
#a is half the aperture diameter.
def GL_coupling_coefficient(r, w_0, a, p):
	#Assumes a corrugated horn waveguide, where the E-field at the horn aperture is known analytically (HE_11 mode).
	#With that, used Goldsmith equation 7.34 to calculate coefficients.
	x = 2.*r**2./w_0**2.
	s = w_0/a
	t = 1.7005
	integral = integrate(integrand,0.,2./s**2., args=(s,t,p))[0]
	c_p = 1.362*s*integral
	return c_p
