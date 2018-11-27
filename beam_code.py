#Title: beam_code.py
#Author: Jason W. Henning 2011/5/6
#Purpose: Calculate horn beam profile using a specified number of Gauss-Laguerre orthonormal modes.
#Directions: Specify frequency (f), horn radius (a), beam waist radius (w_0), number of modes (p), size 
#            of profile box in radius r and optical axis z (length), and spatial resolution (resolution).
#            w_0 is chosen to maximize the size of the p=0 mode coupling coefficient.
#
#Required modules: beam_functions.py, constants.py, numpy, scipy, matplotlib, pickle, time
#Revisions: 2011/5/8 - (JWH) Removed lots of unnecessary for loops which sped up code.


from numpy import *
from matplotlib import *
from constants import *
from beam_functions import *
import pickle as pk
from time import *
import pylab as py

t_start = clock()
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#Parameter definition block.
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#Frequency (in GHz).
f = 145.

#Horn radius (in mm).
a = 0.39*25.4/2.

#Horn opening angle (in degrees).
theta_h = 12.6

#Beam radius at aperture.  This maximizes the power in the primary gaussian mode for
#corrugated cylindrical horns.
w_ap = 0.644*a #in mm. For a corrugated cylindrical waveguide.


#Number of Gauss-Laguerre modes to use.
p = 30 #45 #Integer

#Size of r and z in mm
length = 100

#Resolution, in mm.
resolution = 0.1
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

#--------------------------------------------------------------------------------------------------------
#Main program body: start
#--------------------------------------------------------------------------------------------------------
#wavelength
Lambda = c/(f*10.**9.) * 10. #in mm

#wavenumber
k = 2.*pi/Lambda #in 1/mm

#Get radius of curvature at aperture.
R_ap = a/sin(theta_h*pi/180.)

#Calculate w_0.
w_0 = w_ap / sqrt(1. + (pi**2. * w_ap**4.)/(Lambda**2. * R_ap**2.))

#Confocal distance
z_c = (pi * w_0**2.)/Lambda

#Get z_ap (the distance from the beam waist to the aperture).
z_ap = z_c * (pi*w_ap**2./(Lambda*R_ap))

#Distance from beam origin (beam waist, NOT aperture).
z = array(arange(0,length, resolution))

#Radial distance from center of beam.
r = array(arange(0,length, resolution))

#Radius of curvature
Rc = array(z + z_c**2./z)

#Beam radius
w = array(w_0* sqrt(1. + (z/z_c)**2.))

#Phase shift
phi_0 = array(arctan(z/z_c))

#Define a 'u' variable for calculating Laguerre polynomials.
u = zeros((len(r),len(z)))
for i in range(len(r)):
	u[i,:] = 2.*r[i]**2./w**2.

#First get the Laguerre polynomials of u and mode p.  We're assuming axisymmetric horns, so m=0.
L = []
p_var = range(p+1)
for i in range(p+1):
	L.append(GL_mode(u,p_var[i]))

#Now define the E-fields for each mode.  These form an orthonormal set of functions that can be used to
#decompose the E-field at any distance along the optical axis from the horn aperture.
E = []
for n in range(p+1):
	E_temp = zeros((len(r),len(z)), 'complex')
	for i in range(len(r)):
			#E-field for m=0 (axially symmetric modes) in cylindrical coordinates
			#by Goldsmith equation 2.53.
		E_temp[i,:] = sqrt(2./(pi*w**2.)) * \
			      L[n][i,:] * \
                              exp((-r[i]**2./w**2.) -1j*k*z - (1j*pi*r[i]**2./(Lambda*Rc)) + 1j*(2.*n + 1.)*phi_0)
	E.append(E_temp)

#Calculate coupling coefficients for each mode.  I.e., get each c_p.  This is done by using the known
#E-field illumination at the horn aperture as a boundary condition to determine how much of each
#orthonormal mode is required to reproduce the E-field illumination pattern.
c_p = []
for n in range(p+1):
	q = GL_coupling_coefficient(r, w_0, a, n)
	c_p.append(q)

c_tot = zeros(1)
for i in range(p+1):
	c_tot += c_p[i]**2.

#Sum up the E field modes to get the total E-field illumination pattern.
E_tot = zeros((len(r),len(z)), 'complex')
for n in range(p+1):
	E_tot += c_p[n]*E[n]
	
##Calculate final power.  Simply |E|^2.
P_tot = zeros((len(r),len(z)))
P_tot += (E_tot * E_tot.conj()).real

#Find maximum of P_tot at each z-value, and divide so beam is normalized to r=0 power at all z along optical axis.
for l in range(len(z)):
	P_tot[:,l] /= max(P_tot[:,l])

t_end = clock()
print 'Elapsed time:', t_end - t_start
#---------------------------------------------------------------------------------------------------------------------------
#Main program body: End
#---------------------------------------------------------------------------------------------------------------------------


#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#Save block
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#Save the variables in a pickle file since this takes forever to run.
#File = open(str(f)+'GHz_beam_'+str(p)+'modes_'+str(length)+'mm_res'+str(resolution)+'mm.dat', 'wb')
#pk.dump(f,File)
#pk.dump(a,File)
#pk.dump(p,File)
#pk.dump(k,File)
#pk.dump(w_0,File)
#pk.dump(z_c,File)
#pk.dump(z, File)
#pk.dump(r, File)
#pk.dump(Rc, File)
#pk.dump(w, File)
#pk.dump(phi_0, File)
#pk.dump(u, File)
#pk.dump(L, File)
#pk.dump(E, File)
#pk.dump(c_p, File)
#pk.dump(E_tot, File)
#pk.dump(P_tot, File)
#File.close()
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
#///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Plotting block
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#Make a beam profile plot in the z-r plane.
py.ion()
py.clf()
py.imshow(P_tot, norm = colors.LogNorm(vmin=1e-6, vmax=1.), interpolation='nearest')
ax = py.gca()
py.ylim(ax.get_ylim()[::-1]) #Reverse y-axis
#py.xlim((-floor(z_ap - z_c),50))
py.hot()
#py.colorbar(ticks = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1.])
#py.contour(P_tot, [1e-1,1e-2,1e-3,1e-4,1e-5,1e-6], colors='b')
#py.plot([0./resolution,5./resolution],[2.8/resolution,2.8/resolution], 'k')
#py.plot([2.134/resolution,2.134/resolution],[0./resolution,5./resolution], 'k')
locs = zeros(10)
labels = zeros(10)
for i in range(10):
	locs[i] = i*length/resolution/10.
	labels[i] = i*length/10.	
#locs = range(0,int(length/resolution),int(length/(resolution*10.)))
#labels = range(0,int(length),int(length/len(locs)))
py.yticks(locs, labels)
py.xticks(locs, labels)
py.xlabel('Distance from horn aperture (mm)')
py.ylabel('Beam radius (mm)')
py.savefig('SPTpol_150GHz_beam.pdf')
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
