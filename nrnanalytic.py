# -*- coding: latin-1 -*-
"""
Find numerical values for the analytical expressions presented in:

Monai H, Omori T, Okada M, Inoue M, Miyakawa H, Aonishi T. 
"An analytic solution of the cable equation predicts frequency 
preference of a passive shunt-end cylindrical cable in response 
to extracellular oscillating electric fields."
Biophys J. 2010 Feb 17;98(4):524-33.

and

W. Ying and C. S. Henriquez, "Hybrid finite element method 
for describing the electrical response of biological cells
to applied fields," IEEE Transactions on Biomedical Engineering, 
vol. 54, no. 4, pp. 611-620, Apr. 2007. 

and

Tadej Kotnik, Damijan Miklavcic, Tomaz Slivnik, Time course of 
transmembrane voltage induced by time-varying electric fields 
a method for theoretical analysis and its application, 
Bioelectrochemistry and Bioenergetics, Volume 45, Issue 1, 
March 1998.

and

Bedard C, Kroger H and Destexhe A.  Modeling extracellular field
potentials and the frequency-filtering properties of extracellular
space. Biophysical Journal 86: 1829-1842, 2004.

Corresponding to the response of a cable to a DC step current 
input, the cross section of a cylindrical cell to a homogeneous
field, a sphere to a homogeneous field, and the potential 
for a point source current.

Typical parameters are specified as global variables for simplicity
in the Python module and can be changed if necessary before 
function invocation.

Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)

"""
import numpy as np
from numpy import pi, exp, cos, sqrt

#Typical parameters
d = 1.0e-4 #Cable or cell diameter (cm)
D = 10000.0e-4 #Extracellular space diameter (cm)
L = 80e-4 #Cable length (cm)
Rm = 1000.0 #Specific membrane resistance (Ohm*cm^2)
Cm = 1.0 #Specific membrane capacitance (uF/cm^2)
sigmae = 10.0 #Extracellular conductivity (mS/cm)
sigmai = 10.0 #Intracellular conductivity (mS/cm)
E = 10000.0 #Stimulus electric field (mV/cm)
I = 1000000.0 #Point source current (uA)

def Ri(val):
  """
  Set intracellular resistivity
  """
  global sigmai
  sigmai = 1000.0/val

def constants():
  """
  Return the classical (intracellularly only) time and
  length constants
  """
  Ri = 1000.0/sigmai #Specific intracellular resistivity (Ohm*cm)
  ri = 4.0*Ri/(pi*d**2) #Intracellular resistance per unit length (Ohm/cm)
  rm = Rm/(pi*d) #Membrane resistance per unit length (Ohm*cm^2/cm)
  cm = Cm*pi*d #Membrane capacitance per unit length (uF/cm^2*cm)

  print 'Rm',Rm,'Ri', Ri, 'ri', ri, 'rm', rm, 'cm', cm  
  
  taum = 1e-3*rm*cm #Membrane time constant (ms)
  tauc = 1e-3*Cm*L*Ri #Cell time constant (ms)
  lambd = sqrt(rm/ri) #Length constant (cm)
  
  
  return taum, tauc, lambd
  
def cable(x, time, niter=10000):
  """
  Analytical solution of the cable equation for a homogeneous
  electric field E and with membrane voltage measured at point x.
  Here d is the diameter of the cable
  """
  Re = 1000.0/sigmae #Specific extracellular resistivity (Ohm*cm)
  Ri = 1000.0/sigmai #Specific intracellular resistivity (Ohm*cm)
  ri = 4.0*Ri/(pi*d**2) #Intracellular resistance per unit length (Ohm/cm)
  re = 4.0*Re/(pi*D**2) #Extracellular resistance per unit length (Ohm/cm)
  rm = Rm/(pi*d) #Membrane resistance per unit length (Ohm*cm^2/cm)
  cm = Cm*pi*d #Membrane capacitance per unit length (uF/cm^2*cm)
  lambd = sqrt(rm/(re + ri)) #Length constant (cm)
  tau = 1e-3*rm*cm #Time constant (ms)
  print 'Check why Assaf used Ic = E/re'
  Ic = E/re #Stimulus current (mA) (Not sure why re, 
            #see Plonsey Bioelectricity book sect. 6.1)
            #It might also be re/(ri+re)
  #print re,ri
  
  phi = lambda mu, x: cos(mu*x) #(unitless, not potential!)
  alpha = L/2.0 #(cm)
  
  Vm = np.zeros(np.size(time))
  dVmdt = np.zeros(np.size(time))
  for n in range(niter):
      mu = n*pi/L  #(1/cm)
      kappa = tau/(1.0 + mu**2*lambd**2) #(ms)
      A = re*lambd**2/(tau*alpha)*(phi(mu, L) - phi(mu, 0.0))*phi(mu, x) #(Ohm/ms)
      Vm = Vm + A*Ic*kappa*(1.0 - exp(-time/kappa)) #(Ohm*mA = mV)
      dVmdt = dVmdt + A*Ic*kappa*(exp(-time/kappa)/kappa)
      print kappa, A

  #Find membrane current from capacitive and ionic currents
  Ic = Cm*dVmdt  
  Iion = 1000*Vm/Rm
  Im = Cm*dVmdt + Iion
  
  return Vm, Im, Ic, Iion

def cylinder(r, theta, time):
  """
  time is a vector or float
  Either theta or r can be a vector but not at the same time
  """
  assert(type(r) is np.ndarray or type(theta) is np.ndarray)
  
  n = max(np.size(r),np.size(theta))
  
  ettauip = np.zeros(np.size(time))
  a = np.zeros(np.size(time))
  b = np.zeros(np.size(time))
    
  tauip = 1.0/(1.0/(Cm*Rm) + (2.0*sigmai*sigmae)/(Cm*d*(sigmai + sigmae)))
  ettauip[:] = exp(-time/tauip)
  epsilon = tauip/(Cm*Rm)
  a[:] = ((2.0*sigmae)/(sigmai + sigmae))*(ettauip + (1 - ettauip)*epsilon)
  b[:] = 1.0 - ((2.0*sigmai)/(sigmai + sigmae))*(ettauip + (1 - ettauip)*epsilon)
  
  phie = np.zeros((np.size(time), n))
  phii = np.zeros((np.size(time), n))
  Vm = np.zeros((np.size(time), n))

  for i in range(0, np.size(time)):
    phie[i,:] = -E*r*cos(theta) - b[i]*E*(d**2/(4.0*r))*cos(theta)
    phii[i,:] = -a[i]*E*r*cos(theta)
    Vm[i,:] = E*d*cos(theta)*(1.0 - ettauip[i])*(1.0 - epsilon)

  return phie, phii, Vm

def sphere(theta, time):
  """
  Membrane potential for a sphere under a constant field
  """
  ettauip = np.zeros(np.size(time))
  tauip = 1.0/(1.0/(Cm*Rm) + (4.0*sigmai*sigmae)/(Cm*d*(sigmai + 2.0*sigmae)))
  ettauip[:] = exp(-time/tauip)
  Vm = np.zeros((np.size(time), np.size(theta)))
  for i in range(0, np.size(time)):
    Vm[i,:] = 3.0/2.0*E*d/2.0*cos(theta)*(1.0 - ettauip[i])
  
  return Vm

def totalcurrent(J):
  """  
  Calculate total current (uA) given a current density (uA/cm2) and a 
  sphere of diameter d (module variable)
  """
  #Radius
  a = d/2.0
  I = J*(4*np.pi*a**2)
  return I

def source(r):
  """
  Provides the potential (mV) and current density at d/2 
  (uA/cm2) for a point source current at distance r. This 
  can be used to predict the potential decay for a constant
  current source within a spherical cell.
  After a long time the amount of current coming out of the cell 
  should be homogeneous (radially) and should acount for the 
  total current being injected.
  It is assumed the conductivity of the medium is given by sigmae.
  Also provides electric field (mV/cm)
  """

  #Radius
  a = d/2.0
  
  #Calculate current density at the cell surface (uA/cm2)
  J = I/(4*np.pi*a**2)

  #Calculate potential at r, r can be a numpy array
  phi = (a**2*J/sigmae)*(1.0/r)
  
  #Calculate electric field at r, r can be a numpy array
  E = -(a**2*J/sigmae)*(1.0/r**2)
  
  return phi, E, J

def line(time, Phi0, Phi1):  
  """
  Analytic solution for the complete cell membrane problem in 1D.
  A constant everywhere conductivity is used (sigmae)
  D must be set properly as the domain size
  Phi0 and Phi1 are the dirichlet values at the extremes
  TODO: use E instead of Phi0, Phi1 for consistence with other
  solutions
  Gives the result for Vm in the right hand side of the 1D cell.
  """
  
  gl = 1e3*(1/Rm) #mS/cm2
  PhiD = Phi1 - Phi0
  
  a = (2.0*sigmae)/(Cm*D) + gl/Cm
  b = -(sigmae*PhiD)/(Cm*D)
  
  Vm = -(b/a)*exp(-a*time) + b/a
  Ic = Cm*b*exp(-a*time)
  Iion = gl*Vm
  Im = Ic + Iion
  
  return Vm, Ic, Im, Iion

def lineit():  
  """
  Iterative euler solution for the complete cell membrane problem in 1D.
  A constant everywhere conductivity is used (sigmae)  
  """
  
  #alpha = 1 - 2.0*  
  
  return 0

if __name__=='__main__':
  """
  Run some tests
  """

  #Monai 2010 Parameters
  #d = 1.2e-4 #Cable diameter (cm)
  #D = 10e-4 #Extracellular space diameter (cm)
  #L = 700e-4 #Cable length (cm)
  #Rm = 30000.0 #Specific membrane resistivity (Ohm*cm^2)
  #Cm = 1.5 #Specific membrane capacitance (uF/cm^2)
  #Re = 20.0 #Specific extracellular resistivity (Ohm*cm)
  #Ri = 200.0 #Specific intracellular resistivity (Ohm*cm)
  ##I ~ -5e-9 #Stimulus current (mA)
 
  #Use typical parameters to compare with well known
  #NEURON results (Final voltage at the tip ~39.168 mV)
 
#  t = np.arange(0.0,0.5,0.001)
#  x = 0.0
#  #D = 100*D
#  Vm,_,_,_ = cable(x,t) 
#  figure(1);clf()
#  plot(t,Vm)
#  print Vm[-1]

#  # Use typical parameters to and see effect of number
#  # of iterations
#  taum, tauc, lambd = constants()
#  L = lambd*0.1
#  t = np.arange(0.0,0.5,0.001)
#  x = 0.0
#  Vm0,_,_,_ = cable(x,t,1)
#  Vm1,_,_,_ = cable(x,t,2)
#  Vm2,_,_,_ = cable(x,t,3)
#  Vm3,_,_,_ = cable(x,t,4)
#  #Vm4,_,_,_ = cable(x,t,25)
#  #Vm5,_,_,_ = cable(x,t,26)
#  #VmInf,_,_,_ = cable(x,t,6)
#  figure(1);clf()
#  plot(t,Vm0,label='Vm0')
#  plot(t,Vm1,'-x',label='Vm1',markevery=15)
#  plot(t,Vm2,label='Vm2')
#  plot(t,Vm3,'-x',label='Vm3',markevery=15)
#  plot(t,Vm4,label='Vm4')
#  #plot(t,VmInf)
#  legend()

  #Compare to Ying 2007 figure 4, time = 0.02 us
#  d = 15e-4 #Cell diameter (cm)
#  Rm = 1000.0 #Specific membrane resistivity (Ohm*cm^2)
#  Cm = 1.0 #Specific membrane capacitance (uF/cm^2)
#  sigmae = 20.0 #Extracellular conductivity (mS/cm)
#  sigmai = 5.0 #Intracellular conductivity (mS/cm)
#  E = 100.0 #Stimulus electric field (mV/cm)
#    
#  r = d/2
#  theta = np.arange(0, 2*pi, 15*(pi/180))
#  t = 0.00002
#  
#  phie, phii, Vm = cylinder(r, theta, t)
#
#  figure(2);clf()
#  plot(theta, phie[0,], '-x')
#  plot(theta, phii[0,], '-*')
#  plot(theta, Vm[0,], '-+')
  
  #Produce some informative plots about the Rm/Ri relation
  #in determining the length constant
  #for different diameters
#  ds = [1.0e-4,10.0e-4]  
#  Rms = arange(1000.0,40000.0,1000.0)
#  Ris = [50.0,100.0,200.0]
#
#  figure(3);clf();grid()
#  for d0 in ds:
#    d = d0
#    for Ri0 in Ris:
#      Ri(Ri0)
#      lambd = empty_like(Rms)
#      for i in range(len(Rms)):
#        Rm = Rms[i]
#        _,_,lambd[i] = constants()
#      plot(Rms,lambd,'k')
#      #plot(Rms,sqrt((d*Rms)/(4.0*Ri0))) #DEBUG
#      text(Rms[-1]*1.02,lambd[-1],u'%g ohm*cm,%gµm'%(Ri0,d0*1e4))
#  tight_layout()
#  xlabel('Rm [ohm*cm^2]')
#  ylabel('Length const. [cm]')
#  xlim(0.0,52000.0);ylim(0.0,0.5)


  
  