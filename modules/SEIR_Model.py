# -*- coding: utf-8 -*-
"""
Created on 8/6/2020

@author: bgossage
"""

"""
A simple epidemic model

Mathematical Biology, J.D. Murray
Chapter 19.1

Given a disease in which recovery confers immunity:
 S(t) .=. Susceptibles
 E(t) .=. Exposed
 I(t) .=. Infectives
 R(t) .=. Removed (immune, or isolated)
 D(t) .=. Deaths


(Often called "SIR" models.)

Assumptions:
    1) The rate increase of infectives (r) is proportional to the number of
       infectives and susceptibles
    2) The susceptibles are lost at the same rate r.
    3) The rate of removal of infectives (a) is proportional to the number of Infectives
       that is = a * I
    4) The incubation time is neglibile
    5) The total poplulation is a constant N s.t.  S + E+ I + R = N
    6) dti .=. incubation period

The ODE system is:

   dS/dt = -r * S(t) * I(t)
   dE/dt = r * S(t) * I(t)
   dI/dt = c * E(t) - (a+d) * I(t)
   dR/dt = (a-d) * I(t)
   dD/dt = d * I(t)
"""

import numpy
import scipy.integrate

## A simple epidemic model
#
#     Mathematical Biology, J.D. Murray
#     Chapter 19.1.
#
class SEIR_Model:
## The constructor.
   def __init__(self, name="", sequence=""):

   # Parameters
      self.r = 0.4   ## The infection rate
      self.a = 0.01  ## The removal rate
      self.d = 0.01  ## The death rate

   # Initial values
      self.S0 = 1.0 - self.R0 - self.I0
      self.I0 = 0.0001  ##
      self.R0 = 0.0

      self.D0 = 0.0

   # Time...
      self.tStart = 0.0
      self.tStop =  1.0
      self.tInc = 0.01
      
      self.solution = None

   #end SEIR_Model.__init__ !!!!!!!!!!!!!!!!!!!!!!!!!!!!

## Run the model
   
   def run( self, tstart, tend, dt ):
       
   # Bundle parameters for ODE solver
      params = [ self.r, self.a, self.d ]
      
   # Bundle initial conditions for ODE solver
      y0 = [ self.S0, self.I0, self.R0, self.D0 ]
      
   # Make a time array for solution samples...
      t_0 = numpy.arange( self.tStart, self.tStop, self.tInc)
       
   # Solve using numerical integration...
      self.solution = scipy.integrate.odeint( deriv, y0, t_0, args=(params,) )
   
   #end SEIR_Model.run ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   
# end class SEIR_Model
   
   
# EOF
