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
    5) The total poplulation is a constant N s.t.  S + I + R = N
    6) dti .=. incubation period

The ODE system is:

   dS/dt = -r * S(t) * I(t)
   dI/dt = r * S(t) * I(t) - (a+d) * I(t)
   dR/dt = (a-d) * I(t)
   dD/dt = d * I(t)
"""
class SEIR_Model:

   def __init__(self, name="", sequence=""):

   # Parameters
      self.r = 0.4   # The infection rate
      self.a = 0.01  # The removal rate

   # Initial values
      self.I0 = 0.0001
      self.R0 = 0.0
      self.S0 = 1.0 - self.R0 - self.I0
      self.D0 = 0.0

   # Time...
      self.tStart = 0.0
      self.tStop =  0.0

   #end SEIR_Model.__init__

