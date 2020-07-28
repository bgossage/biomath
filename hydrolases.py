# -*- coding: utf-8 -*-
"""
Created on Mon Jul 24 17:30:13 2017

@author: bgossage
"""

"""
Enzyme Kinetics (Hydrolases)

E is the enzyme and S the substrate.  ES is the enzyme-substrate complex.
P1 is the product

Chemical equations:
E + S -> ES 
E + S <- ES
ES -> E + P1

In terms of concenrations [ES] = c, [E] = e0 - c, [S] = s;
where e0 is the initial enzyme concentration.

Assuming mass-action (The reaction rate is proportional to the concentrations.),
the differential equations are:

ds/dt = -kp1 * [(e0 - c) * s)] + km1 * c 

dc/dt = -(km1+kp2) * c + kp1 * [(e0 - c) * s)]

where:
   kp1, km1, and kp1 are the mass-action rate constants.


 Below is the solution in python using numerical integration. 
 
"""

import numpy
import scipy.integrate

import matplotlib

matplotlib.use('TkAgg')  ## set the back end (a wart)
import matplotlib.pyplot

# Parameters
kp1 = 2.0  # The E + S to ES reaction rate
kp2 = 1.9  # The ES to E + P1 reaction rate
km1 = 1.5  # The ES to E + S reaction rate
e0 = 2.0   # The initial enzyme concentration

#
# Define a function that computes the derivatives at time t
# given the current concentrations s and c
#
def deriv( y, t, params ):

    s, c = y      # unpack current values

    kp1, kp2, km1, e0 = params  # unpack parameters
    
    derivs = [ -kp1 * (e0 - c) * s + km1 * c,
               -(km1+kp2) * c + kp1 * (e0 - c) * s  ]
    
    return derivs

# end function derivs ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Initial values
C0 = 0.0
S0 = 10.0

# Bundle parameters for ODE solver
params = [ kp1, kp2, km1 , e0 ]

# Bundle initial conditions for ODE solver
y0 = [ S0, C0 ]

# Make a time array for solution samples...
tStop = 10.0
tInc = tStop / 300.0
t = numpy.arange(0., tStop, tInc)

# Solve using numerical integration...
psoln = scipy.integrate.odeint( deriv, y0, t, args=(params,) )

Ssoln = psoln[:,0]
Csoln = psoln[:,1]

# Estimate the max value of c(t)
km = (km1 + kp2) / kp1
Cmax = e0 * S0 / (S0 + km) # Assuming S0 >> e0

print("Cmax = ", Cmax )

# Plot the solution...
matplotlib.pyplot.plot( t, Csoln, label="Enzyme/Substrate (ES)" )
matplotlib.pyplot.plot( t, Ssoln, label="Substrate (S)" )

matplotlib.pyplot.axhline( y=Cmax, color='r', linestyle='-', label="max(ES)" )

matplotlib.pyplot.title( "Enzyme Kinetics" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.xlabel( "t (sec)" )
matplotlib.pyplot.ylabel( "Concentration" )

matplotlib.pyplot.show()


# EOF
 