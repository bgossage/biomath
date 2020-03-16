# -*- coding: utf-8 -*-
"""
Created on Wed March 14, 2018

@author: bgossage
"""

"""
A simple epidemic model

Mathematical Biology, J.D. Murray
Chapter 19.1

Given a disease in which recovery confers immunity:
 S(t) = susceptibles
 I(t) = infectives
 R(t) = Removed (died, immune, or isolated)

(Often called "SIR" models.)

Assumptions:
    1) The rate increase of infectives (r) is proportional to the number of
       infectives and susceptibles
    2) The susceptibles are lost at the same rate r.
    3) The rate of removal of infectives (a) is proportional to the number of Infectives
       that is = a * I
    4) The incubation time is neglibile
    5) The poplulation is a constant N s.t.  S + I + R = N

The ODE system is:

   dS/dt = -r * S * I
   dI/dt = r * S * I - a * I
   dR/dt = a * I

Challenge: Show that an epidemic can only occur if S(0) > a/r.
Epidemic <= for some t, I(t) > I0

 Below is the solution in python using numerical integration.

"""

import numpy
import scipy.integrate

import matplotlib

#matplotlib.use('TkAgg')  ## set the back end (a wart)
import matplotlib.pyplot

# Parameters
    r = 0.25  # The infection rate (social distancing)
a = 0.1  # removal rate  (individual quarantine)


#
# Define a function that computes the derivatives at time t
# given the current populations of S, I, and R
#
def deriv( y, t, params ):

    S, I, R = y      # unpack current values

    r, a = params  # unpack parameters

    derivs = [ -r * S*I,
               r * S*I - a*I,
               a * I ]

    return derivs

# end function derivs ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Initial values
S0 = 1.0
I0 = 0.0001
R0 = 0.0

# Bundle parameters for ODE solver
params = [ r, a ]

# Bundle initial conditions for ODE solver
y0 = [ S0, I0, R0 ]

# Make a time array for solution samples...
tStop = 100.0
tInc = 0.001
t = numpy.arange(0., tStop, tInc)

print(t)

# Solve using numerical integration...
psoln = scipy.integrate.odeint( deriv, y0, t, args=(params,) )

Ssoln = psoln[:,0]
Isoln = psoln[:,1]
Rsoln = psoln[:,2]


# Plot the solution...
#matplotlib.pyplot.plot( t, Ssoln, label="Susceptibles" )
matplotlib.pyplot.plot( t, Isoln, label="Infectives" )
#matplotlib.pyplot.plot( t, Rsoln, label="Removed" )


matplotlib.pyplot.title( "Epidemic" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.xlabel( "t (days)" )
matplotlib.pyplot.ylabel( "Population" )

matplotlib.pyplot.show()


# EOF
