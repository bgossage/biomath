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
 S(t) .=. Susceptibles
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
    5) The poplulation is a constant N s.t.  S + I + R = N
    6) dti .=. incubation period

The ODE system is:

   dS/dt = -r * S(t) * I(t)
   dI/dt = r * S(t) * I(t) - (a+d) * I(t)
   dR/dt = (a-d) * I(t)
   dD/dt = d * I(t)

Challenge: Show that an epidemic can only occur if S(0) > a/r.
Epidemic <= for some t, I(t) > I0

 Below is the solution in python using numerical integration.

"""

import numpy
import scipy.integrate
from datetime import datetime

import matplotlib

#matplotlib.use('TkAgg')  ## set the back end (a wart)
import matplotlib.pyplot

# Parameters
r = 0.2  # The infection rate
a = 0.2  # removal rate


distancing_factor = 0.0 # 0.5

recovery_rate = 0.065
death_rate = 0.01
quarantine_rate = 0.0  # 0.1

r *= 1.0 - distancing_factor;
a = quarantine_rate + recovery_rate;
d = death_rate

S0 = a / r

print( "recovery rate = ", a )
print( "infection rate = ", recovery_rate )
print( "R0 = ", S0 )

conditions = "no intervention"
#
# Define a function that computes the derivatives at time t
# given the current populations of S, I, and R
#
def deriv( y, t, params ):

    S, I, R, D = y      # unpack current values

    r, a, d = params  # unpack parameters

    derivs = [ -r * S*I,
               r * S*I - a*I,
               (a-d) * I,
               d * I ]

    return derivs

# end function derivs ~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Initial values
S0 = 1.0
I0 = 0.00001
R0 = 0.0
D0 = 0.0

# Bundle parameters for ODE solver
params = [ r, a, d ]

# Bundle initial conditions for ODE solver
y0 = [ S0, I0, R0, D0 ]

# Make a time array for solution samples...
tStop = 200.0
tInc = 0.0010
t = numpy.arange(0., tStop, tInc)


# Solve using numerical integration...
psoln = scipy.integrate.odeint( deriv, y0, t, args=(params,) )

t += 14

Ssoln = psoln[:,0]
Isoln = psoln[:,1]
Rsoln = psoln[:,2]
Dsoln = psoln[:,3]


# Plot the solution...
#matplotlib.pyplot.plot( t, Ssoln, label="Susceptibles" )
matplotlib.pyplot.figure(0)
matplotlib.pyplot.plot( t, Isoln, label="Infectives, " + conditions, linestyle='dashed' )
matplotlib.pyplot.plot( t, Rsoln, label="Immune", linestyle=':', color='g' )
matplotlib.pyplot.plot( t, Dsoln, label="Deaths", color='r'  )

matplotlib.pyplot.title( "Epidemic" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.xlabel( "t (days)" )
matplotlib.pyplot.ylabel( "Population Fraction" )

matplotlib.pyplot.show()

# Read data...
#matplotlib.pyplot.figure(1)
total_pop = 8.74488E6

ny_data = numpy.genfromtxt('coronavirus-data/case-hosp-death.csv',
                           delimiter=',',dtype=None,encoding="utf8" )

datestrings = ny_data[1:,0]
deaths = ny_data[1:,3]

date0 = datetime.strptime(datestrings[0], '%m/%d/%y')

print( date0 )

size = datestrings.size
case_data = numpy.zeros(shape=(size,2))

row = 0
cum_deaths = 0
for datestr in datestrings:
    date = datetime.strptime( datestr, '%m/%d/%y')
    timedelta = date - date0
    case_data[row,0] = timedelta.days
    if deaths[row]:
        cum_deaths += float(deaths[row])
        case_data[row,1] = cum_deaths
    else:
         case_data[row,1] = 0.0
    row += 1

matplotlib.pyplot.title( "NY Data" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.xlabel( "t (days)" )
matplotlib.pyplot.ylabel( "Number" )

case_data[:,1] /= total_pop

matplotlib.pyplot.plot( case_data[:,0], case_data[:,1], label="Data" )
matplotlib.pyplot.show()

# EOF
