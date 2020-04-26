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
ny_data = numpy.genfromtxt('coronavirus-data/case-hosp-death.csv',
                           delimiter=',',dtype=None,encoding="utf8" )
import matplotlib

#matplotlib.use('TkAgg')  ## set the back end (a wart)
import matplotlib.pyplot

# Parameters
r = 0.7 # The infection rate
a = 0.01  # removal rate

distancing_factor = 0.6 # 0.5

recovery_rate = 0.1
death_rate = 0.04
quarantine_rate = 0.0  # 0.1

#r *= 1.0 - distancing_factor
a = quarantine_rate + recovery_rate;
d = death_rate


ld_start = 15.0 # start of lockdown
ld_end = 45.0   # end of lockdown

#S0 = a / r

print( "recovery rate = ", a )
print( "infection rate = ", r )
#print( "R0 = ", S0 )

conditions = "no intervention"
#
# Define a function that computes the derivatives at time t
# given the current populations of S, I, and R
#
def deriv( y, t, params ):

    S, I, R, D = y      # unpack current values

    r, a, d = params  # unpack parameters

    derivs = [ -r * S*I,            # Suseptibles
               r * S*I - (a+d)*I,   # Infectives
               (a) * I,           # Survived/Immune
               d * I ]              # Deaths

    return derivs

# end function derivs ~~~~~~~~~~~~~~~~~~~~~~~~~~~



# Initial values
I0 = 0.002
R0 = 0.2
S0 = 1.0 - R0 - I0
D0 = 0.0

# Bundle parameters for ODE solver
params = [ r, a, d ]

# Bundle initial conditions for ODE solver
y0 = [ S0, I0, R0, D0 ]

# Make a time array for solution samples...
tStart = 0.0
tStop = ld_start
tInc = 0.010
t_0 = numpy.arange(tStart, tStop, tInc)

# Solve using numerical integration...
psoln_0 = scipy.integrate.odeint( deriv, y0, t_0, args=(params,) )

tStart = ld_start
tStop = ld_end
r_start = r;
r *= (1.0 - distancing_factor)
params = [ r, a, d ]
t_1 = numpy.arange(tStart, tStop, tInc)
# Bundle initial conditions for ODE solver
last = len(psoln_0[:,0])-1
S0 = psoln_0[last,0]
I0 = psoln_0[last,1]
R0 = psoln_0[last,2]
D0 = psoln_0[last,3]
y0 = [ S0, I0, R0, D0 ] 

# Solve using numerical integration...
psoln_1 = scipy.integrate.odeint( deriv, y0, t_1, args=(params,) )

tStart = ld_end
tStop = 100.0
r = r_start
params = [ r, a, d ]
t_2 = numpy.arange(tStart, tStop, tInc)
# Bundle initial conditions for ODE solver
last = len(psoln_1[:,0])-1
S0 = psoln_1[last,0]
I0 = psoln_1[last,1]
R0 = psoln_1[last,2]
D0 = psoln_1[last,3]
y0 = [ S0, I0, R0, D0 ] 

# Solve using numerical integration...
psoln_2 = scipy.integrate.odeint( deriv, y0, t_2, args=(params,) )


Ssoln_0 = psoln_0[:,0]
Isoln_0 = psoln_0[:,1]
Rsoln_0 = psoln_0[:,2]
Dsoln_0 = psoln_0[:,3]

Ssoln_1 = psoln_1[:,0]
Isoln_1 = psoln_1[:,1]
Rsoln_1 = psoln_1[:,2]
Dsoln_1 = psoln_1[:,3]

Ssoln_2 = psoln_2[:,0]
Isoln_2 = psoln_2[:,1]
Rsoln_2 = psoln_2[:,2]
Dsoln_2 = psoln_2[:,3]


# Plot the solution...
#matplotlib.pyplot.plot( t, Ssoln, label="Susceptibles" )
matplotlib.pyplot.figure( 0, figsize=(8,7) )

#matplotlib.pyplot.plot( t, Isoln, label="Infectives, " + conditions, linestyle='dashed' )
matplotlib.pyplot.plot( t_0, Rsoln_0, label="Immune Before", linestyle=':', color='r' )
matplotlib.pyplot.plot( t_1, Rsoln_1, label="Immune Shut", linestyle='--', color='b' )
matplotlib.pyplot.plot( t_2, Rsoln_2, label="Immune Restart", linestyle='-.', color='g' )

matplotlib.pyplot.plot( t_0, Dsoln_0, label="Total Death Before", linestyle=':', color='r'  )
matplotlib.pyplot.plot( t_1, Dsoln_1, label="Total Deaths Shut", linestyle='--', color='b'  )
matplotlib.pyplot.plot( t_2, Dsoln_2, label="Total Deaths Restart", linestyle='-.', color='g'  )


matplotlib.pyplot.plot( t_0, Isoln_0, label="Infectives Before", linestyle=':', color='r'  )
matplotlib.pyplot.plot( t_1, Isoln_1, label="Infectives Shut", linestyle='--', color='b'  )
matplotlib.pyplot.plot( t_2, Isoln_2, label="Infectives Restart", linestyle='-.', color='g'  )

matplotlib.pyplot.title( "Epidemic" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.xlabel( "t (days)" )
matplotlib.pyplot.ylabel( "Population Fraction" )

matplotlib.pyplot.show()

# Read data...
#matplotlib.pyplot.figure(1)
total_pop = 1585873

total_pop = 1.0E5

ny_data = numpy.genfromtxt('coronavirus-data/case-hosp-death.csv',
                           delimiter=',',dtype=None,encoding="utf8" )

datestrings = ny_data[1:,0]
deaths = ny_data[1:,3]

date0 = datetime.strptime(datestrings[0], '%m/%d/%y')

print( date0 )

last = len(Dsoln_2)-1
print( "Final Death Toll: ", Dsoln_2[last])

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
matplotlib.pyplot.ylabel( "Fraction" )

case_data[:,1] /= total_pop


matplotlib.pyplot.plot( case_data[:,0], case_data[:,1], label="Deaths" )
matplotlib.pyplot.legend( loc='best' )
matplotlib.pyplot.show()

# EOF
