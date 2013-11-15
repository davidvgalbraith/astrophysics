import math

#Astronomy 160 Saha equation Python Homework, by David Galbraith

"""
The following code, when run, yields: 

Electron density is 4.79151711261e+11 electrons/cm^3
Ratio of ionized to neutral hydrogen is 0.00477501985223
Ratio of ionized to neutral iron is 747.863225685
Pressure in container is 505.8179 barye

"""

# The Saha equation solved for electron density, returns the
# ionization fraction nH+/nHtotal of ionized to ground-state atoms
# given electron density ne, temperature t in kelvins, partition functions
# u1 and u2, and ionization energy ei of ion i from the fundamental level in eV
# (only works for hydrogen or if you can neglect ionizations after the first)

def saha(u1, u2, ei, t, ne):
    k = .00008617385 #Boltzmann constant in eV/K
    constant = 1.8*10**10 #2 * pi * electron mass * k / h^2
    ratio = (1/ne) * ((constant * t)**1.5) * ((2 * u2) / u1) * (math.e ** (-1 * ei/(k*t))) #equation 1.38
    return ratio / (1 + ratio) #equation 1.40

# The other form of the Saha equation, returns the ratio of
# atoms in ionization state i+1 to state i given partition
# functions u1 = U(i), u2 = U(i+1), ionization energy ei of
# the higher level in eV, temperature t in kelvins, and 
# electron density ne

def saharatio(u1, u2, ei, t, ne):
    k = .00008617385 #Boltzmann constant in eV/K
    constant = 1.8*10**10 #2 * pi * electron mass * k / h^2
    ratio = (1/ne) * ((constant * t)**1.5) * ((2 * u2) / u1) * (math.e ** (-1 * ei/(k*t)))
    return ratio

# Solves the given problem for arbitrary box volumes, numbers of atoms,
# and temperatures using the saha equations above

def astrophysics(vol, atoms, t):
    habundance = .9097/(.9097+.00003154) #abundance of hydrogen
    feabundance = 1-habundance #abundance of iron
    na = atoms/vol #density of atoms
    ne = na/500 #initial guess for electron density
    ne2 = 0 #Other initial guess for electron density
    iteration = 0
    while abs((ne-ne2)/ne) > .001 and iteration < 100000: #Error tolerance 1 part in 1000
        hionizationratio = saha(2.0, 1.0, 13.6, t, ne)
        feionizationratio = saha(10 ** 1.9, 10 ** 1.9, 7.9, t, ne)
        ne2 = ne #save old ne for comparison
        ne = hionizationratio * na * habundance + feionizationratio * na * feabundance #update ne
        iteration += 1 #converges after 617 iterations
    if iteration == 100000:
        print "Sorry but it didn't converge"
        return
    print "Electron density is %s electrons/cm^3" %ne
    hratio = saharatio(2.0, 1.0, 13.6, t, ne)
    print "Ratio of ionized to neutral hydrogen is %s" %hratio
    fe1ratio = saharatio(10 ** 1.9, 10 ** 1.9, 7.9, t, ne) #~750
    fe2ratio = saharatio(10 ** 1.9, 10 ** 1.9, 16.2, t, ne) * fe1ratio #~.04
    totalferatio = fe1ratio + fe2ratio #fe1/fe0 + fe2/fe0 = (fe1 + fe2)/fe0
    print "Ratio of ionized to neutral iron is %s" %totalferatio #higher energy levels abundance << .04
    pressure = na * 8.617 * (10**-16) * t #ideal gas law
    print "Pressure in container is %s barye" %pressure
    

astrophysics(1000000.0, 10.0**20 5870.0)
