# This example simulates the CRN that computes the function
#
#       [Z]_\infty = f([A]_0, [B]_0, [C]_0)
#
#       where f(a, b, c) = c + a - min(a, b)
#

# import everything from the crn module
from crn import *

# define the species that will be present in your reactions
a, a1, a2, b, c, t, z = species("A A1 A2 B C T Z")

# define a CRN, aka the system of reactions
# reactions by default have a reaction constant of 1
sys = CRN(
    a >> a1 + a2,
    a1 + b >> t,
    c >> z,
    # to change the reaction constant, use this syntax
    (a2 >> z).k(2.5),
    z + t >> 0)

# simulate for the initial concentrations [A]_0 = 1.5 [B]_0 = 2.0.
# concentrations that are omitted are assumed to be zero
# t is an optional parameter for how far to carry out the simulation.
# if it is omitted, it's defaulted to 20.
sim = sys.simulate({a: 2.5, b: 2.0, c: 1.5}, t=5)

# you can access the specific time series of each species by accessing
# `sim` like you would a dictionary
z_time_series = sim[z]

# save a plot of the simulation with it's title to the file "sim.png".
# if the filename is omitted, the plot is shown on a new window.
sim.plot("sim.png", title="Example Simulation")


# The stochastic equivalent of the simulation & plotting above
sim = sys.stoch_simulate({a: 25, b: 20, c: 15}, t=5)
sim.plot("stoch_sim.png", title="Example Stochastic Simulation")

