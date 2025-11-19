NEBULINK - Monte Carlo Neutron Transport Simulator
==================================================

WHAT THIS PROGRAM DOES
----------------------
This program simulates how neutrons travel through materials to determine
if a nuclear configuration is "critical" - meaning it can sustain a chain
reaction.

Think of it like tracking thousands of tiny particles bouncing around
inside different materials. Each neutron can:
  - Get absorbed (disappear)
  - Scatter (bounce off an atom and change direction)
  - Cause fission (split a uranium atom, releasing more neutrons)
  - Leak out (escape the system)

THE KEY OUTPUT: k-effective
---------------------------
The main result is a number called "k-effective" (k-eff):

  k-eff < 1.0  = SUBCRITICAL  - Reaction dies out (safe)
  k-eff = 1.0  = CRITICAL     - Reaction sustains itself (reactor operation)
  k-eff > 1.0  = SUPERCRITICAL - Reaction grows (dangerous without control)

HOW IT WORKS
------------
The program uses the "Monte Carlo method" - it simulates many random neutron
paths and averages the results. This is the same approach developed at
Los Alamos in the 1940s for the Manhattan Project.

For each neutron:
1. Start at a random position with random direction
2. Travel until something happens (collision or boundary)
3. Randomly decide what happens based on physics probabilities
4. Repeat until neutron is absorbed or escapes
5. Count how many new neutrons are created vs. lost

THE TEST PROBLEM
----------------
The default simulation models a sphere of Uranium-235 (8.7 cm radius)
surrounded by a water reflector (10 cm thick). Water "moderates" (slows
down) neutrons, making them more likely to cause fission when they bounce
back into the uranium.

ENERGY GROUPS
-------------
Neutrons are tracked at three energy levels:
  - FAST:       Just born from fission, very high speed
  - EPITHERMAL: Slowing down, medium speed
  - THERMAL:    Fully slowed, most likely to cause fission in U-235

MATERIALS INCLUDED
------------------
  - Uranium-235: Fissile fuel
  - Water:       Neutron moderator (slows neutrons)
  - Graphite:    Alternative moderator
  - Lead:        Reflector/shielding
  - Concrete:    Shielding

TO RUN
------
Execute mcnp_mingw.exe - no input required. The program will run 50
generations of 1000 neutrons each and report k-effective with statistics.

EDUCATIONAL PURPOSE
-------------------
This is a simplified educational tool demonstrating nuclear physics concepts.
Real reactor analysis uses much more sophisticated codes with detailed
nuclear data libraries.
