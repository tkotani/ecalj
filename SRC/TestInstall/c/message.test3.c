         The C test tests the code's implementation of homogeneous background mode
         It checks that program lmf generates the correct total energy and 
         ionization potential for a single C atom.

         Note that:

	*The total energy of the neutral atom computed by lmf (-74.996 Ry) is very close
         to the free-atom energy computed by lmfa (-74.995 Ry), and that the free-atom 
         potential is approximately self-consistent; 

	*Also the SUM of total energy of the ionized system (-74.443) and the
         interaction of a point charge with a homogeneous background E^I,
         estimated from E^I=9/5R, with 4*pi*R^3/3=vol => R=6.20 and E^I = 0.290 Ry
         evaluates to -74.443+0.290 = -74.153 Ry, is close to the total energy of
         the positively charged free ion as computed by lmfa (-74.171 Ry).


        We should have; ---------------------
         Energy of charged system         = .5524587
         Estat energy q*q/9/5/<r>         = 0.2902

           Corrected charged system energy  =  0.842659  <---
           Energy of neutral system         = -.0012168  <---
         ----------------------------------------------         
          difference                        = 0.843876
