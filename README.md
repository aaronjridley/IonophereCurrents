# IonophereCurrents

This is a code that will take an amie-style output file, read it, and make a bunch of calculations to derive the currents from it.  
It uses the potential and the Hall and Pedersen conductances to create altitude profiles of these quantities, then derives the Hall and Pedersen 
currents that would be associated with those files. Then it creates the field-aligned component by calculating the divergence.  It outputs all of
this to a gitm-style output file.
