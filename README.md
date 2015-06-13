# nrnanalytic

Find numerical values for the analytical expressions presented in the papers:

- Monai H, Omori T, Okada M, Inoue M, Miyakawa H, Aonishi T. 
"An analytic solution of the cable equation predicts frequency 
preference of a passive shunt-end cylindrical cable in response 
to extracellular oscillating electric fields."
Biophys J. 2010 Feb 17;98(4):524-33.

- W. Ying and C. S. Henriquez, "Hybrid finite element method 
for describing the electrical response of biological cells
to applied fields," IEEE Transactions on Biomedical Engineering, 
vol. 54, no. 4, pp. 611-620, Apr. 2007. 

- Tadej Kotnik, Damijan Miklavcic, Tomaz Slivnik, Time course of 
transmembrane voltage induced by time-varying electric fields 
a method for theoretical analysis and its application, 
Bioelectrochemistry and Bioenergetics, Volume 45, Issue 1, 
March 1998.

- Bedard C, Kroger H and Destexhe A.  Modeling extracellular field
potentials and the frequency-filtering properties of extracellular
space. Biophysical Journal 86: 1829-1842, 2004.

The expressions correspond to the response of a cable to a DC step current 
input, the cross section of a cylindrical cell to a homogeneous
field, a sphere to a homogeneous field, and the potential 
for a point source current.

Typical parameters are specified as global variables for simplicity
in the Python module and can be changed if necessary before 
function invocation.

Copyright Andres Agudelo-Toro (https://sites.google.com/site/aagudelotoro/)
