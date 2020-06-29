# Cepstral Identification of Autoregressive Systems
Code for numerical illustrations accompanying the paper Cepstral Identification of Autoregressive Systems.
## Summary
The (Matlab) code in this repository recreates the numerical examples in the manuscript "Cepstral Identification of Autoregressive Systems", a preprint of which can be found on [arXiv](placeholder), by Oliver Lauwers, Christof Vermeersch and Bart De Moor, which develops a novel, extremely efficient system identification technique, to estimate the difference equation of a system, starting from the power cepstra of input and output signals. The algorithm then consists of an exact solution to a set of difference equations.
This manuscript has been sent in to be considered for publication to Automatica.

Several numerical examples are available:
 - A simple numerical illustration, consisting of a synthetic model that is fully known. We simulate input and output data, and identify the system, to compare to the original system.
 - An example with synthetic data, to show robustness and convergence behaviour. 
 - An example showing the improved bias behaviour of the cepstral system identification technique to identify representative dynamics of a cluster of signals coming from the same system.
 - A historical, canonical real-life data set of Yule's sunspot numbers.
 - A real-life application on structural health monitoring via autoregressive models.

Bear in mind that this code is not meant as a fully working software package, but serves merely as an illustration accompanying the manuscript mentioned earlier.

## Reference
When using this code or discussing results of cepstral system identification, please refer to [this paper](placeholder).

The data for the structural health monitoring application was provided by the Los Alamos National Laboratory.
It comes from the "Bookshelf Frame Structure -DSS 2000" data set, at https://www.lanl.gov/projects/national-security-education-center/engineering/software/shm-data-sets-and-software.php

