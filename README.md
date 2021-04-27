# Harmonic Oscillator

This repo contains some codes to work with **quantum** and **classical harmonic oscillator**.
Folders contain the main results and the programs are in the scope of the repo.

## Code descriptions

- `analytical.f`: computes the probability of finding a particle for the classical harmonic oscillator (analitically), together with the
density probability of the quantum harmonic oscillator in order to compare each order and test the Correspondence Principle;
- `plot-analytical.py`: plots the comparison: classical probability X quantum density of probability, using analytical results. Besides, this program moves the data files created in the main directory by `analytical.f` to the folder `analytical/` 
- `energies.f`: computes eigenvalues for the energy of the quantum harmonic oscillator using Numerov's method;
- `hermite.f`: computes the Hermite polynomials until order 5;
- `plot_hermite.py`: plots the Hermite polynomials until order 5. Besides, this program moves the data files created in the main directory by `hermite.f` to the folder `hermite/`;
- `numerov.f`: computes the probability density for the quantum harmonic oscillator using Numerov's method;
in order to prove the Correspondence Principle;
- `plot-numerov.py`: plots the comparison Numerov's X analytical probability densities. Besides, this program moves the data files created in the main directory by `numerov.f` to the folder `numerov/`;

An interesting explanation about the theory and plots of this repo is available in my blog: [natalidesanti](https://natalidesanti.github.io/blog-post-3/).