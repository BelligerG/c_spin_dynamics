# c_spin_dynamics

The class "SpinDynamics" will allow the generation of objects that can be used to calculate the yields to measure the MFEs for different orientations and numbers of radicals. By providing the object with the spins, the class can handle spins other than 1/2, which can provide hyperfine interactions of non-standard spin nuclei.

## Getting Started

For an indepth guide, see the tutorial below.  For a quick start guide see the annotated code below:

###Tutorial

This guide is intended for use by those wanting to run spin dynamics calculations, but don't have much C++ experience.

The first step is to instantiate the class, this can be done by a single lines e.g.:

SpinDynamics SPD;

Following this, we want to set the number of electrons in our system.

SPD.SetNumberOfElectrons(3);

So in our example system there are 3 electrons.  After this we need to set our spin operators. This is possible to set from outside the object, but the easiest way is just to provide the spins as an array:

float spins [4] = {0.5, 0.5, 0.5, 0.5};

Here we specify that we have 4 spin particles (3e- and 1n+) all at spin 1/2.  The code will use the first n spins (where n is the number of electrons - in our case 3) as the electrons.  Generating the spin operators is easy with:

SPD.spinsToSpinOperators(spins, 4);

We pass the spins and number of spins into the function and we set internally the spin operators.

