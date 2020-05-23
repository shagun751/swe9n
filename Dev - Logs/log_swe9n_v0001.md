## Initial development log

1. [Feature : Removed 'global' variables in basicVars [2020-05-23]](#log_swe9n_v0001_1)
1. [Feature : Include meshFreeMod [2020-05-23]](#log_swe9n_v0001_2)
1. [Feature : Include airyWaveModule [2020-05-23]](#log_swe9n_v0001_3)

### Attempting
- Clean the code
- Include meshfree module for 1st derivative


### List of Work
- [x] Removed "global" variables from basicVars
- [x] Removed rpimModule
- [x] Added meshFreeMod from Bsnq
	- [ ] Wet-dry in sweMFree
- [x] Airy wave module

-----------------------------------------------


<a name = 'log_swe9n_v0001_3' />

### Feature : Include airyWaveModule [2020-05-24]
- Module _airyWaveModule_ inside _modsInletBC.f90_ is added from Bsnq for tidal input.
- Will mostly be used only to calculate &eta;, as p = sqrt(gH) &eta;

-----------------------------------------------

<a name = 'log_swe9n_v0001_2' />

### Feature : Include meshFreeMod [2020-05-23]
- The wave-absorption require calculation of d(&eta;)/dn.
- This is done using the _meshFreeMod_ (_modsMFree.f90_) and _findNeiFromLinkList.f90_ developed in _Bsnq_.
- Slight modifications were made as the storing system in _swe_ and _Bsnq_ are different.
- A new subroutine _sweMFree_ was written based on the subroutine _setMFree_ inside _Bsnq_ to set up the mesh free derivatives.
- **Note: For wet-dry _sweMFree_ will have to be re-run.**
	- **Note: I have not yet incorporated wet-dry in _sweMFree_**
- The calculation was **tested and verified** using simple expression `eta = 2*coorx + 3*coory`, where the x and y derivatives will be 2 and 3 respectively.
	- Could not test using Paraview because the output in Paraview is in lon-lat.
- The coef for radius is set to 0.65 => 0.65 of max distance within the directly attached neighs of the nodes.
	- For the vertex nodes it covers the near 9 nodes
	- For the edge centres it should be more than sqrt(5)/sqrt(2) of the cell side => 0.65 x max distance
	- However for the 9th nodes which is the cell centre, all the neighs are the nodes within that one lement only. Therefore radius = 0.65 x max distance will lead to a no neighs except for itself.
		- This is solved by making radius = 2 x radius if numNeigh &lt; 3.

- Also **removed _rpimModule_** which was developed in the old Bsnq code for the same purpose of interpolation and derivatives calculation and the new _meshFreeMod_ is done much better using MLS.

-----------------------------------------------

<a name = 'log_swe9n_v0001_1' />

### Feature : Removed 'global' variables in basicVars [2020-05-23]
- The idea of having 'global' variables that are commonly used in subroutines, such as i, j, k, tmpr1, etc. was an idea from the older version of Bsnq development.
- However in the modular version of Bsnq I realised that this is quite a bad idea, as you might modify the variable inside a nfunction call thus messing up everything.
- Its also a bad idea when you want to run two objects of the same module (mainly a issue in Bsnq and not here).
- Hence I removed all non paramter type variables from _basicVars_.

-----------------------------------------------

## References
