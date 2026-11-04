The settings for this directory are those of the fiducial model in the paper (Mstar = 10^6 Msun, 'fixed' opacity, alpha_MLT = 1, no SR).

The physical settings are mainly set in inlist_common, with other inlists changing them when necessary for stability. ./run_all will run all the stages from start to finish.

The stages are as follows:
 - inlist_initial: Initializes the quasi-star with a static black hole at its center, along with the inner boundary conditions for the quasi-star model. Output: initial.mod
 - inlist_expand: Keeps the black hole static (non-accreting), waiting for the envelope to expand and stabilize in response to the change in inner BCs. Output: expanded.mod
 - inlist_evolve: After the envelope has settled, the black hole growth begins and the quasi-star evolves until the end of its lifetime. Output: evolved.mod
