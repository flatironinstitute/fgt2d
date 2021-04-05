Fast Gauss transform in two and three dimensions - point case.

2d/sx: contains the SX scheme, that is, form Hermite, Hermite to SOE and SOE/X,
translate SOE and SOE/X (mp to mp, mp to loc, loc to loc), SOE and SOE/X to local,
local eval.

2d/planewave: contains the planewave scheme, that is, form Hermite, Hermite to plane wave,
translate plane wave (mp to mp, mp to loc, and loc to loc), plane wave to local, local
eval.

2d/nufft: contains plane wave only scheme, that is, use type 1 NUFFT to form plane wave
expansion at the cutoff level, translate PW exp (mp to loc only), use type 2 NUFFT to
evaluate PW expansion.
---------------------------------------------------------------------
3d/planewave: planewave scheme for 3D.

3d/nufft: nufft scheme for 3D.
---------------------------------------------------------------------

To run, type "make lib", then "make test".


For run nufft schemes, install finufft, and change
the line "FINUFFT = /home/shidong/finufft" to the directory where you
install the finufft dynamic library.

