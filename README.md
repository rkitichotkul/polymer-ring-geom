# polymer-ring-geom
Geometry of twistable ring polymer and nucleosome-bound DNA
This repository contains MATLAB scripts which generate initial configurations of twisted rings and nucleosome-bound DNA used in polymer simulations.

## twistedRing
This code is used to generate simple twisted rings. The twist is specified by the normal binormal vectors. Edit/run mainTwistedRing.m to use this script.
1. **Specify parameters**: Parameters of the ring can be specified at the top of mainTwistedRing.m
    - **L0**: the distance between beads (arbitrary unit)
    - **points**: the number of beads
    - **Tw**: the total amount of twist (radian)
2. **Run the script**
    - Run the whole script in order to ...
        - Obtain positions and orientation vectors of all beads
        - Plot all beads and the reference circle
        - Export the initial configuration
            - position **r0**
            - tangent vectors **u0**
            - binormal vectors **v0**
            - If the normal vectors **f0** are needed, simply add the following line to the end of the script:
            ```
            dlmwrite('f0', n_out, 'delimiter', '\t', 'precision', '%.12f')
            ```
    - The plot can be skipped. Simply run every section in mainTwistedRing.m except the **%% Plot** section.

## cyclicNucChain
This code is used to generate configurations of beads representing a ring of nucleosome-bound DNA. Each propagator traces the path forms by a nucleosomal helix and a linker. Each bead is placed at the entrance of the nucleosomal helix. In this script, the linear superhelix of nucleosome-bound DNA is looped into a ring by the parallel transport of the end-to-end axis of the superhelix. The orientations of beads at two ends may not match. More detailed explanation can be found in NucChainGeometry.pdf.

1. **Specify parameters**: Parameters of the ring can be specified at the top of cyclicNucChain.m
    - **L**: the contour length of the nucleosomal helix
    - **T**: the number of times DNA is wrapped around a nucleosome (radian)
    - **h**: the pitch of the nucleosomal helix
    - **lL**: the length of the linker
    - **tauDNA**: the intrinsic twist of DNA
    - **steps**: the number of beads
2. **Run the script**
    - Run the whole script in order to ...
        - Obtain positions and orientation vectors of all beads
        - Plot all beads and the reference circle
        - Export the initial configuration
            - position **r0**
            - tangent vectors **u0**
            - normal vectors **v0**
     - The plot can be skipped. Simply run every section in mainTwistedRing.m except the **%% Plot** section.
