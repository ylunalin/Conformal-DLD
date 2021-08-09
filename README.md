# Conformal-DLD

This is a discrete computation model to simulation diffusion-limited dissolution.
Coming soon.

Authors: Yuexia Luna Lin, Chris H. Rycroft

_John A. Paulson School of Engineering and Applied Sciences, Harvard University_
## Software/library requirements:

- GNUMake
- GNUGCC compiler
- GSL (GNU Scientific Library)
- Perl         (For image processing only)
- Povray 3.7   (For image processing only)

## Get started:
1. Download or <code> git clone </code> Conformal-DLD repository.
2. Modify the config.mk file to configure which compiler to use and basic compiler and linker flags. The default config.mk file provide Mac Ports installed GCC 9 on Macs, another sample config for Linux systems is provided in the configs/ folder.
3. After these changes, type <code> make </code> in the commandline in Conformal-DLD/ directory, which will build executables flat and radial. These executables will use OpenMP to parallelize the computation if it is available.
4. To run an executable, simulation parameters can be passed in on the command line. These command-line arguments can be printed by
<code> ./flat -h </code>
or 
<code> ./radial -h </code>

More complete documentation is under development.
In the meantime, for any questions, feel free to contact the authors of this repository.

## Acknowledgement
This work has been partially supported by the Applied Mathematics Program of the U.S. DOE Office of Science Advanced Scientific Computing Research under contract number DE-AC02-05CH11231, the Department of Energy Computational Science Graduate Fellowship.

