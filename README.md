# Figure package

## Installation

Requires:

  * [Birch.Cairo](https://github.com/lawmurray/Birch.Cairo)

To build and install, use:

    birch build
    birch install
    
## Usage
    
The run, use:

    birch draw [options]

This outputs frames as PDF images into the `figs/` directory, showing
iterations of conditional Sequential Monte Carlo (SMC). Options include

  * `--width` (default 1024): width, in pixels, of the images.
  * `--height` (default 512): height, in pixels, of the images.
  * `-M` (default 4): Number of iterations.
  * `-T` (default 20): Number of time steps.
  * `-N` (default 12): Number of particles.
  * `--seed`: Pseudorandom number seed.
