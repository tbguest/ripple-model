"""
A simple model for the formation dynamics of ripples by wind-blown sand, based
on the model by Nishimori and Ouchi (1992). 

The model incorporates the elementary sub-aerial sand transport processes of
grain saltation (i.e. near-bed 'downwind' transport) and surface creep 
(i.e., diffusion) in a nonlinear reaction-diffusion-type process on a lattice. 

Tristan Guest
8 Jan 2016
"""

import os
from constants import NSTEPS, WRITE_FRAMES
from initialize import initialize
from processes import saltate, diffuse
from output import to_image, show_image, write_gif
from pathlib import Path

def main():
    # Initialize the bed height (h) on the lattice 
    h = initialize()

    # Time steps
    for t in range(NSTEPS):

        # Saltation step (grains "jumping" downwind)
        h = saltate(h)

        # Diffusion step (representing "surface creep")
        h = diffuse(h)

        if WRITE_FRAMES:
            img = to_image(h)
            Path("./imgs").mkdir(parents=True, exist_ok=True)
            img.save(os.path.join("imgs", f"{t}.png"))

    # Show output
    show_image(h)


if __name__ == "__main__":
    main()
    
    if WRITE_FRAMES:
        write_gif()
