"""
A simple model for the formation dynamics of ripples by wind-blown sand, based
on the model by Nishimori and Ouchi (1992). 

The model incorporates the elementary sub-aerial sand transport processes of
grain saltation (i.e. near-bed 'downwind' transport) and surface creep 
(included herein as a localized diffusive process). 

Tristan Guest
8 Jan 2016
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image


def to_image(arr, maxVal=2):
    """Array to uint8, where maxVal approximates max(|arr|)"""
    img = (arr / (2 * maxVal)) + 0.5
    img[img > 1] = 1
    img[img < 0] = 0

    img = img * 255

    return np.uint8(img)


def main():

    # define lattice size
    nx = 100
    ny = 100
    x = np.linspace(0, nx - 1, nx)
    y = np.linspace(0, ny - 1, ny)

    # amplitude of the random initial bed state
    eps = 0.1
    h0 = eps * np.random.randn(ny, nx)

    # q is the transferred height of saltated grains (i.e. grain diameter).
    # should be related to eps
    q0 = 1
    q = q0 * eps

    # Nishimori-Ouchi saltation control parameters
    # b is proportional to the mean wind velocity a grain experiences in flight
    # l0 is proportional to the wind force
    L0 = 5.0
    b = 2.0

    # creep (i.e. dispersion) coefficient
    D = 0.25

    # timesteps
    nstep = 100

    # initialize
    hAll = np.zeros((ny, nx))
    h = h0

    # time step loop
    for kk in range(1, nstep):

        # compute (non-negative) saltation distance at each lattice point
        L = L0 + b * h
        L[L < 0] = 0

        # integerize
        Lint = np.floor(L)

        # saltation step
        for jj in range(0, nx):

            for ii in range(0, ny):

                h[ii, jj] = h[ii, jj] - q

                jump = int(Lint[ii, jj])

                if jj + jump <= (nx - 1):
                    h[ii, jj + jump] = h[ii, jj + jump] + q

                else:
                    wrap = jj + jump - nx
                    h[ii, wrap] = h[ii, wrap] + q

        # initialize
        NNsum = np.zeros((ny, nx))

        # creep (i.e. dispersion) step
        # 1) boundaries
        for ii in range(1, ny - 1):
            NNsum[ii, 0] = (
                h[ii, 1]
                + h[ii, -1]
                + h[ii + 1, 1] / 2
                + h[ii + 1, 0]
                + h[ii + 1, -1] / 2
                + h[ii - 1, 1] / 2
                + h[ii - 1, 0]
                + h[ii - 1, -1] / 2
            )
            NNsum[ii, -1] = (
                h[ii, 0]
                + h[ii, -2]
                + h[ii + 1, 0] / 2
                + h[ii + 1, -1]
                + h[ii + 1, -2] / 2
                + h[ii - 1, 0] / 2
                + h[ii - 1, -1]
                + h[ii - 1, -2] / 2
            )

        for jj in range(1, nx - 1):
            NNsum[0, jj] = (
                h[1, jj]
                + h[-1, jj]
                + h[1, jj + 1] / 2
                + h[0, jj + 1]
                + h[-1, jj + 1] / 2
                + h[1, jj - 1] / 2
                + h[0, jj - 1]
                + h[-1, jj - 1] / 2
            )
            NNsum[-1, jj] = (
                h[0, jj]
                + h[-2, jj]
                + h[0, jj + 1] / 2
                + h[-1, jj + 1]
                + h[-2, jj + 1] / 2
                + h[0, jj - 1] / 2
                + h[-1, jj - 1]
                + h[-2, jj - 1] / 2
            )

        # 2) corners
        NNsum[0, 0] = (
            h[0, 1]
            + h[1, 0]
            + h[-1, 0]
            + h[0, -1]
            + h[1, 1] / 2
            + h[-1, -1] / 2
            + h[-1, 1] / 2
            + h[1, -1] / 2
        )
        NNsum[0, -1] = (
            h[0, 0]
            + h[-1, -1]
            + h[0, -2]
            + h[1, -1]
            + h[1, 0] / 2
            + h[1, -2] / 2
            + h[-1, -2] / 2
            + h[-1, 0] / 2
        )
        NNsum[-1, 0] = (
            h[-1, 1]
            + h[-1, -1]
            + h[-2, 0]
            + h[0, 0]
            + h[0, 1] / 2
            + h[0, -1] / 2
            + h[-2, 1] / 2
            + h[-2, -1] / 2
        )
        NNsum[-1, -1] = (
            h[-1, -2]
            + h[-1, 0]
            + h[-2, -1]
            + h[0, -1]
            + h[0, 0] / 2
            + h[0, -2] / 2
            + h[-2, 0] / 2
            + h[-2, -2] / 2
        )

        # 3) remainder of lattice
        for jj in range(1, nx - 1):
            for ii in range(1, ny - 1):
                NNsum[ii, jj] = (
                    h[ii + 1, jj - 1] / 2
                    + h[ii, jj - 1]
                    + h[ii - 1, jj - 1] / 2
                    + h[ii + 1, jj]
                    + h[ii - 1, jj]
                    + h[ii + 1, jj + 1] / 2
                    + h[ii + 1, jj]
                    + h[ii + 1, jj - 1] / 2
                )

        h = h + D * (NNsum / 6 - h)

        # to create a gif
        img = to_image(h)
        im = Image.fromarray(img)
        im.save(os.path.join("imgs", f"{kk:03d}.png"))

        hAll = h

    # plot output
    fig, ax = plt.subplots()
    cax = ax.imshow(hAll, cmap=cm.coolwarm)

    # Add colorbar
    cbar = fig.colorbar(cax)

    plt.show()


if __name__ == "__main__":
    main()
