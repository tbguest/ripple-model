import numpy as np
from constants import NX, NY, Q, L0, B, D


# saltation step
def saltate(h):
    # Precompute the saltation distances at each lattice point for this timestep
    L = L0 + B * h
    L[L < 0] = 0

    # Integerize
    L_integer = np.floor(L)

    for j in range(0, NX):
        for i in range(0, NY):
            h[i, j] = h[i, j] - Q
            jump = int(L_integer[i, j])
            if j + jump <= (NX - 1):
                h[i, j + jump] = h[i, j + jump] + Q
            else:
                wrap = j + jump - NX
                h[i, wrap] = h[i, wrap] + Q
    return h


# Diffusion step
def diffuse(h):
    # initialize (NN is "nearest neighbour")
    NNsum = np.zeros((NY, NX))

    # Creep (i.e. dispersion) step
    # 1) boundaries
    for i in range(1, NY - 1):
        NNsum[i, 0] = (
            h[i, 1]
            + h[i, -1]
            + h[i + 1, 1] / 2
            + h[i + 1, 0]
            + h[i + 1, -1] / 2
            + h[i - 1, 1] / 2
            + h[i - 1, 0]
            + h[i - 1, -1] / 2
        )
        NNsum[i, -1] = (
            h[i, 0]
            + h[i, -2]
            + h[i + 1, 0] / 2
            + h[i + 1, -1]
            + h[i + 1, -2] / 2
            + h[i - 1, 0] / 2
            + h[i - 1, -1]
            + h[i - 1, -2] / 2
        )

    for j in range(1, NX - 1):
        NNsum[0, j] = (
            h[1, j]
            + h[-1, j]
            + h[1, j + 1] / 2
            + h[0, j + 1]
            + h[-1, j + 1] / 2
            + h[1, j - 1] / 2
            + h[0, j - 1]
            + h[-1, j - 1] / 2
        )
        NNsum[-1, j] = (
            h[0, j]
            + h[-2, j]
            + h[0, j + 1] / 2
            + h[-1, j + 1]
            + h[-2, j + 1] / 2
            + h[0, j - 1] / 2
            + h[-1, j - 1]
            + h[-2, j - 1] / 2
        )

    # 2) Corners
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

    # 3) Remainder of lattice
    for j in range(1, NX - 1):
        for i in range(1, NY - 1):
            NNsum[i, j] = (
                h[i + 1, j - 1] / 2
                + h[i, j - 1]
                + h[i - 1, j - 1] / 2
                + h[i + 1, j]
                + h[i - 1, j]
                + h[i + 1, j + 1] / 2
                + h[i, j + 1]
                + h[i - 1, j + 1] / 2
            )

    h = h + D * (NNsum / 6 - h)
    return h
