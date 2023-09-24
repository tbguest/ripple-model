import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from PIL import Image
import imageio
import os

def to_image(arr, maxVal=2):
    """Array to uint8, where maxVal approximates max(|arr|)"""
    img = (arr / (2 * maxVal)) + 0.5
    img[img > 1] = 1
    img[img < 0] = 0

    img = np.uint8(img * 255)
    img = Image.fromarray(img)

    return img


def show_image(h):
    # plot output
    fig, ax = plt.subplots()
    cax = ax.imshow(h, cmap=cm.coolwarm)

    # Add colorbar
    fig.colorbar(cax)

    plt.show()


def write_gif():
    filenames = sorted(os.listdir("imgs"))
    images = []
    for filename in filenames:
        images.append(imageio.imread(os.path.join("imgs", filename)))
    imageio.mimsave(os.path.join("imgs", "ripples.gif"), images)

