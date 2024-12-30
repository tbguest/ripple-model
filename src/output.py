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

    # Calculate new dimensions
    width, height = img.size
    scale_factor = 4
    new_size = (width * scale_factor, height * scale_factor)

    # Resize using high-quality interpolation
    img = img.resize(new_size, resample=Image.Resampling.NEAREST)

    return img


def show_image(h):
    # plot output
    fig, ax = plt.subplots()
    cax = ax.imshow(h, cmap=cm.coolwarm)

    # Add colorbar
    fig.colorbar(cax)

    plt.show()


def write_gif():
    # Natural sort the filenames (e.g., frame1.png, frame2.png, frame10.png)
    filenames = sorted(
        os.listdir("imgs"), key=lambda x: int("".join(filter(str.isdigit, x)))
    )
    images = []
    for filename in filenames:
        images.append(imageio.imread(os.path.join("imgs", filename)))
    imageio.mimsave(
        "ripples.gif",
        images,
        # duration=0.05,
        optimize=True,  # Enable optimization
        subrectangles=True,  # Only update changed pixels between frames
    )
