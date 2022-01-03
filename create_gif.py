import imageio
import os


filenames = sorted(os.listdir("imgs"))

images = []
for filename in filenames:
    images.append(imageio.imread(os.path.join("imgs", filename)))
imageio.mimsave(os.path.join("imgs", "ripples.gif"), images)

