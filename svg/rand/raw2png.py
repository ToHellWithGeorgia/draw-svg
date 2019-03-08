import os
import sys
import numpy as np

fn = "rand_000.raw"

with open(fn) as fh:
    raw = fh.readlines()
    # Read the sizes of the image
    sizes = raw[0].split(" ")[:2]
    sizes = [int(size) for size in sizes]
    sizes.append(3)

    # Create a numpy array to store the image
    raw = raw[1:]
    raw = [pix.split(' ')[:3] for pix in raw]
    arr = [[int(rgb) for rgb in pix] for pix in raw]
    nparr = np.array(arr).reshape(sizes)
    # print(nparr.shape)

