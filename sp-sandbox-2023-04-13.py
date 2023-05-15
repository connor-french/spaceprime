# %%
#---
#title: "spaceprime sandbox"
#format: html
#---

# %% [markdown]
# # A place to toy around with spaceprime stuff
# Notes are in my Obsidian vault. 
# Note to self- I'm using the spaceprime conda environment.

# %% [markdown]
# ## raster import with rasterio

# %%
import rasterio

filepath = "scratch-data/iheringii.tif"

with rasterio.open(filepath) as src:
    print(src.profile)
    array = src.read(1, masked=True)

# %% [markdown]
# Plot the raster with matplotlib

# %%
#from matplotlib import pyplot

#pyplot.imshow(array, cmap='pink')

#pyplot.show() 

# %% [markdown]
# Manipulating raster data

# %%
import numpy as np
np.min(array)

# %% [markdown]
# Now to write a basic function! This will read in a single raster and create an msprime demography object with that single raster. I'm importing Jerome Kelleher's version of the 2D stepping stone model to work with.

# %%
import msprime
import numpy as np

def stepping_stone2d(initial_size, rate):
    assert len(initial_size.shape) == 2
    n, m = initial_size.shape


    N = n * m
    model = msprime.Demography.isolated_model(initial_size.reshape(N))
    M = model.migration_matrix
    for j in range(n):
        for k in range(m):
            index = j * m + k
            model.populations[index].name = f"pop_{j}_{k}"
            M[index, index - 1] = rate
            M[index, (index + 1) % N] = rate
            M[index, index - m] = rate
            M[index, (index + m) % N] = rate

            M[index - 1, index] = rate
            M[(index + 1) % N, index] = rate
            M[index - m, index] = rate
            M[(index + m) % N, index] = rate
    return model


# %%
def raster_to_demography(filepath):
    # open file and read in data, including mask
    with rasterio.open(filepath) as src:
        print(src.profile)
        a = src.read(1, masked=True)

    # if the nodata value is less than zero, set it to zero
    a[a<0] = 0

    a = np.round(a * 100)

    m = stepping_stone2d(a, 0.01)

    return(m)

    
    

# %%
mod = raster_to_demography("scratch-data/iheringii.tif")

# %%
print(mod[10000])

# %%
a = msprime.sim_ancestry(200)

# %%
print(a)

# %%

ts = msprime.sim_ancestry(samples={"pop_14_480": 2}, demography=mod)



# %%
#from IPython.display import SVG
#SVG(ts.draw_svg(y_axis=True))

print(ts)

# %%



