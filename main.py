import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval, ImageNormalize
from regions import CirclePixelRegion, PixCoord
import numpy as np
import os

fits_dir = r"C:\Users\chend\Desktop\Projects\CepheidVariables\data"
fits_list = os.listdir(fits_dir)
reg_dir = r"C:\Users\chend\Desktop\Projects\CepheidVariables\reg"
reg_list = os.listdir(reg_dir)
reg_list = [file for file in reg_list if file.endswith('.txt')]
fits_list = [file for file in fits_list if file.endswith('.fits')]

#hdulist = [fits.open(file) for file in fits_list]
hdulist = fits.open(fits_list[0])
data = hdulist[0].data
wcs = WCS(hdulist[0].header)

hdulist.close()
h = fits.getheader(fits_list[0])
# Read the text file
# file_path = '2x1.txt'  # Replace with the actual file path
with open(os.path.join(reg_dir, reg_list[0]), 'r') as file:
    lines = file.readlines()

# Initialize empty arrays for x, y, and radius
x_array = []
y_array = []
radius_array = []

# Parse each line and extract coordinates and radius
for line in lines:
    if line.startswith('circle'):
        parts = line.split('(')[1].split(')')[0].split(',')
        x_array.append(float(parts[0]))
        y_array.append(float(parts[1]))
        radius_array.append(float(parts[2]))

# Convert arrays to NumPy arrays
x_array = np.array(x_array)
y_array = np.array(y_array)
radius_array = np.array(radius_array)

flux_list = []

# lombsceargle

for x, y, r in zip(x_array, y_array, radius_array):
    center_pix = PixCoord(x-1, y-1)
    radius_pix = 0.75*r  # Specify the radius in pixels
    circle_region = CirclePixelRegion(center=center_pix, radius=radius_pix)

    # Extract data within the circular region
    region_data = circle_region.to_mask(mode='subpixels').cutout(data)
    region_uncertainty = circle_region.to_mask(mode='subpixels').cutout(data)
    zscale = ZScaleInterval()
    norm = ImageNormalize(region_data, interval=zscale)

    plt.imshow(region_data, origin='lower', cmap='viridis', norm=norm)
    plt.title((x, y))
    plt.show()
    print(np.sum(region_data))
    flux_list.append(region_data)

flux_pairs = []
for i in range(0, len(flux_list), 2):
    flux_pairs.append((flux_list[i], flux_list[i+1]))


def background(c1, c2):
    if np.sum(c2) > np.sum(c1):
        bg = (np.sum(c2)-np.sum(c1))/(np.size(c2)-np.size(c1))
        return bg
    elif np.sum(c2) < np.sum(c1):
        bg = np.mean(c2)
        return bg
    else:
        bg = 0


background_set = [background(i[0], i[1]) for i in flux_pairs]

star_flux_set = []

# for i in range(len(background_set)):
#  star_flux.append
