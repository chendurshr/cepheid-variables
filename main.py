# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 21:55:23 2023

@author: chend
"""

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.visualization import ZScaleInterval, ImageNormalize
from regions import CirclePixelRegion, PixCoord
import numpy as np
import os
from scipy.signal import lombscargle
pixel_area = 0.099595926**2


def region_coord(reg_list):
    fits_region_coord = []
    for region in reg_list:

        with open(os.path.join(reg_dir, region)) as region_file:
            number = 0 if region.startswith('2') else 7

            reg_array = []
            x_array = []
            y_array = []
            radius_array = []

            lines = region_file.readlines()

            for line in lines:

                if line.startswith('circle'):
                    number += 0.5
                    star_number = np.ceil(number)
                    parts = line.split('(')[1].split(')')[0].split(',')
                    reg_array.append(f"Star {star_number:1.0f}")
                    x_array.append(float(parts[0]))
                    y_array.append(float(parts[1]))
                    radius_array.append(float(parts[2]))
            reg_array = np.array(reg_array)
            x_array = np.array(x_array)
            y_array = np.array(y_array)
            radius_array = np.array(radius_array)
        fits_region_coord.append(
            list(zip(reg_array, x_array, y_array, radius_array)))
    return fits_region_coord


def data_gen(fits_list):
    data = []
    wcs = []
    header = []

    for fit in fits_list:
        fit_path = os.path.join(fits_dir, fit)
        hdulist = fits.open(fit_path)
        data.append(hdulist[0].data)
        wcs.append(WCS(hdulist[0].header))
        header.append(fits.getheader(fit))

    data = np.array(data)
    return data, wcs, header


def region_data(reg_data, region_coord, image_date):
    count = 0
    previous_region_data = np.array([0])
    region_data = np.array([1])
    regional_flux_list = []
    for image in range(reg_data.shape[0]):
        image_regions = region_coord[image]

        current_data_set = reg_data[image]
        n = image_regions[:, 0]
        x = image_regions[:, 1]
        y = image_regions[:, 2]
        r = image_regions[:, 3]
        previous_region_data = None

        for n_val, x_val, y_val, r_val in zip(n, x, y, r):
            x_val = float(x_val)
            y_val = float(y_val)
            r_val = float(r_val)

            max_iterations = 100
            for iteration in range(max_iterations):
                radius_pix = 0.72 * r_val

                center_pix = PixCoord(x_val - 1, y_val - 1)

                circle_region = CirclePixelRegion(
                    center=center_pix, radius=radius_pix)

                region_data = circle_region.to_mask(
                    mode='subpixels').cutout(current_data_set)

                if not np.array_equal(region_data, previous_region_data):
                    break

                r_val *= 1.1

            regional_flux_list.append((str(n_val), region_data))
            previous_region_data = np.array(region_data)
    return regional_flux_list


def pair(flux_list):
    flux_pairs = []
    for i in range(0, len(flux_list), 2):
        flux_pairs.append(
            (flux_list[i][0], flux_list[i][1], flux_list[i+1][1]))
    return flux_pairs


def plot(region_data, image_date):
    count = -.5
    for index, (n_val, region) in enumerate(region_data):
        count += 0.5 if n_val == "Star 1" else 0
        count += 0.5 if n_val == "Star 8" else 0

        zscale = ZScaleInterval()
        norm = ImageNormalize(region, interval=zscale)

        plt.imshow(region, origin='lower', cmap='viridis', norm=norm)
        plt.title(
            n_val + f" Day: {image_date[int(count)]:1.0f}")

        plt.show()


def background(c1, c2):
    if np.sum(c2) > np.sum(c1):
        bg = (np.sum(c2)-np.sum(c1))/(np.size(c2)-np.size(c1))
        return bg
    elif np.sum(c2) < np.sum(c1):
        bg = np.mean(c2)
        return bg
    else:
        bg = 0


fits_dir = r"C:\Users\chend\Desktop\Projects\CepheidVariables"
fits_list = os.listdir(fits_dir)
reg_dir = r"C:\Users\chend\Desktop\Projects\CepheidVariables\reg"
reg_file_list = os.listdir(reg_dir)

chip2_reg_list = [file for file in reg_file_list if file.endswith(
    '.txt') and file.startswith("2")]
chip2_reg_list = sorted(chip2_reg_list, key=lambda x: (
    int(x.split('x')[1].split('.')[0]), x))

chip3_reg_list = [file for file in reg_file_list if file.endswith(
    '.txt') and file.startswith("3")]
chip3_reg_list = sorted(chip3_reg_list, key=lambda x: (
    int(x.split('x')[1].split('.')[0]), x))


chip2_fits_list = [file for file in fits_list if file.endswith('2.fits')]
chip2_fits_list = sorted(chip2_fits_list, key=lambda x: (x.split('_')[1], x))

chip3_fits_list = [file for file in fits_list if file.endswith('3.fits')]
chip3_fits_list = sorted(chip3_fits_list, key=lambda x: (x.split('_')[1], x))

chip2_data, chip2_wcs, chip2_header = data_gen(chip2_fits_list)
chip2_dates = [np.ceil(h['EXPSTART'] - 50923) for h in chip2_header]

chip3_data, chip3_wcs, chip3_header = data_gen(chip3_fits_list)
chip3_dates = [np.ceil(h['EXPSTART'] - 50923) for h in chip3_header]


chip2_region_coord = region_coord(chip2_reg_list)
chip2_region_coord = np.array(chip2_region_coord)

chip3_region_coord = region_coord(chip3_reg_list)
chip3_region_coord = np.array(chip3_region_coord)

chip2_region_data = region_data(chip2_data, chip2_region_coord, chip2_dates)
chip3_region_data = region_data(chip3_data, chip3_region_coord, chip2_dates)

chip2_pairs = pair(chip2_region_data)
chip3_pairs = pair(chip3_region_data)

chip2sums = [(np.sum(c1), np.sum(c2)) for n1, c1, c2 in chip2_pairs]
chip2_background_adjusted_data = [(n1, (np.mean(c1) -
                                   background(c1, c2))/pixel_area) for n1, c1, c2 in chip2_pairs]
chip3_background_adjusted_data = [(n1, (np.mean(c1) -
                                   background(c1, c2))/pixel_area) for n1, c1, c2 in chip3_pairs]

chip2_app_mag = [(n1, (np.log(1000/count) + 22.57))
                 for n1, count in chip2_background_adjusted_data]
chip3_app_mag = [(n1, (np.log(1000/count) + 22.57))
                 for n1, count in chip3_background_adjusted_data]


# =============================================================================
# lombscargle
# =============================================================================
