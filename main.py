# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 21:55:23 2023

@author: chend
"""

import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.io.fits import getdata
from astropy.wcs import WCS
import astropy.units as u
from astropy.visualization import ZScaleInterval, ImageNormalize
from regions import CirclePixelRegion, PixCoord
import numpy as np
import os
from astropy.timeseries import LombScargle
from scipy.fft import fft, fftfreq
#from scipy.signal import lombscargle
pixel_area = 0.099595926**2
expected_periods = [20.5, 17.6, 17.0, 21.4, 16.0, 15.2,
                    10.4, 11.7, 15.2, 17.3, 31.4, 18.2, 12.3, 13.6, 13.7]


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

    previous_region_data = np.array([0])
    region_data = np.array([1])
    regional_flux_list = []
    for image in range(reg_data.shape[0]):
        image_regions = region_coord[image]

        current_data_list = reg_data[image]
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
                    mode='subpixels').cutout(current_data_list)

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
chip2_dates = [(h['EXPSTART'] - 50923) for h in chip2_header]

chip3_data, chip3_wcs, chip3_header = data_gen(chip3_fits_list)
chip3_dates = [(h['EXPSTART'] - 50923) for h in chip3_header]


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

chip2_app_mag = [list((n1, np.log(1000/count) + 22.57))
                 for n1, count in chip2_background_adjusted_data]
chip3_app_mag = [list((n1, np.log(1000/count) + 22.57))
                 for n1, count in chip3_background_adjusted_data]

star_list = []
total_star_list = []


def splitter(app_mag_array, star_count):
    star_list = []
    total_star_list = []
    for n in range(star_count):
        star_list = [i[1] for i in app_mag_array[n::star_count]]
        total_star_list.append(star_list)
    return total_star_list


number = 3
app_mag_list = np.array(splitter(chip2_app_mag, 7) +
                        splitter(chip3_app_mag, 8))
mean_app_mag_list = [(np.max(star_list) +
                     np.min(star_list))/2 for star_list in app_mag_list]

dates2 = chip2_dates[:-1]
dates3 = chip3_dates[:-1]
calculated_period_list = []
for n, i in enumerate(app_mag_list):
    time = dates2 if n < 7 else dates3
    i = i[:-1]
    magnitude = i
# =============================================================================
#     plt.plot(time, magnitude)
#     plt.xlabel('Time (days)')
#     plt.ylabel('Magnitude')
#     plt.title('Sine Wave')
#     plt.show()
#
# =============================================================================
    frequency, power = LombScargle(time, magnitude).autopower()
    top_n_indices = np.argsort(power)[-200:]
# =============================================================================
#
#     # Plot the periodogram
#     plt.plot(1 / frequency, power)
#     plt.xlabel('Period (days)')
#     plt.ylabel('Lomb-Scargle Power')
#     plt.show()
#
# =============================================================================

    period = 1 / frequency[top_n_indices]

    best_period_index = np.argmin(abs(expected_periods[n] - period))
    best_period = period[best_period_index]

    calculated_period_list.append(best_period)

    print(
        f'Star: {n+1} Calculated period: {best_period:2.3f} days  '
        f'Expected period: {expected_periods[n]} days  '
        f'Difference: {abs(best_period-expected_periods[n]):2.3f}')


absolute_magnitude_list = np.array([(-3.07*((np.log(P)/np.log(10)) - 1) - 4.89)
                                   for P in calculated_period_list])

distance_list = 10**((mean_app_mag_list - absolute_magnitude_list+5)/5)

plt.plot((np.log(calculated_period_list)/np.log(10)),
         absolute_magnitude_list)
plt.gca().invert_yaxis()
plt.scatter((np.log(calculated_period_list)/np.log(10)),
            absolute_magnitude_list, color='r')
plt.xlabel("Log(P)")
plt.ylabel('M')
plt.title("Absolute Magnitude-Period Graph for Cepheid Variables")

print(
    f"Mean Distance to galaxy = {(np.mean(distance_list)/10**6):3.2f} Megaparsec \n"
    f"Associated Hubble Constant of {(440/(np.mean(distance_list)/10**6)):3.2f} Km/s")
# =============================================================================
# absolute_magnitude_list -> absolute magnitude for all 15 stars.
# app_mag_list -> apparent magnitude for all 15 stars on each 12 days.
# mean_app_mag_list -> weighted mean apparent magnitude for all 15 stars.
#
# calculated_period_list -> Period for all 15 stars.
# chip2/3_data -> 800x800 data grid for both chips on each 12 days.
# chip2/3_region_data -> data grid of cepheid regions and associated background
#       region
# chip2/3_background_adjusted_data -> data grid of cepheid regions accounting
#       for background
#
# =============================================================================
