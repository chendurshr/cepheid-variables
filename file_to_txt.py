# -*- coding: utf-8 -*-
"""
Created on Mon Dec 11 21:19:02 2023

@author: chend
"""
import os
import shutil


def convert_files(input_folder, output_folder):
    # Ensure the output folder exists
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Get a list of all files without extensions in the input folder
    files_to_convert = [f for f in os.listdir(
        input_folder) if '.' not in f and f.startswith("3")]

    # Iterate over each file and convert it
    for file_name in files_to_convert:
        input_path = os.path.join(input_folder, file_name)
        output_path = os.path.join(output_folder, file_name + ".txt")

        # Copy the file and change the extension
        shutil.copyfile(input_path, output_path)

        print(f"Converted: {input_path} -> {output_path}")


# Specify the input and output folders
input_folder = r"C:\Users\chend\Desktop\Projects\CepheidVariables\reg"
output_folder = r"C:\Users\chend\Desktop\Projects\CepheidVariables\reg"

# Call the function to convert files
convert_files(input_folder, output_folder)
