# Copyright (c) 2026 Flavia Araujo
# SPDX-License-Identifier: MIT

import csv
import os
from pathlib import PureWindowsPath
from datetime import datetime


# Global path to save the Lynx CSV files
LYNX_ROOT_PATH = PureWindowsPath('Y:\\MethodManager4\\Momentum_Input\\')
# LYNX_ROOT_PATH = ('/Users/flavia/PycharmProjects/growth_profile_watcher/output/') # used for testing


def create(plate_data, plate_info, target_OD, target_volume):
    """
    Create a Lynx CSV file to dispense volumes based on normalized OD values.
    :param plate_data: List of dictionaries containing plate growth data.
    :param plate_info: List of dictionaries containing plate_type, plate_id, num_rows, and num_columns.
    :return: Path to the created Lynx CSV file.
    """
    # Extract plate information
    plate_id = plate_info[0]['plate_id']
    num_rows = plate_info[0]['num_rows']
    num_columns = plate_info[0]['num_columns']
    num_wells = num_rows * num_columns
    num_rows = 8
    num_columns = 12

    main_plate_id = plate_id.split("_")[0]

    # Define output CSV file name
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    output_filename_sample = f"Sample_{main_plate_id}_{timestamp}.csv"
    output_filename_media = f"Media_{main_plate_id}_{timestamp}.csv"
    output_filepath_sample = os.path.join(LYNX_ROOT_PATH, output_filename_sample)
    output_filepath_media = os.path.join(LYNX_ROOT_PATH, output_filename_media)

    # Extract the last entry of plate data for OD values
    last_entry = plate_data[-1]
    wells_growth = last_entry['wells_growth']

    # Define parameters
    # target_od = 3 # add as function parameter
    # target_volume = 125 # add as function parameter
    max_dispense_volume = 300 # in µL
    min_dispense_volume = 10 # in µL

    # Lists to hold calculated volumes
    dispense_volumes_sample = []
    dispense_volumes_media = []

    for od in wells_growth:
        try:
            od_value = float(od)
            if od_value == 0:
                sample_volume = max_dispense_volume
            else:
                sample_volume = (target_volume / od_value) * target_OD
                sample_volume = max(min(sample_volume, max_dispense_volume), min_dispense_volume)

            media_volume = max_dispense_volume - sample_volume
            dispense_volumes_sample.append(sample_volume)
            dispense_volumes_media.append(media_volume)
        except ValueError:
            dispense_volumes_sample.append(max_dispense_volume)
            dispense_volumes_media.append(0.0)

    # --- Prepare and Write Sample Dispense CSV ---
    header = [f"VI;{num_columns};{num_rows}"] + [str(i) for i in range(1, num_columns + 1)]
    lynx_data_sample = [header]

    for row_index in range(num_rows):
        row_label = chr(65 + row_index)
        row_data = [row_label]
        for col in range(num_columns):
            if num_wells == 96 or (row_index % 2 == 0 and col % 2 == 0):  # Adjust for 24-well plate
                if num_wells == 96:
                    sample_volume = dispense_volumes_sample[row_index * num_columns + col]
                else:  # 24-well plate
                    sample_volume = dispense_volumes_sample[(row_index // 2) * (num_columns // 2) + (col // 2)]
            else:
                sample_volume = 0.0
            row_data.append(f"{sample_volume:.2f}")
        lynx_data_sample.append(row_data)

    with open(output_filepath_sample, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(lynx_data_sample)

    # --- Prepare and Write Media Dispense CSV ---
    lynx_data_media = [header]

    for row_index in range(num_rows):
        row_label = chr(65 + row_index)
        row_data = [row_label]
        for col in range(num_columns):
            if num_wells == 96 or (row_index % 2 == 0 and col % 2 == 0):
                if num_wells == 96:
                    media_volume = dispense_volumes_media[row_index * num_columns + col]
                else: # 24-well plate
                    media_volume = dispense_volumes_media[(row_index // 2) * (num_columns // 2) + (col // 2)]
            else:
                media_volume = 0.0
            row_data.append(f"{media_volume:.2f}")
        lynx_data_media.append(row_data)

    with open(output_filepath_media, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(lynx_data_media)

    print(f"Lynx Sample CSV file created at: {output_filepath_sample}")
    print(f"Lynx Media CSV file created at: {output_filepath_media}")

    return output_filename_sample, output_filename_media