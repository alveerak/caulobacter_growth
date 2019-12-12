from os.path import join
import glob
import os

import numpy as np
import pandas as pd

import bokeh.io
import holoviews as hv
import bokeh_catplot
import panel as pn

import skimage.feature
import skimage.filters
import skimage.filters.rank
import skimage.io
import skimage.morphology
import skimage.segmentation

import scipy.ndimage

import colorcet

import bebi103

def zero_crossing_filter(im, thresh):
    """
    Returns image with 1 if there is a zero crossing and 0 otherwise.

    thresh is the the minimal value of the gradient, as computed by Sobel
    filter, at crossing to count as a crossing.
    """
    # Square structuring element
    selem = skimage.morphology.square(3)

    # Do max filter and min filter
    im_max = scipy.ndimage.filters.maximum_filter(im, footprint=selem)
    im_min = scipy.ndimage.filters.minimum_filter(im, footprint=selem)

    # Compute gradients using Sobel filter
    im_grad = skimage.filters.sobel(im)

    # Return edges
    return (((im >= 0) & (im_min < 0)) | ((im <= 0) & (im_max > 0))) & (
        im_grad >= thresh
    )

def mark_growth_events(filename):
    # Read in data
    df = pd.read_csv(filename)

    # Add empty column for growth event label
    df["Growth Event"] = ""

    # Add first growth event value
    label_cnt = 1
    df.iloc[0, df.columns.get_loc("Growth Event")] = label_cnt

    # Loop through remaining rows and compare difference in area to label growth event
    for index in range(1, len(df.index)):
        # Large difference in area signifies new growth event
        if (
            df.iloc[index - 1, df.columns.get_loc("area (µm²)")]
            - df.iloc[index, df.columns.get_loc("area (µm²)")]
            > 0.5
        ):
            # Update to new growth event
            label_cnt += 1
        # Add growth event to row
        df.iloc[index, df.columns.get_loc("Growth Event")] = label_cnt

    # Separate bacteria and reindex growth event for bacteria 2 from 1
    df_bac1 = df.loc[df["bacterium"] == 1]
    df_bac2 = df.loc[df["bacterium"] == 2]
    df_bac2['Growth Event'] = df_bac2['Growth Event'] - df_bac2['Growth Event'].min() + 1

    # Write modified data back to original dataframe
    df = df_bac1.append(df_bac2, ignore_index=True)

    # Write data to CSV
    output_file = join("..", "data", "hw6.1_caulobacter_growth_event_results.csv")
    df.to_csv(path_or_buf=output_file)
    return (df_bac1, df_bac2)

def get_time_diff(df, df_tdiff):
    """Return dataframe with time differences of cell growth"""
    # Loop through growth events
    for index in range(1, len(df_tdiff.index) + 1):
        # Find starting time of each growth event
        t1 = int(df.loc[df["Growth Event"] == index, ["time (min)"]].min())
        t2 = int(df.loc[df["Growth Event"] == index + 1, ["time (min)"]].min())

        # Take difference in time and append to empty column in row
        t_diff = t2 - t1
        df_tdiff.at[index - 1, "Time Difference (min)"] = t_diff

    return df_tdiff
