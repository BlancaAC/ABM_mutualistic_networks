#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 12:14:42 2024

@author: galeanojav
"""

import pandas as pd
import numpy as np

def freq_plant(plot: str, month: int, df_plant: pd.DataFrame) -> tuple:
    """Calculate frequency of plant species for a given plot and month.

    Parameters:
    plot (str): Plot to be analyzed.
    month (int): Month to be analyzed.
    df_plant (pd.DataFrame): DataFrame containing plant data.

    Returns:
    df_plant_pm (pd.DataFrame): Filtered DataFrame for specified plot and month.
    table_plants (pd.DataFrame): Frequency table of plant species.
    """

    df_plant_pm = df_plant.query('Plot == @plot and Month == @month')

    # Count the number of each species
    counts = df_plant_pm['Plant_sp_complete'].value_counts().rename('counts')

    # Calculate the frequency of each species
    freqs = counts / counts.sum()

    # Combine counts and frequencies into one DataFrame
    table_plants = pd.DataFrame({'counts': counts, 'Freq': freqs}).reset_index().rename(columns={'index': 'Plant_sp_complete'}).round(3)

    return df_plant_pm, table_plants

def area_plot(df_plantPM, margin=2):
    """Calculate minima, maxima and diagonal in our Plot."""
    
    min_vals = df_plantPM[['X', 'Y']].min() - margin
    max_vals = df_plantPM[['X', 'Y']].max() + margin
    
    dx, dy = max_vals - min_vals
    diag = np.hypot(dx, dy)
    
    return (*min_vals, *max_vals, diag)
