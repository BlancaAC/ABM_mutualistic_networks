#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 11:59:45 2024

@author: galeanojav
"""

# visualization.py

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns


def plot_abundances(dist_pol, name_file, bins=30, save=False):
    """
    Plot a histogram of species abundances with high quality settings suitable for publication.
    
    Parameters:
    - dist_pol: DataFrame or Series containing abundance data.
    - name_file: Title for the plot and filename if saved.
    - bins: Number of bins in the histogram (default 30).
    - save: Boolean, if True, saves the plot to a file.
    """
    # Set the style of the plot to 'whitegrid' for better visibility
    sns.set(style="whitegrid")
    
    # Create the figure with a specific size and high DPI for publication quality
    plt.figure(figsize=(10, 6), dpi=300)
    
    # Create histogram
    dist_pol.hist(bins=bins, color='skyblue', edgecolor='black')
    plt.xlabel('P(Abundancies)', fontsize=12)
    plt.ylabel('Number of Species', fontsize=12)
    plt.title(name_file, fontsize=14)
    
    
    # Save the plot to a file if requested, with high resolution
    if save:
        plt.savefig(f"{name_file}.png", format='png', dpi=300, bbox_inches='tight')  # Saves the plot as a PNG file
        
    # Show the plot
    plt.show()
    
    # Clear the figure to free memory and avoid interference with other plots
    plt.clf()
    plt.close()
    
def plot_priors(prior_e,prior_g,name_priors):
    prior_e.hist(alpha=0.5, bins=100,label='prior_esp', density=True)
    prior_g.hist(alpha=0.5, bins=100,label='prior_esp', density=True)
    plt.legend()
    plt.title('Priors gamma'+name_priors)
    plt.show()
    


def plot_agents(env_plant, env_pol, plot, month):
    """Plot the locations of plant and pollinator agents in the environment.

    Parameters:
    env_plant (Environment_plant): The environment containing the plant agents.
    env_pol (Environment_pol): The environment containing the pollinator agents.

    """
    plt.figure(figsize=(10, 10))

    for plant in env_plant.plant_list:
        plt.scatter(plant.x, plant.y, color='green', label='Plants')

    for pol in env_pol.pol_list:
        plt.scatter(pol.x, pol.y, color='red', label='Pollinators')

    # Remove duplicates in legend
    handles, labels = plt.gca().get_legend_handles_labels()
    by_label = dict(zip(labels, handles))
    plt.legend(by_label.values(), by_label.keys())

    plt.xlim([xmin, xmax])
    plt.ylim([ymin, ymax])
    plt.title(f'Plant and Pollinator Locations: Plot {plot} month {month}')
    plt.xlabel('X coordinate')
    plt.ylabel('Y coordinate')
    plt.show()
    


def comparative_plot(deg_pol, deg_pla, text, plot_function, **kwargs):
    fig, axes = plt.subplots(1, 2)
    axes[0].set_title(text + ' of Pollinators')
    plot_function(ax=axes[0], data=deg_pol, **kwargs)
    axes[1].set_title(text + ' of Plants')
    plot_function(ax=axes[1], data=deg_pla, **kwargs)
    plt.show()

def box_plot_comp(deg_pol, deg_pla, text):
    comparative_plot(deg_pol, deg_pla, text, sns.boxplot)

def KDE_plot_comp(deg_pol, deg_pla, text):
    comparative_plot(deg_pol, deg_pla, text, sns.kdeplot)

def Cumulative_comp(deg_pol, deg_pla, text):
    comparative_plot(deg_pol, deg_pla, text, sns.histplot, bins=50, stat="density", 
                     element="step", fill=False, cumulative=True, common_norm=False)

def LogLog_comp(deg_pol, deg_pla, text):
    # Apply log transformation to data
    log_deg_pol = np.log(deg_pol)
    log_deg_pla = np.log(deg_pla)
    # Use the comparative_plot function to plot the log-log distributions
    comparative_plot(log_deg_pol, log_deg_pla, text, sns.histplot, bins=50, stat="density", 
                     element="step", fill=False, common_norm=False)

def LogLog_comp_scatter(deg_pol, deg_pla, text):
    # Apply log transformation to data
    log_deg_pol = np.log(deg_pol)
    log_deg_pla = np.log(deg_pla)

    # Use the comparative_plot function to plot the log-log distributions
    comparative_plot(log_deg_pol, log_deg_pla, text, sns.scatterplot)

