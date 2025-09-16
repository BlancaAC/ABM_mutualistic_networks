#!/usr/bin/env python3
"""
Main script for running plant-pollinator interaction simulations.
"""

# Standard libraries
import os
import warnings

# Third-party libraries
import numpy as np
import pandas as pd
from networkx.convert_matrix import to_pandas_adjacency

# Custom modules
from classes import Environment_plant, Environment_pol
from data_analysis import area_plot
from modelling import (
    initial_network,
    update_totalinks,
    remove_zero,
    degree_dist,
    initial_pollinators_random,
)
from visualization import plot_abundances, plot_priors
import nxproperties as nxp

# Ignore all warnings (use with caution)
warnings.filterwarnings('ignore')


# Directory where output files will be saved
output_dir = "Slice_analysisv2"
os.makedirs(output_dir, exist_ok=True)

# --- Simulation parameters ---

# Number of draws from the prior distribution
n_draws = 100_000

# Plot and time slice selection
selected_plot = 'A'
slice_index = 3
interaction_scale = 5  # Number of pollinators per plant

# Construct the output filename
output_filename = f"{selected_plot}{slice_index}_nd{n_draws}"

# --- Load interaction network for the selected plot and slice ---

network_path = f"Data/Temporal_slices_individuals/network_plot_{selected_plot}_slice_{slice_index}.csv"
df_net = pd.read_csv(network_path, sep=';')

# Extract plant IDs as strings
plant_list = df_net["Plant_id"].astype(str).tolist()

# --- Load plant coordinates and filter for current plant list ---

df_plant = pd.read_csv("Data/coords_plot_month.csv", sep=';')
df_plant_pm = df_plant.query("Plant_id in @plant_list").drop_duplicates(subset="Plant_id")

# --- Load total number of visits for the selected plot and slice ---

df_visit = pd.read_csv("Data/n_visits_slice_plot.csv", sep=';')
nt_links = df_visit.query("Slice == @slice_index and Plot == @selected_plot")["Frequency"].values[0]

# --- Calculate spatial boundaries from plant data ---
xmin, ymin, xmax, ymax, diag = area_plot(df_plant_pm)

# --- Initialize plant environment ---
# Creates plant agents with spatial positions (randomly distributed within bounds),
# species identity, and unique IDs

envp = Environment_plant(
    df_plant_pm,
    random_position=True,
    xmin=xmin,
    xmax=xmax,
    ymin=ymin,
    ymax=ymax
)

#envp = Environment_plant(df_plant_pm)
#envp = Environment_plant(df_plant_pm, regular_position=True, xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax)

# Number of plant agents
n_plants = len(envp.plant_list)

# Number of pollinator agents (based on scale)
n_pollinators = interaction_scale * n_plants

# --- Compute relative abundances of pollinator species ---

# Proportional abundance per species (excluding the first column, assumed to be Plant_id)
dist_pol = df_net.iloc[:, 1:].sum()
dist_pol = dist_pol / dist_pol.sum()
dist_pol = np.round(dist_pol, decimals=4)

# Sort species by abundance to ensure last adjustment is on the largest value
dist_pol = dist_pol.sort_values()

# Adjust final value to ensure exact sum to 1 (fix for rounding error)
dist_pol.iloc[-1] += (1.0 - dist_pol.sum())

# --- Count number of unique species ---

n_plant_species = df_plant_pm["Plant_sp_complete"].nunique()
n_pollinator_species = len(dist_pol.index)

# --- Report summary statistics ---

print(f"Number of plant agents: {n_plants}")
print(f"Number of pollinator agents: {n_pollinators}")
print(f"Number of plant species: {n_plant_species}")
print(f"Number of pollinator species: {n_pollinator_species}")

# --- Save pollinator abundance distribution ---

dist_pol_path = os.path.join(output_dir, f"{output_filename}_species_dist_pol_model2.csv")
dist_pol.to_csv(dist_pol_path)

# Optional: plot and save species abundance distribution
# plot_abundances(dist_pol, output_filename, save=True)

# --- Step 4: Generate priors (r distributions for specialists and generalists) ---

# Parameters for gamma distributions
shape_specialist = 2
shape_generalist = 2
scale_specialist = 2
scale_generalist = 2

# Name suffix for saving prior plots or files
prior_name = f"_shapes_eg_{shape_specialist}_{shape_generalist}_scales{scale_specialist}_{scale_generalist}"

# Generate gamma-distributed priors
prior_specialist = pd.Series(np.random.gamma(shape_specialist, scale_specialist, size=n_draws))
# prior_generalist = pd.Series(np.random.gamma(shape_generalist, scale_generalist, size=n_draws))

# Optional: alternative exponential distribution priors
# exp_specialist = 10
# exp_generalist = 20
# prior_name = f"_scales_exp_{exp_specialist}_{exp_generalist}"
# prior_specialist = pd.Series(np.random.exponential(scale=exp_specialist, size=n_draws))
# prior_generalist = pd.Series(np.random.exponential(scale=exp_generalist, size=n_draws))

# Optional: plot prior distributions
# plot_priors(prior_specialist, prior_generalist, prior_name)

# Option to save ABM matrices for each iteration
save_abm_matrices = False  # Set to True to enable saving

# --- Step 5: Run simulation over multiple draws ---

all_deg_pol = pd.DataFrame()
all_deg_pla = pd.DataFrame()

for i in range(n_draws):
    
    # Generate initial pollinators with random distribution
    generalists, df_polPM = initial_pollinators_random(
        dist_pol, n_pollinator_species, n_pollinators,
        xmin, xmax, ymin, ymax
    )
    
    #generalistas,df_polPM = initial_pollinators(dist_pol, n_pollinator_species, n_pollinators, xmin, xmax, ymin, ymax)
    #generalistas,df_polPM = initial_pollinators_random(dist_pol, n_pollinator_species, n_pollinators, xmin, xmax, ymin, ymax,random_distribution=True)

    # Assign interaction radius from prior (same for both types here)
    df_polPM.loc[df_polPM['Tipo'] == 'Especialista', 'Radius'] = prior_specialist[i]
    df_polPM.loc[df_polPM['Tipo'] == 'Generalista', 'Radius'] = prior_specialist[i]
    #df_polPM.loc[df_polPM['Tipo'] == 'Generalista', 'Radius'] = prior_generalist[i]

    # Initialize pollinator environment
    envpol = Environment_pol(df_polPM)

    # Create initial bipartite network (pollinator × plant)
    B = initial_network(df_polPM['Pol_id'].tolist(), df_plant_pm['Plant_id'].tolist())

    # Assign interactions based on distance constraints
    update_totalinks(nt_links, envpol, envp, B, xmin, xmax, ymin, ymax)

    # Remove agents with no links
    remove_zero(B)

    # Extract plant-pollinator submatrix
    B_df = to_pandas_adjacency(B)
    df_ABM = B_df.iloc[:-len(df_plant_pm), len(df_polPM):]  # [pollinators, plants]

    # Annotate ABM matrix with species
    df_ABM['Sp_pol'] = df_ABM.index.map(df_polPM.set_index('Pol_id')['Specie'])
    df_ABM.loc['Sp_plant'] = df_ABM.columns.map(df_plant_pm.set_index('Plant_id')['Plant_sp'])

    # Optionally save ABM matrix for this iteration
    if save_abm_matrices:
        abm_path = os.path.join(output_dir, f"{i}_{output_filename}_{prior_name}_ABM.csv")
        df_ABM.to_csv(abm_path)

    # Aggregate pollinator links by species
    df_species = df_ABM.groupby('Sp_pol').sum()

    # Transpose to get plant × pollinator species matrix (model)
    df_model = df_species.T

    # Real interaction matrix (based on input network)
    df_real = df_net.set_index('Plant_id')

    # --- Degree distributions ---

    pol_d_model, plant_d_model = degree_dist(df_model)
    pol_d_real, plant_d_real = degree_dist(df_real)

    # Merge plant degree distributions
    plant_d_model.index = plant_d_model.index.map(str)
    plant_d_real.index = plant_d_real.index.map(str)
    index_union = plant_d_real.index.union(plant_d_model.index)

    deg_plant = pd.DataFrame(index=index_union)
    deg_plant['Model'] = plant_d_model
    deg_plant['Real'] = plant_d_real
    deg_plant = deg_plant.fillna(0)

    # Merge pollinator degree distributions
    deg_pol = pd.concat([pol_d_model, pol_d_real], axis=1, keys=['Model', 'Real'])
    deg_pol = deg_pol.fillna(0)

    # Add metadata
    deg_pol['iteration'] = i
    deg_plant['iteration'] = i
    deg_pol['r_esp'] = prior_specialist[i]
    deg_pol['r_gen'] = prior_specialist[i]
    #deg_pol['r_gen'] = prior_generalist[i]

    #box_plot_comp(deg_pol, deg_pla, 'Degree')
    #KDE_plot_comp(deg_pol, deg_pla, 'Degree')

    # Add pollinator type
    deg_pol['Tipo'] = deg_pol.index.map(lambda j: 'Generalista' if j in generalists else 'Especialista')

    # Accumulate results
    all_deg_pol = pd.concat([all_deg_pol, deg_pol])
    all_deg_pla = pd.concat([all_deg_pla, deg_plant])

# --- Save aggregated degree distributions ---

deg_pol_path = os.path.join(output_dir, f"{output_filename}_{prior_name}_species_deg_pol_model2.csv")
deg_pla_path = os.path.join(output_dir, f"{output_filename}_{prior_name}_species_deg_pla_model2.csv")

all_deg_pol.to_csv(deg_pol_path)
all_deg_pla.to_csv(deg_pla_path)




    