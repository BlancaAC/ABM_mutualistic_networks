#Importing libraries

import glob
import os
import pandas as pd # Importing the pandas library for data manipulation
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import jensenshannon
from scipy.stats import ks_2samp, wasserstein_distance

###### Functions #################################

def hellinger_distance(p, q):
    """
    Compute the Hellinger distance between two probability distributions.

    Parameters:
        p (array-like): First probability distribution.
        q (array-like): Second probability distribution.

    Returns:
        float: Hellinger distance between p and q.
    """
    return np.sqrt(0.5 * np.sum((np.sqrt(p) - np.sqrt(q))**2))

def compute_metrics(df):
    """
    Compute a set of statistical distance metrics comparing 'Model' and 'Real' distributions
    for each simulation iteration in a DataFrame.

    Metrics include:
    - Jensen-Shannon divergence
    - Kolmogorov–Smirnov distance
    - Wasserstein distance
    - Hellinger distance
    - Relative Mean Absolute Error (MAE)
    - Absolute MAE
    - Kullback-Leibler divergence

    Parameters:
        df (pandas.DataFrame): Input DataFrame containing columns ['Model', 'Real', 'iteration', 'r_esp', 'r_gen'].

    Returns:
        pandas.DataFrame: A DataFrame where each row corresponds to a simulation iteration and its associated metrics.
    """
    results = []
    epsilon = 1e-10

    for iteration, group in df.groupby("iteration"):
        real = group["Real"].values
        model = group["Model"].values

        # Normalize distributions
        real_dist = real / real.sum() if real.sum() > 0 else np.ones_like(real) / len(real)
        model_dist = model / model.sum() if model.sum() > 0 else np.ones_like(model) / len(model)

        # Avoid log(0) and division by zero
        real_dist = np.clip(real_dist, epsilon, 1)
        model_dist = np.clip(model_dist, epsilon, 1)

        # Compute distance metrics
        js_div = jensenshannon(real_dist, model_dist, base=np.e) ** 2
        ks_stat, _ = ks_2samp(real, model)
        wass_dist = wasserstein_distance(real, model)
        hell_dist = hellinger_distance(real_dist, model_dist)
        rel_mae = np.mean(np.abs(model - real) / np.clip(real, epsilon, None))
        abs_mae = np.mean(np.abs(model - real))
        kl_div = np.sum(real_dist * np.log(real_dist / model_dist))

        r_esp = group["r_esp"].iloc[0]
        r_gen = group["r_gen"].iloc[0]

        results.append({
            "iteration": iteration,
            "JS_divergence": js_div,
            "KS_distance": ks_stat,
            "Wasserstein_distance": wass_dist,
            "Hellinger_distance": hell_dist,
            "Relative_MAE": rel_mae,
            "Absolute_MAE": abs_mae,
            "KL_divergence": kl_div,
            "r_esp": r_esp,
            "r_gen": r_gen
        })

    return pd.DataFrame(results)


def top_metric(df_metrics, metric, top_fraction=0.1):
    """
    Select the top-performing entries from a metrics DataFrame according to a specific metric.

    Parameters:
        df_metrics (pandas.DataFrame): DataFrame containing metric columns.
        metric (str): The metric name to sort by.
        top_fraction (float): Fraction of top entries to return (default is 0.1, i.e., top 10%).

    Returns:
        pandas.DataFrame: Top entries sorted by the given metric in descending order.
    """
    df_sorted = df_metrics.sort_values(by=metric, ascending=False)
    top_n = max(1, int(len(df_sorted) * top_fraction))
    return df_sorted.head(top_n)

def resume_error(df):
    """
    Compute error statistics per species in a simulation output DataFrame.

    Automatically detects if species are identified in the index or in a column
    (e.g., 'Unnamed: 0').

    Parameters:
        df (pandas.DataFrame): DataFrame with columns ['Model', 'Real'] and either a species index
                               or a column containing species names.

    Returns:
        pandas.DataFrame: Summary statistics including error and z-scores per species.
    """
    df = df.copy()

    # Detect species identifier column
    if "Unnamed: 0" in df.columns:
        species_column = "Unnamed: 0"
    elif df.index.name:
        df = df.reset_index()
        species_column = df.columns[0]
    else:
        raise ValueError("Could not determine the species identifier in the DataFrame.")

    # Group and calculate statistics
    grouped = df.groupby(species_column).agg(
        real_value=('Real', 'mean'),
        model_mean=('Model', 'mean'),
        model_std=('Model', 'std'),
        model_min=('Model', 'min'),
        model_max=('Model', 'max')
    ).reset_index()

    grouped['abs_error'] = np.abs(grouped['real_value'] - grouped['model_mean'])

    grouped['z_score'] = np.where(
        grouped['model_std'] > 0,
        grouped['abs_error'] / grouped['model_std'],
        np.nan
    )

    grouped['model_range'] = grouped['model_max'] - grouped['model_min']
    grouped['z_range'] = np.where(
        grouped['model_range'] > 0,
        grouped['abs_error'] / grouped['model_range'],
        np.nan
    )

    grouped = grouped.rename(columns={species_column: 'species'})
    return grouped

def resume_errorv2(df):
    """
    Compute error statistics for each species in a simulation output DataFrame.

    For each species (identified by the index column, usually representing species names),
    this function calculates:
    - Mean real value (assumed constant per species)
    - Mean, standard deviation, min and max of model predictions
    - Absolute error between model mean and real value
    - Traditional z-score (abs_error / std)
    - Z-score based on value range (abs_error / (max - min))

    Parameters:
        df (pandas.DataFrame): Input DataFrame with columns ['Real', 'Model'], and index or a column
                               identifying species (typically unnamed after CSV import).

    Returns:
        pandas.DataFrame: DataFrame with summary statistics and error measures per species.
    """
    df = df.copy()

    # Attempt to identify the species column (e.g., first unnamed column from CSV)
    #if "Unnamed: 0" in df.columns:
    #    species_column = "Unnamed: 0"
    #else:
    #    raise ValueError("Expected a column 'Unnamed: 0' containing species identifiers.")

    # Aggregate statistics by species
    grouped = df.groupby("Unnamed: 0").agg(
        real_value=('Real', 'mean'),  # debería ser único por especie
        model_mean=('Model', 'mean'),
        model_std =('Model', 'std'),
        model_min=('Model', 'min'),
        model_max=('Model', 'max')
    ).reset_index()

    # Compute absolute error
    grouped['abs_error'] = np.abs(grouped['real_value'] - grouped['model_mean'])

    # Compute traditional z-score (avoid division by zero)
    grouped['z_score'] = np.where(
        grouped['model_std'] > 0,
        grouped['abs_error'] / grouped['model_std'],
        np.nan
    )

    # Compute range-based z-score and range
    grouped['model_range'] = grouped['model_max'] - grouped['model_min']
    grouped['z_range'] = np.where(
        grouped['model_range'] > 0,
        grouped['abs_error'] / grouped['model_range'],
        np.nan
    )

    # Rename species column for clarity
    #grouped = grouped.rename(columns={species_column: 'species'})
    grouped = grouped.rename(columns={'Unnamed: 0': 'especie'})

    return grouped


def plot_degree_comp(data, plot_label):
    """
    Plot the comparison between real and simulated (model) degrees per species,
    including error bands and summary error metrics (Absolute Error, RMSE).

    Parameters:
        data (pandas.DataFrame): DataFrame containing columns 'Model', 'Real', 'iteration',
                                 and a species identifier column (typically 'Unnamed: 0').
        plot_label (str): Label to include in the plot title and output filename.

    Returns:
        float: Total absolute error across all species.
    """

    # Check required columns
    required_cols = {'Model', 'Real', 'iteration', 'Unnamed: 0'}
    if not required_cols.issubset(data.columns):
        raise ValueError(f"DataFrame must include the columns: {required_cols}")

    # Group data by species to compute summary statistics
    summary = data.groupby('Unnamed: 0').agg(
        mean_model=('Model', 'mean'),
        std_model=('Model', 'std'),
        min_model=('Model', 'min'),
        max_model=('Model', 'max'),
        real_degree=('Real', 'first')  # Real values are constant across iterations
    ).reset_index()

    # Build matrix: rows = iterations, columns = species
    model_matrix = data.pivot_table(index='iteration', columns='Unnamed: 0', values='Model', aggfunc='mean')
    species = summary['Unnamed: 0'].values
    real = summary['real_degree'].values
    mean_model = summary['mean_model'].values

    # Compute errors
    abs_error = np.abs(mean_model - real)
    rmse = np.sqrt(np.mean((model_matrix[species].values - real) ** 2, axis=0))
    total_abs_error = abs_error.sum()

    # Plot
    plt.figure(figsize=(14, 6))
    plt.plot(species, real, marker='o', linestyle='-', label='Real', linewidth=2)
    plt.plot(species, mean_model, marker='s', linestyle='--', label='Model Mean', linewidth=1)
    plt.fill_between(species, summary['min_model'], summary['max_model'],
                     color='gray', alpha=0.3, label='Model Min-Max Range')

    # Error lines
    plt.plot(species, abs_error, linestyle='dotted', marker='x', color='red', label='Abs. Error')
    plt.plot(species, rmse, linestyle='dashdot', marker='d', color='purple', label='RMSE')

    # Title and formatting
    plot_title = f"{plot_label} — Total Abs. Error: {total_abs_error:.2f}"
    plt.title(plot_title)
    plt.xticks(rotation=90)
    plt.xlabel('Species')
    plt.ylabel('Degree')
    plt.legend()
    plt.tight_layout()

    # Optional: save figure
    # output_file = f"Top10_{plot_label}.png"
    # plt.savefig(output_file, dpi=300)
    # plt.close()

    return total_abs_error

import os
import numpy as np
import matplotlib.pyplot as plt

def plot_degree_comp(data, plot_label, save=False, output_dir="."):
    """
    Plot the comparison between real and simulated (model) degrees per species.
    
    The plot includes:
    - Real vs. mean simulated degree
    - Shaded area for min-max simulated values
    - Absolute error and RMSE per species
    
    Parameters:
        data (pandas.DataFrame): DataFrame with columns 'Model', 'Real', 'iteration',
                                 and species identifiers in column 'Unnamed: 0'.
        plot_label (str): Identifier for the plot (e.g., plot-month-model combination).
        save (bool): Whether to save the plot as a PNG file (default: False).
        output_dir (str): Directory where the plot will be saved (default: current directory).
    
    Returns:
        float: Total absolute error across all species.
    """
    # Group by species to compute summary statistics
    summary = data.groupby('Unnamed: 0').agg(
        mean_model=('Model', 'mean'),
        std_model=('Model', 'std'),
        min_model=('Model', 'min'),
        max_model=('Model', 'max'),
        real_degree=('Real', 'first')  # Assumes real value is constant across iterations
    ).reset_index()

    # Construct model matrix (rows: iterations, columns: species)
    model_matrix = data.pivot_table(
        index='iteration', columns='Unnamed: 0', values='Model', aggfunc='mean'
    )

    species = summary['Unnamed: 0'].values
    real = summary['real_degree'].values
    mean_model = summary['mean_model'].values

    # Compute errors
    abs_error = np.abs(mean_model - real)
    rmse = np.sqrt(np.mean((model_matrix[species].values - real) ** 2, axis=0))
    total_abs_error = abs_error.sum()

    # Plot
    plt.figure(figsize=(14, 6))
    plt.plot(species, real, marker='o', linestyle='-', label='Real', linewidth=2)
    plt.plot(species, mean_model, marker='s', linestyle='--', label='Model Mean', linewidth=1)
    plt.fill_between(
        species,
        summary['min_model'],
        summary['max_model'],
        color='gray',
        alpha=0.3,
        label='Model Min-Max Range'
    )
    plt.plot(species, abs_error, linestyle='dotted', marker='x', color='red', label='Abs. Error')
    plt.plot(species, rmse, linestyle='dashdot', marker='d', color='purple', label='RMSE')

    # Formatting
    plot_title = f"Top10_{plot_label} — Total Abs. Error: {total_abs_error:.2f}"
    plt.title(plot_title)
    plt.xticks(rotation=90)
    plt.xlabel('Species')
    plt.ylabel('Degree')
    plt.legend()
    plt.tight_layout()

    # Save plot if requested
    if save:
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, f"Top10_{plot_label}.png")
        plt.savefig(output_path, dpi=300)
        plt.close()
    else:
        plt.show()

    return total_abs_error



###### Main #################


# Define the pattern to find matching files
file_pattern = 'Slice_analysisv2/*species_deg_pol_model2.csv'


#file_pattern = 'Slice_analysis/*random_deg_pol_randomproof.csv'
#file_pattern = 'Slice_analysis/A*species_deg_pol_model2.csv'
#file_pattern = 'Slice_analysis/*regular_deg_pol.csv'
#file_pattern = 'Slice_analysis/*XY_deg_pol.csv'
#file_pattern = 'Slice_analysis/*random_deg_pol.csv'

# Find all files that match the pattern
file_list = glob.glob(file_pattern)
#print("Matched files:", file_list)

# Initialize a list to store individual DataFrames
abs_errors = []  # List to store absolute errors (or other results)

for file_path in file_list:
    print(f"Processing: {file_path}")

    # Extract filename from full path
    file_name = os.path.basename(file_path)

    # Extract metadata from filename (e.g., plot/month and model type)
    plot_month = file_name.split('_')[0]         # e.g., 'A3'
    model_type = file_name.split('_')[-1]        # e.g., 'model2.csv'

    # Load the degree distribution data
    data = pd.read_csv(file_path)
    
    # Compute distribution comparison metrics for the loaded data
    metrics_df = compute_metrics(data)

    # List of metrics to extract top-performing agents (or species)
    metrics = [
        "JS_divergence", 
        "KS_distance", 
        "Wasserstein_distance", 
        "Hellinger_distance", 
        "Relative_MAE", 
        "Absolute_MAE",
        "KL_divergence"
    ]

    # Get top-performing entries (e.g., best model fits) for each metric
    top_dfs = {
        metric: top_metric(metrics_df, metric, top_fraction=0.1)
        for metric in metrics
    }

    # Dictionary mapping metric names to shorter labels for plotting and filenames
    metric_labels = {
        'Absolute_MAE': 'MAE',
        'Relative_MAE': 'RMAE',
        'JS_divergence': 'JS',
        'KS_distance': 'KS',
        'Wasserstein_distance': 'WASS',
        'Hellinger_distance': 'HELL',
        'KL_divergence': 'KL'
    }

    # Iterate over each metric and its corresponding label
    for metric, label in metric_labels.items():
        if metric in top_dfs:
            # Build a name for the plot/output based on the metric, plot/month, and model
            name = f"{label}_{plot_month}_{model_type}"

            # Get the top-performing iterations for this metric
            top_iterations = top_dfs[metric]["iteration"].values

            # Filter original data to only those top iterations
            filtered_data = data[data["iteration"].isin(top_iterations)].copy()

            # Compute error summary for these filtered iterations
            summary = resume_error(filtered_data)
            summary["metric"] = label
            summary["plot"] = name

            # Accumulate summary in global list
            abs_errors.append(summary)

            # Optional plotting (commented out)
            #plot_degree_comp(data, plot_label="A3_model2", save=False) #Show Pictures in the screen
            plot_degree_comp(filtered_data, plot_label=name, save=True, output_dir="Figures")


# Combine all error summaries into a single DataFrame
df_abs_errors = pd.concat(abs_errors, ignore_index=True)

# Save to CSV
output_path = f"abserrors_zscore_{model_type}"
df_abs_errors.to_csv(output_path, index=False)