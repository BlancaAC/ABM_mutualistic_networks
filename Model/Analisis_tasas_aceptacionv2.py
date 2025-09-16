# Análisis de las tasas de aceptaciones

import numpy as np
import pandas as pd

import os
import matplotlib.pyplot as plt

def plot_acceptance_rates(rate_data, labels, title="Comparison of Acceptance Rates Across Models",
                          save=False, filename="acceptance_rates.png", output_dir="."):
    """
    Plot a boxplot comparing acceptance rates across different model configurations.

    Parameters:
        rate_data (list of array-like): List of arrays or Series with acceptance rate values.
        labels (list of str): List of labels corresponding to each model.
        title (str): Title of the plot.
        save (bool): Whether to save the plot as a PNG file (default: False).
        filename (str): Name of the output file (used only if save=True).
        output_dir (str): Directory where the figure will be saved (default: current directory).

    Returns:
        None
    """
    plt.figure(figsize=(10, 6))
    plt.boxplot(rate_data, labels=labels)
    plt.title(title)
    plt.ylabel("Acceptance Rate")
    plt.xticks(rotation=15)
    plt.tight_layout()

    if save:
        os.makedirs(output_dir, exist_ok=True)
        path = os.path.join(output_dir, filename)
        plt.savefig(path, dpi=300)
        plt.close()
        print(f"Plot saved to: {path}")
    else:
        plt.show()


# Load data
df1 = pd.read_csv("resultados_tasa_B1_nd100000__shapes_eg_2_2_scales2_2random_deg_pol_randomproof.csv")
df2 = pd.read_csv("resultados_tasa_D1_nd100000__shapes_eg_2_2_scales2_2species_deg_pol_model2.csv")
df3 = pd.read_csv("resultados_tasa_B1_nd100000__shapes_eg_2_2_scales5_1random_deg_pol.csv")
df4 = pd.read_csv("resultados_tasa_A3_nd100000__shapes_eg_2_2_scales5_1regular_deg_pol.csv")
df5 = pd.read_csv("resultados_tasa_C2_nd100000__shapes_eg_2_2_scales5_1XY_deg_pol.csv")

rate_data = [df1['tasa'], df2['tasa'], df3['tasa'], df4['tasa'], df5['tasa']]
labels = ["All Random", "Real_D", "Random", "Regular", "XY"]

# Show 
plot_acceptance_rates(rate_data, labels)

# Save
plot_acceptance_rates(rate_data, labels, save=True, filename="acceptance_rates_comparison.png", output_dir="figures")

