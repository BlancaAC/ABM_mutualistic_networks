"""
Script to analyze and visualize the mean absolute error across different distance metrics
for various models based on summary CSVs (e.g., 'abserrors_*').
"""

import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import os


def mean_metric(df, metric):
    """
    Compute the mean absolute error for a specific metric.

    Parameters:
        df (pandas.DataFrame): DataFrame with metric results.
        metric (str): Short label of the metric (e.g., 'MAE', 'KL').

    Returns:
        float: Mean value of 'total_abs_error' for the given metric.
    """
    metric_df = df[df['metric'] == metric]
    return metric_df['abs_error'].mean()


def analyze_absolute_errors(
    file_pattern='abserrors_*',
    output_csv='resultados_metricas_v2.csv',
    output_plot='abserror_by_metrics.png',
    save_plot=True
):
    """
    Analyze and plot mean absolute error per metric across models.

    Parameters:
        file_pattern (str): Glob pattern to match input CSV files.
        output_csv (str): Output CSV file with summarized results.
        output_plot (str): Output filename for the barplot.
        save_plot (bool): Whether to save the plot (True) or display it (False).

    Returns:
        pandas.DataFrame: DataFrame with summarized mean errors.
    """
    # Define the set of metric codes
    metrics = ['HELL', 'JS', 'KL', 'KS', 'MAE', 'RMAE', 'WASS']

    results = []

    # Find all matching files
    file_list = glob.glob(file_pattern)

    for file_name in file_list:
        model_name = os.path.splitext(os.path.basename(file_name))[0].split('_')[-1]
        try:
            data = pd.read_csv(file_name)
        except Exception as e:
            print(f"Error reading {file_name}: {e}")
            continue

        for metric in metrics:
            value = mean_metric(data, metric)
            results.append({
                'Metric': metric,
                'Model': model_name,
                'Mean_abs_error': value
            })

    # Create summary DataFrame and save
    df_results = pd.DataFrame(results)
    df_results.to_csv(output_csv, index=False)

    # Plot
    plt.figure(figsize=(12, 8))
    sns.barplot(x='Metric', y='Mean_abs_error', hue='Model', data=df_results)
    plt.title('Mean Absolute Error Across Metrics by Model')
    plt.ylabel('Mean Absolute Error')
    plt.xlabel('Metric')
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()

    if save_plot:
        plt.savefig(output_plot, dpi=300)
        plt.close()
        print(f"Plot saved to: {output_plot}")
    else:
        plt.show()

    return df_results

if __name__ == "__main__":
    analyze_absolute_errors(
        file_pattern='abserrors_*',
        output_csv='resultados_metricas_model2.csv',
        output_plot='Figures/abserror_by_metrics.png',
        save_plot=True
    )