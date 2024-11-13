from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Get the absolute path to the project directory
project_dir = Path(__file__).resolve().parent.parent

# Create a directory for results if it doesn't exist
result_dir = project_dir / '1_samples_properties'

# Read the CSV file
otu_table = pd.read_csv(project_dir / '1_dada2' /
                        'otu_table_treeholes_dada2_renamed_clean_red_t.tsv',
                        sep='\t')
otu_table = otu_table.rename(columns={'Unnamed: 0': 'sampleid'})
print(f'otu_table shape: {otu_table.shape}')

# Read the TSV file
metadata = pd.read_csv(project_dir / '1_samples_properties' /
                       'Metadata_singleTree.csv', sep='\t')
metadata = metadata.rename(columns={'ID.Sample': 'sampleid'})
print(f'metadata shape: {metadata.shape}')

# Merge the two dataframes
merged_data = pd.merge(metadata, otu_table,
                       left_on='sampleid', right_on='sampleid')
print(f'merged_data shape: {merged_data.shape}')

# Calculate the total number of reads per sample
merged_data['total_reads'] = merged_data.iloc[:, 35:].sum(axis=1)

# Create a new dataframe with the required information
plot_data = merged_data[['Date', 'Location', 'total_reads']].copy()

# Convert Date to datetime
plot_data['Date'] = pd.to_datetime(plot_data['Date'], format='%d.%m.%y')

# Sort the dataframe by date
plot_data = plot_data.sort_values('Date')


def create_plot(x, evenly_spaced=False):

    fig, ax = plt.subplots(figsize=(16, 10))

    min_reads = plot_data['total_reads'].min() - 10

    if evenly_spaced:
        unique_dates = plot_data['Date'].unique()
        num_dates = len(unique_dates)
        
        # Create evenly spaced x-coordinates
        x_coords = np.arange(num_dates)
        
        # Set the x-ticks to these evenly spaced coordinates
        plt.xticks(x_coords, [d.strftime('%Y-%m-%d') for d in unique_dates],
                   rotation=90, ha='center', fontsize=8)
        
        # Remap the x-values for the scatter plot
        date_to_coord = dict(zip(unique_dates, x_coords))
        x_values = np.array([date_to_coord[date] for date in plot_data['Date']])
        
        scatter = ax.scatter(x=x_values, y=plot_data['Location'],
                             s=(plot_data['total_reads'] - min_reads) / 10,
                             c=plot_data['total_reads'], cmap='viridis',
                             alpha=0.6)
    else:
        scatter = ax.scatter(x=plot_data['Date'], y=plot_data['Location'],
                             s=(plot_data['total_reads'] - min_reads) / 10,
                             c=plot_data['total_reads'], cmap='viridis',
                             alpha=0.6)
        ax.xaxis.set_major_locator(mdates.DayLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d'))
        plt.xticks(plot_data['Date'], rotation=90, ha='center', fontsize=6)

    plt.title('Number of Reads per Tree Hole over Time', fontsize=16)
    plt.xlabel('Date', fontsize=12)
    plt.ylabel('Tree Hole', fontsize=12)

    cbar = plt.colorbar(scatter)
    cbar.set_label('Number of Reads', fontsize=12)

    plt.tight_layout()
    return fig


# Create and save the time-based plot
fig_time = create_plot('Date', evenly_spaced=False)
fig_time.savefig(result_dir / 'reads_per_treehole_over_time.png',
                 bbox_inches='tight', dpi=300)
plt.close(fig_time)

# Create and save the evenly-spaced plot
fig_even = create_plot('Date', evenly_spaced=True)
fig_even.savefig(result_dir / 'reads_per_treehole_evenly_spaced.png',
                 bbox_inches='tight', dpi=300)
plt.close(fig_even)

print("Two figures have been created and saved in the 'results' directory:")
print("1. 'reads_per_treehole_over_time.png' (time-based x-axis with all dates)")
print("2. 'reads_per_treehole_evenly_spaced.png' (evenly spaced x-axis)")
