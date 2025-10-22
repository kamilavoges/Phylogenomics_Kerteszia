#!/usr/bin/env python3.9
"""
plot_ADMIXTURE_v0.py     runs make ADMIXTURE plots .  18oct2025

usage: 
# Basic usage
python plot_ADMIXTURE_v0.py --input admixture_results.txt

# Group by species with custom order
python plot_ADMIXTURE_v0.py --input admixture_results.txt --group species --order D,B,A

# Group by multiple columns, save as PDF
python plot_ADMIXTURE_v0.py --input admixture_results.txt --group species,population1 --graph pdf --output my_plot.pdf

# Add separation lines between groups
python plot_ADMIXTURE_v0.py --input admixture_results.txt --group species --group_sep

# Use population1 for bar labels and ggplot2-style plot
python plot_ADMIXTURE_v0.py --input admixture_results.txt --bar_label population1 --graph_style 2

# Use custom height and X-label size
python plot_ADMIXTURE_v0.py --input admixture_results.txt --graph_style 2 --height 6 --Xlabel_size 10
"""

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import os
from matplotlib import gridspec

def read_admixture_file(filename):
    """Read ADMIXTURE file using pandas with regex delimiter (like awk)"""
    # Use pandas with regex to split on any whitespace
    df = pd.read_csv(filename, sep=r'\s+', engine='python')
    return df

def get_component_columns(df):
    """Extract V1, V2, V3... component columns"""
    component_cols = [col for col in df.columns if re.match(r'^V\d+$', col)]
    # Sort by the number after V
    component_cols.sort(key=lambda x: int(x[1:]))
    return component_cols

def plot_admixture_1(df, component_cols, group_cols, group_order, output_format, output_file, input_file, group_sep, bar_label, xlabel_size):
    """Create ADMIXTURE stacked bar plot (original style)"""
    
    # Create a copy to avoid modifying original
    plot_df = df.copy()
    
    # Create grouping column
    if group_cols:
        # Verify all group columns exist
        missing_cols = [col for col in group_cols if col not in plot_df.columns]
        if missing_cols:
            print(f"Error: Group columns {missing_cols} not found in data.")
            print(f"Available columns: {list(plot_df.columns)}")
            return
        
        plot_df['_group'] = plot_df[group_cols].astype(str).agg('_'.join, axis=1)
        
        # Determine group order
        if group_order:
            # Use specified order for the first group column
            group_mapping = {group: i for i, group in enumerate(group_order)}
            plot_df['_group_order'] = plot_df[group_cols[0]].map(group_mapping)
            # Sort by the specified order and then by the group columns
            plot_df = plot_df.sort_values(['_group_order'] + group_cols)
        else:
            # Use alphabetical order
            plot_df = plot_df.sort_values(group_cols)
    
    # Prepare data for plotting
    components = plot_df[component_cols].values
    sample_names = plot_df[bar_label].values
    
    # Create plot
    fig, ax = plt.subplots(figsize=(max(12, len(sample_names) * 0.3), 8))
    
    # Create stacked bar plot
    bottom = np.zeros(len(sample_names))
    colors = plt.cm.Set3(np.linspace(0, 1, len(component_cols)))
    
    for i, col in enumerate(component_cols):
        ax.bar(range(len(sample_names)), components[:, i], bottom=bottom, 
               label=col, color=colors[i], width=1.0)
        bottom += components[:, i]
    
    # Add group separation lines if requested
    if group_sep and group_cols:
        # Get the first group column
        first_group_col = group_cols[0]
        
        # Find positions where group changes
        group_changes = []
        current_group = None
        
        for i, group_val in enumerate(plot_df[first_group_col]):
            if group_val != current_group:
                if current_group is not None:  # Not the first group
                    group_changes.append(i - 0.5)  # Position between bars
                current_group = group_val
        
        # Add vertical lines at group boundaries
        for pos in group_changes:
            ax.axvline(x=pos, color='black', linestyle=':', alpha=0.7, linewidth=1)
    
    # Customize plot
    ax.set_xlabel('Samples')
    ax.set_ylabel('Ancestry Proportion')
    ax.set_title('ADMIXTURE Analysis')
    
    # Set x-axis ticks with sample names but remove the tick marks
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right', fontsize=xlabel_size)
    # Remove x-axis tick marks
    ax.tick_params(axis='x', which='both', length=0)
    
    # Add legend
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    
    # Add grid for better readability
    ax.grid(True, axis='y', alpha=0.3)
    
    # Adjust layout
    plt.tight_layout()
    
    # Save or show plot
    if output_format != 'none':
        if not output_file:
            # Generate output file name from input file
            base_name = os.path.splitext(input_file)[0]
            output_file = f"{base_name}.{output_format}"
        else:
            # Check if output_file already has an extension
            if not any(output_file.endswith(ext) for ext in ['.png', '.svg', '.pdf']):
                output_file = f"{output_file}.{output_format}"
        
        plt.savefig(output_file, format=output_format, bbox_inches='tight')
        print(f"Plot saved as {output_file}")
    else:
        # Only show plot if not saving
        plt.show()

def plot_admixture_2(df, component_cols, group_cols, group_order, output_format, output_file, input_file, group_sep, bar_label, height, xlabel_size):
    """Create ggplot2-style ADMIXTURE plot with facets"""
    
    # Create a copy to avoid modifying original
    plot_df = df.copy()
    
    # Use first group column for faceting if available, otherwise use bar_label
    facet_col = group_cols[0] if group_cols else bar_label
    
    # Create grouping column for sorting
    if group_cols:
        # Verify all group columns exist
        missing_cols = [col for col in group_cols if col not in plot_df.columns]
        if missing_cols:
            print(f"Error: Group columns {missing_cols} not found in data.")
            print(f"Available columns: {list(plot_df.columns)}")
            return
        
        plot_df['_group'] = plot_df[group_cols].astype(str).agg('_'.join, axis=1)
        
        # Determine group order
        if group_order:
            # Use specified order for the first group column
            group_mapping = {group: i for i, group in enumerate(group_order)}
            plot_df['_group_order'] = plot_df[group_cols[0]].map(group_mapping)
            # Sort by the specified order and then by the group columns
            plot_df = plot_df.sort_values(['_group_order'] + group_cols)
        else:
            # Use alphabetical order
            plot_df = plot_df.sort_values(group_cols)
    else:
        # If no group columns, just sort by bar_label
        plot_df = plot_df.sort_values(bar_label)
    
    # Get unique facet groups and samples in order
    facet_groups = plot_df[facet_col].unique() if facet_col in plot_df else [bar_label]
    
    # Create figure with subplots for each facet
    n_facets = len(facet_groups)
    
    # Calculate figure dimensions: height controls bar area, add space for legend and labels
    bar_area_height = height
    legend_space = 0.8  # Space for legend at top
    label_space = 0.5   # Space for population labels
    total_height = bar_area_height + legend_space + label_space
    
    fig = plt.figure(figsize=(max(12, n_facets * 4), total_height))
    
    # Create main gridspec with dedicated spaces
    gs_main = gridspec.GridSpec(3, 1, height_ratios=[legend_space, bar_area_height, label_space],
                               left=0.07, right=0.95, top=0.98, bottom=0.02, hspace=0)
    
    # Legend area (top)
    ax_legend = fig.add_subplot(gs_main[0])
    ax_legend.axis('off')  # Hide the axis, we only want the legend
    
    # Bar area (middle) - create sub-gridspec for the bars
    gs_bars = gridspec.GridSpecFromSubplotSpec(1, n_facets, subplot_spec=gs_main[1],
                                              width_ratios=[len(plot_df[plot_df[facet_col] == group]) if facet_col in plot_df else len(plot_df) 
                                                          for group in facet_groups],
                                              wspace=0.1)
    
    # Label area (bottom) - for x-axis labels if needed
    ax_labels = fig.add_subplot(gs_main[2])
    ax_labels.axis('off')  # We'll use individual subplot labels
    
    # Define colors for components (using ggplot-like colors)
    if len(component_cols) == 2:
        colors = ['#e41a1c', '#377eb8']  # Red and blue for K=2
    elif len(component_cols) == 3:
        colors = ['#e41a1c', '#377eb8', '#4daf4a']  # Red, blue, green for K=3
    else:
        colors = plt.cm.Set3(np.linspace(0, 1, len(component_cols)))
    
    # Plot each facet
    axs = []  # Store axis references
    for i, group in enumerate(facet_groups):
        if facet_col in plot_df:
            group_data = plot_df[plot_df[facet_col] == group]
        else:
            group_data = plot_df
        
        ax = fig.add_subplot(gs_bars[i])
        axs.append(ax)
        
        # Prepare data for this facet - REVERSE the order of components
        components = group_data[component_cols].values
        # Reverse the components so V1 is on top
        components_reversed = components[:, ::-1]
        sample_names = group_data[bar_label].values
        
        # Create stacked bar plot with reversed components
        bottom = np.zeros(len(sample_names))
        # Use reversed colors to match reversed components
        colors_reversed = colors[::-1]
        
        for j, col in enumerate(component_cols):
            # Plot in reverse order so V1 is on top
            ax.bar(range(len(sample_names)), components_reversed[:, j], bottom=bottom, 
                   label=col, color=colors_reversed[j], width=1.0, edgecolor='white', linewidth=0.5)
            bottom += components_reversed[:, j]
        
        # Customize this subplot
        ax.set_title(str(group), fontsize=14, pad=10)
        ax.set_ylim(0, 1)
        
        # Only show y-axis for first subplot
        if i == 0:
            ax.set_ylabel('Ancestry Proportion', fontsize=14)
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
        else:
            # Remove y-axis completely for non-first subplots
            ax.set_ylabel('')
            ax.set_yticks([])  # Remove ticks completely
            # Remove y-axis spine
            ax.spines['left'].set_visible(False)
        
        # Customize x-axis - remove the thin black line and tick marks
        ax.set_xticks(range(len(sample_names)))
        ax.set_xticklabels(sample_names, rotation=90, ha='center', fontsize=xlabel_size)
        ax.set_xlabel('')
        
        # Remove the x-axis spine (thin black line)
        ax.spines['bottom'].set_visible(False)
        # Remove x-axis tick marks
        ax.tick_params(axis='x', which='both', length=0)
        
        # Add grid - only show for leftmost plot
        if i == 0:
            ax.grid(True, axis='y', alpha=0.3, linestyle='-')
        else:
            # For other plots, add grid but make it less visible or remove
            ax.grid(False, axis='y')
        
        ax.set_axisbelow(True)
        
        # Remove top and right spines (ggplot style)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
    
    # Add a common legend in the dedicated legend area (center top)
    legend_labels = component_cols
    legend_handles = [plt.Rectangle((0,0),1,1, color=colors[i]) for i in range(len(component_cols))]
    
    ax_legend.legend(legend_handles, legend_labels, loc='center', 
                     ncol=len(component_cols), frameon=False, fontsize=12)
    
    # No title as requested
    
    # Save or show plot
    if output_format != 'none':
        if not output_file:
            # Generate output file name from input file
            base_name = os.path.splitext(input_file)[0]
            output_file = f"{base_name}.{output_format}"
        else:
            # Check if output_file already has an extension
            if not any(output_file.endswith(ext) for ext in ['.png', '.svg', '.pdf']):
                output_file = f"{output_file}.{output_format}"
        
        plt.savefig(output_file, format=output_format, bbox_inches='tight')
        print(f"Plot saved as {output_file}")
    else:
        # Only show plot if not saving
        plt.show()

def main():
    parser = argparse.ArgumentParser(description='Create ADMIXTURE ancestry plots')
    parser.add_argument('--input', required=True, help='Input ADMIXTURE file')
    parser.add_argument('--group', help='Comma-separated column names for grouping (e.g., species,population1)')
    parser.add_argument('--order', help='Comma-separated order for groups (e.g., A,B,C,D)')
    parser.add_argument('--graph', choices=['none', 'svg', 'png', 'pdf'], 
                       default='png', help='Output format (default: png)')
    parser.add_argument('--output', help='Output file name (optional)')
    parser.add_argument('--group_sep', action='store_true', 
                       help='Add dotted lines separating groups (uses first group column)')
    parser.add_argument('--bar_label', default='sample', 
                       help='Column to use for bar labels (default: sample)')
    parser.add_argument('--graph_style', type=int, choices=[1, 2], default=1,
                       help='Plot style: 1=original, 2=ggplot2-style with facets (default: 1)')
    parser.add_argument('--height', type=float, default=4.0,
                       help='Bar area height for graph_style 2 (default: 4.0)')
    parser.add_argument('--Xlabel_size', type=float, default=8.0,
                       help='Font size for X-axis labels (sample names) (default: 8.0)')
    
    args = parser.parse_args()
    
    # Read input file
    df = read_admixture_file(args.input)
    print(f"Loaded {len(df)} samples from {args.input}")
    print(f"Columns found: {list(df.columns)}")
    
    # Verify bar_label column exists
    if args.bar_label not in df.columns:
        print(f"Error: Bar label column '{args.bar_label}' not found in data.")
        print(f"Available columns: {list(df.columns)}")
        return
    
    # Get component columns
    component_cols = get_component_columns(df)
    print(f"Found component columns: {component_cols}")
    
    # Parse group columns
    group_cols = []
    if args.group:
        group_cols = [col.strip() for col in args.group.split(',')]
    
    # Parse group order
    group_order = None
    if args.order:
        group_order = [item.strip() for item in args.order.split(',')]
    
    # Create plot with selected style
    if args.graph_style == 1:
        plot_admixture_1(df, component_cols, group_cols, group_order, args.graph, args.output, args.input, args.group_sep, args.bar_label, args.Xlabel_size)
    else:
        plot_admixture_2(df, component_cols, group_cols, group_order, args.graph, args.output, args.input, args.group_sep, args.bar_label, args.height, args.Xlabel_size)

if __name__ == '__main__':
    main()