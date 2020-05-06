"""
Data visualization functions
"""
print('Plotting function moves to ALLCools, this is deprecated, will be deleted.')

# preprocessing related plots
from .preprocessing import \
    plot_on_plate, \
    cutoff_vs_cell_remain, \
    simple_violin, \
    success_vs_fail, \
    plot_dispersion

# color related
from .color import \
    level_one_palette, \
    level_two_palette, \
    get_kv_dict, \
    palplot

# cell scatter deprecated_plot and related functions
from .cell import \
    density_based_sample, \
    categorical_scatter, \
    continuous_scatter

# cluster pie and sunbrust
from .pie import \
    sunbrust

# cluster tree
from .tree import \
    generate_tree_dict, \
    dendrogram

from .group import \
    stack_bar_plot

from .sankey import \
    Sankey
