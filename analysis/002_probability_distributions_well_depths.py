# Probability Distributions

import plotly.figure_factory as ff

fig = ff.create_distplot([
    well_depth_minimums_hf_charged,
    well_depth_minimums_wb97g_charged,
    well_depth_minimums_mp26_charged,
], [
    'HF  6-31G* Charged',
    'WB97x-D3 6-31G* Charged',
    'MP2 6-31G* Charged'
], bin_size=.5, show_hist=False, show_rug=False, show_curve=True, 
colors=['blue', 'hotpink', 'orange'])

fig.update_traces(line=dict(width=10))

fig.update_layout(legend=dict(itemsizing='constant'))
fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,  
        font = dict(family = "Arial", size = 10),
        bordercolor="LightSteelBlue",
        borderwidth=2,
    ),
    legend_title = dict(font = dict(family = "Arial", size = 20))
)

fig.update_xaxes(
                 ticks="outside",
                 tickwidth=3,
                 tickcolor='black',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 title_text='Minimum Energy (Kcal/mol)',
                 ticklen=15,
                 range=[-140, 140],
)

fig.update_yaxes(
                 ticks="outside", 
                 tickwidth=3,
                 tickcolor='black', 
                 title_text='Probability Density',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 ticklen=15,
                range=[0, 0.008],
                tickvals=[0, 0.002, 0.004, 0.006, 0.008]
)    

fig.update_layout(
    template='simple_white',
    xaxis_tickformat = 'i',
    bargap=0.2,
    height=600,
    width=1000,
    showlegend=False,
)    
fig.show()

fig.write_image('well_depths_probability_distributions_charged.png')

import plotly.figure_factory as ff
import numpy as np

# First filter and normalize each dataset
def filter_and_normalize_data(data_list):
    filtered_data = []
    for dataset in data_list:
        # Convert to numpy array if it isn't already
        dataset = np.array(dataset)
      
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(dataset, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        filtered_data.append(filtered / area)
    return filtered_data

# Filter and normalize the data
filtered_normalized_datasets = filter_and_normalize_data([
    well_depth_minimums_hf_neutral_ion,
    well_depth_minimums_wb97g_neutral_ion,
    well_depth_minimums_mp26_neutral_ion,
])

# Create the distribution plot with filtered normalized data
fig = ff.create_distplot(
    filtered_normalized_datasets,
    [
    'HF  6-31G* Neutral Ion',
    'WB97x-D3 6-31G* Neutral Ion',
    'MP2 6-31G* Neutral Ion'
    ],
    bin_size=.5,
    show_hist=False,
    show_rug=False,
    show_curve=True,
    colors=['blue', 'hotpink', 'orange']
)

fig.update_traces(line=dict(width=5))
fig.update_layout(legend=dict(itemsizing='constant'))
fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,  
        font = dict(family = "Arial", size = 10),
        bordercolor="LightSteelBlue",
        borderwidth=2,
    ),
    legend_title = dict(font = dict(family = "Arial", size = 20))
)

fig.update_xaxes(
    ticks="outside",
    tickwidth=3,
    tickcolor='black',
    tickfont=dict(family='Arial', color='black', size=40),
    title_font=dict(size=46, family='Arial'),
    title_text='Minimum Energy (Kcal/mol)',
    ticklen=15,
    range=[-40, 40],
)

fig.update_yaxes(
    ticks="outside", 
    tickwidth=3,
    tickcolor='black', 
    title_text='Probability Density',
    tickfont=dict(family='Arial', color='black', size=40),
    title_font=dict(size=46, family='Arial'),
    ticklen=15,
    tickvals=[0, 0.02, 0.04, 0.06, 0.08]
)    

fig.update_layout(
    template='simple_white',
    xaxis_tickformat = 'i',
    bargap=0.2,
    height=600,
    width=1000,
    showlegend=False,
)    

fig.show()
fig.write_image('well_depths_probability_distributions_charged.png')

# Probability Distributions

import plotly.figure_factory as ff

# First normalize each dataset
def normalize_data(data_list):
    normalized_data = []
    for dataset in data_list:
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(dataset, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        normalized_data.append(dataset / area)
    return normalized_data
    
well_depth_minimums_hf_neutral_ion = [i for i in well_depth_minimums_hf_neutral_ion if i < 200]
well_depth_minimums_wb97g_neutral_ion = [i for i in well_depth_minimums_wb97g_neutral_ion if i < 200]
well_depth_minimums_mp26_neutral_ion = [i for i in well_depth_minimums_mp26_neutral_ion if i < 200]

# Normalize the data
normalized_datasets = normalize_data([
    well_depth_minimums_hf_neutral_ion,
    well_depth_minimums_wb97g_neutral_ion,
    well_depth_minimums_mp26_neutral_ion,
])

fig = ff.create_distplot(
    normalized_datasets, [
    'HF  6-31G* Neutral Ion',
    'WB97x-D3 6-31G* Neutral Ion',
    'MP2 6-31G* Neutral Ion'
], bin_size=.5, show_hist=False, show_rug=False, show_curve=True, 
colors=['blue', 'hotpink', 'orange'])

fig.update_traces(line=dict(width=10))

fig.update_layout(legend=dict(itemsizing='constant'))
fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,  
        font = dict(family = "Arial", size = 10),
        bordercolor="LightSteelBlue",
        borderwidth=2,
    ),
    legend_title = dict(font = dict(family = "Arial", size = 20))
)

fig.update_xaxes(
                 ticks="outside",
                 tickwidth=3,
                 tickcolor='black',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 title_text='Minimum Energy (Kcal/mol)',
                 ticklen=15,
#                  range=[-50, 50],
)

fig.update_yaxes(
                 ticks="outside", 
                 tickwidth=3,
                 tickcolor='black', 
                 title_text='Probability Density',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 ticklen=15,
#                  tickvals=[0, 0.025, 0.05, 0.075, 0.1, 0.125]
)    

fig.update_layout(
    template='simple_white',
    xaxis_tickformat = 'i',
    bargap=0.2,
    height=600,
    width=1000,
    showlegend=False,
    
)    
fig.show()

fig.write_image('well_depths_probability_distributions_neutral_ion.png')

import plotly.figure_factory as ff
import numpy as np

# First filter and normalize each dataset
def filter_and_normalize_data(data_list):
    filtered_data = []
    for dataset in data_list:
        # Convert to numpy array if it isn't already
        dataset = np.array(dataset)
        # Filter out values outside the range
        filtered = dataset[(dataset >= -140) & (dataset <= 140)]
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(filtered, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        filtered_data.append(filtered / area)
    return filtered_data

# Filter and normalize the data
filtered_normalized_datasets = filter_and_normalize_data([
    well_depth_minimums_hf,
    well_depth_minimums_wb97g,
    well_depth_minimums_mp26,
])

# Create the distribution plot with filtered normalized data
fig = ff.create_distplot(
    filtered_normalized_datasets,
    [
    'HF 6-31G*',
    'WB97x-D3 6-31G*',
    'MP2 6-31G*'
    ],
    bin_size=.5,
    show_hist=False,
    show_rug=False,
    show_curve=True,
    colors=['blue', 'hotpink', 'orange']
)

fig.update_traces(line=dict(width=5))
fig.update_layout(legend=dict(itemsizing='constant'))
fig.update_layout(legend=dict(
        orientation="h",
        yanchor="bottom",
        y=1.02,
        xanchor="right",
        x=1,  
        font = dict(family = "Arial", size = 10),
        bordercolor="LightSteelBlue",
        borderwidth=2,
    ),
    legend_title = dict(font = dict(family = "Arial", size = 20))
)

fig.update_xaxes(
    ticks="outside",
    tickwidth=3,
    tickcolor='black',
    tickfont=dict(family='Arial', color='black', size=40),
    title_font=dict(size=46, family='Arial'),
    title_text='Minimum Energy (Kcal/mol)',
    ticklen=15,
    range=[-20, 20],
)

fig.update_yaxes(
    ticks="outside", 
    tickwidth=3,
    tickcolor='black', 
    title_text='Probability Density',
    tickfont=dict(family='Arial', color='black', size=40),
    title_font=dict(size=46, family='Arial'),
    ticklen=15,
    tickvals=[0, 0.04, 0.08, 0.12]
)    

fig.update_layout(
    template='simple_white',
    xaxis_tickformat = 'i',
    bargap=0.2,
    height=600,
    width=1000,
    showlegend=False,
)    

fig.show()
fig.write_image('well_depths_probability_distributions.png')
