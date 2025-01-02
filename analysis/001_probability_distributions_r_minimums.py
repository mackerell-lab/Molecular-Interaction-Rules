# Probability Distributions
import plotly.figure_factory as ff
import numpy as np

def normalize_data(data_list):
    normalized_data = []
    for dataset in data_list:
        # Convert to numpy array
        dataset = np.array(dataset)
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(dataset, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        normalized_data.append(dataset / area)
    return normalized_data

# Normalize the data
normalized_datasets = normalize_data([
    r_minimums_hf,
    r_minimums_wb97g,
    r_minimums_mp26
])

fig = ff.create_distplot(
    normalized_datasets,
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
    title_text='Minimum Distance (Å)',
    ticklen=15,
    range=[1, 8]
)

fig.update_yaxes(
    ticks="outside", 
    tickwidth=3,
    tickcolor='black', 
    title_text='Probability Density',
    tickfont=dict(family='Arial', color='black', size=40),
    title_font=dict(size=46, family='Arial'),
    ticklen=15,
    tickvals=[0, 0.2, 0.4, 0.6, 0.8, 1.0]  # Updated for normalized scale
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

# Probability Distributions

def normalize_data(data_list):
    normalized_data = []
    for dataset in data_list:
        # Convert to numpy array
        dataset = np.array(dataset)
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(dataset, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        normalized_data.append(dataset / area)
    return normalized_data

# Normalize the data
normalized_datasets = normalize_data([
    r_minimums_hf_charged,
    r_minimums_wb97g_charged,
    r_minimums_mp26_charged,
])

fig = ff.create_distplot(
    normalized_datasets,
    [
    'HF 6-31G* Charged',
    'WB97x-D3 6-31G* Charged',
    'MP2 6-31G* Charged'
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
                 title_text='Minimum Distance (Å)',
                 ticklen=15,
                 range=[1, 8]
)

fig.update_yaxes(
                 ticks="outside", 
                 tickwidth=3,
                 tickcolor='black', 
                 title_text='Probability Density',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 ticklen=15,
                 tickvals=[0, 0.1, 0.2, 0.3]
#                  range=[0, 0.7]
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

# Probability Distributions

import plotly.figure_factory as ff

def normalize_data(data_list):
    normalized_data = []
    for dataset in data_list:
        # Convert to numpy array
        dataset = np.array(dataset)
        # Calculate the area under the curve
        hist, bin_edges = np.histogram(dataset, bins='auto', density=True)
        bin_width = bin_edges[1] - bin_edges[0]
        area = np.sum(hist * bin_width)
        # Normalize by dividing by the area
        normalized_data.append(dataset / area)
    return normalized_data

# Normalize the data
normalized_datasets = normalize_data([
    r_minimums_hf_neutral_ion,
    r_minimums_wb97g_neutral_ion,
    r_minimums_mp26_neutral_ion,
])

fig = ff.create_distplot(
    normalized_datasets,
    [
    'HF  6-31G* Neutral Ion',
    'WB97x-D3 6-31G* Neutral Ion',
    'MP2 6-31G* Neutral Ion'
], bin_size=.5, show_hist=False, show_rug=False, show_curve=True, colors=['blue', 'hotpink', 'orange'])

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
                 title_text='Minimum Distance (Å)',
                 ticklen=15,
                 range=[1, 8]
)

fig.update_yaxes(
                 ticks="outside", 
                 tickwidth=3,
                 tickcolor='black', 
                 title_text='Probability Density',
                 tickfont=dict(family='Arial', color='black', size=40),
                 title_font=dict(size=46, family='Arial'),
                 ticklen=15,
                 tickvals=[0, 0.2, 0.4, 0.6]
#                  range=[0, 0.7]
)    

fig.update_layout(
    template='simple_white',
    xaxis_tickformat = 'i',
    bargap=0.2,
    height=600,
    width=1000,
    showlegend=False
)    
fig.show()
