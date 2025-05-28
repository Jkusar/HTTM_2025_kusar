# -*- coding: utf-8 -*-
"""
Created on Wed Aug  9 13:36:57 2023

@author: Jernej Kušar
"""

import sys
sys.path.append("C:/Users/Jernej Kusar/Documents/LFDT splošno/Dodiplomska/Main")
sys.path.append("C:/Users/Jernej Kusar/Documents/LFDT splošno/Dodiplomska/Enacbe")


import numpy as np
import datalib as jdl
import matplotlib.pyplot as plt
from collections import defaultdict
import control_values as ctrl
import program_values as pval
import seaborn as sns
import random
from matplotlib.ticker import MultipleLocator, ScalarFormatter
import matplotlib.markers as markers
from matplotlib.colors import Normalize
from matplotlib.patches import Patch
import copy

#progress bar
from tqdm import tqdm



def random_pleasant_color():
    """
    Returns a random pleasant color. Meant to be used for graphs
    
    Output parameters:
        - random_color [str]
    """
    # Choose a color palette name from seaborn
    palette_name = random.choice(list(sns.palettes.SEABORN_PALETTES.keys()))
    # Get the palette colors
    palette = sns.palettes.SEABORN_PALETTES[palette_name]
    # Choose a random color from the palette
    random_color = random.choice(palette)
    
    return random_color


def random_markers(size=1):
    """
    Generate random symbols from Matplotlib supported markers suitable for graphs.

    Input parameters:
        - size [int] - number of generated symbols

    Output parameters:
        - random_markers [list or str] - list of random markers if size > 1, otherwise a single random marker as a string
    """
    # List of markers suitable for graphs (excluding lines and empty markers)
    suitable_markers = [
        'o', 's', '^', 'D', 'v', 'p', 'h', '8', '<', '>', 
        'H', 'X', 'P'
    ]
    
    random_markers = np.random.choice(suitable_markers, size=size)
    
    if size == 1:
        return random_markers[0]
    else:
        return random_markers.tolist()

def get_random_marker():
    """
    Returns a random marker from matplotlibs available marker styles.
    
    Output parameters:
        - random_marker [str] - symbol that matplotlib functions recognise
    """
    # Get a list of available marker symbols
    available_markers = markers.MarkerStyle.markers.keys()
    # Convert the keys to a list and shuffle it
    marker_list = list(available_markers)
    random.shuffle(marker_list)
    # Selects a random marker and returns it
    return random.choice(marker_list)

def get_marker_easy(symbol):
    """
    Returns a marker that matplotlib functions will recognise.
    
    Input parameters:
        - symbol [str] - name of the marker that you wish to get, eg. 'point' returns '.'
    
    Output parameters:
        -marker [str] - symbol for a marker that matplotlib recognises
    """
    marker = markers.MarkerStyle(symbol)
    return marker.get_path().vertices

def relative_differences_plot(base_names, nozzle_parameters, values_a, values_b, x_label, y_label, title, **kwargs):
    """
    Plots graphs showing relative differences
    
    Input parameters:
        - base_names [list] - names of nozzles
        - nozzle_parameters [dict] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
        - values_a [dict] - calculated values. Dictionary is writen; {base_name : value}
        - values_b [dict] - meassured values. Dictionary is writen; {base_name : value}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - title [str] - Title of the graph
        - kwargs:
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values eg. RPx.1, RPx.2,...
            - iterations_split [str, optional] - character that differs iterations of data eg. RPx_1, RPx_2,... 
            - symbol_palett [dict, optional] - dictionary from which symbols are defined for nozzles
            - color_palett [dict, optional] - dictionary in which colors are defined to certain nozzles
    """
    nozzles_split = kwargs.get("nozzles_split", ".")
    iterations_split = kwargs.get("iterations_split", "_")
    symbol_palett = kwargs.get("symbol_palett", None)
    color_palett = kwargs.get("color_palett", None)
    both = jdl.compare_lists(list1=list(values_a), list2=list(values_b))
    difference = {}
    for nozzle_name in base_names:
        both_in = jdl.filter_names(both, include=[f"{nozzle_name}-"])
        for i in both_in:
            diff = (values_a.get(i) - values_b.get(i))
            difference.update({i : diff})
    #Plots graphs showing relative differences in calculated and meassured jet diameter for different nozzles
    values = {}
    for nozzle_name in base_names:
        local_values = {}
        names_int = np.array(jdl.filter_names(list(difference), include = [f"{nozzle_name}-"]))
        names_int = sorted(names_int, key=lambda item: (int(item.split('-')[1][1:]), int(item.split('-')[2][1:])))
        ratio = nozzle_parameters.get(nozzle_name.split(nozzles_split)[0].split(iterations_split)[0])[-1]
        for name in names_int:
            value = difference.get(name)/values_a.get(name)*100
            point = name.strip(f"{nozzle_name}-")
            local_values.update({point: value})
        values.update({f"{nozzle_name} ; {ratio}": local_values})
    # Create a plot
    plt.figure(figsize=(8, 6))
    # Use defaultdict to store unique x-values and their corresponding y-values
    x_y_values = defaultdict(list)
    labels = set()  # Use a set to store unique labels
    for label, data in values.items():
        for x, y in data.items():
            x_y_values[x].append((y, label))  # Store the y-value and label together
            labels.add(label)  # Add the label to the set
    if symbol_palett is None:
        symbol_palett = {}
        for label in labels:
            symbol_palett.update({label : random_markers()})
    if color_palett is None:
        color_palett = {}
        for label in labels:
            color_palett.update({label: random_pleasant_color()})
    # Sort x_y_values by its keys (x-values) based on G and L integrals
    sorted_x_y_values = sorted(x_y_values.items(), key=lambda item: (int(item[0][1: item[0].index('-')]), int(item[0][item[0].index('L') + 1:])))
    # Create a dictionary to map each unique label to a specific color and symbol
    color_mapping = {label: color_palett.get(label.split(" ; ")[0]) for label in labels}
    symbol_mapping = {label: symbol_palett.get(label.split(" ; ")[0]) for label in labels}
    # Plot all points with the same x-value using the same color and add to legend only once per label
    legend_labels = []
    for x_value, y_label_pairs in sorted_x_y_values:
        y_values, label = zip(*y_label_pairs)  # Unzip the y-value and label pairs
        for y_value, label in zip(y_values, label):
            if label not in legend_labels:
                plt.scatter(x_value, y_value, label=label, color=color_mapping.get(label), marker = symbol_mapping.get(label), s=100)
                legend_labels.append(label)
            else:
                plt.scatter(x_value, y_value, color=color_mapping.get(label), marker = symbol_mapping.get(label), s=100)
    # Customize the plot
    plt.xlabel(x_label, fontsize=30)
    plt.ylabel(y_label, fontsize=30)
    plt.title(title, fontsize=50)
    handles, labels = plt.gca().get_legend_handles_labels()
    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: t[0]))
    plt.legend(handles, labels, fontsize=20)
    plt.grid()
    plt.xticks(rotation=90, fontsize=15)
    plt.yticks(fontsize=20)
    # Show the plot
    plt.show()

def plot_two_sets_of_data(nozzles, values_a, values_b, x_label, y_label, **kwargs):
    """
    Plots a plot consisting of two sets of data (eg. pressure data) for specified nozzles
    
    Input parameters:
        - nozzles [list] - names of nozzles
        - values_a [dict] - calculated values. Dictionary is writen; {base_name : value}
        - values_b [dict] - meassured values. Dictionary is writen; {base_name : value}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - kwargs:
            - plot_vertical_lines [bool, optional] - option to plot vertical lines in plot
            - v_line_color [list, optional] - colors of vertical lines
            - v_line_style [list, optional] - styles of vertical lines
            - v_line_coords [list, optional] - list of x coordinates for ploting vertical lines
            - nozzle_parameters [dict, optional] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
            - x_values [dict, optional] - values on x axis. By default x-axis values are keys of values a and b.
            - x_axis_spacing [int, optional] - Spacing between x ticks and axis label
            - y_axis_spacing [int, optional] - Spacing between y ticks and axis label
            - label_a [str, optional] - Label of values a (default it is values a)
            - label_b [str, optional] - Label of values b (default it is values b)
            - color_a [str, optional] - color of a points. Default it is red
            - color_b [str, optional] - color of b points. Default it is green
            - figure_size [tuple, optional] - size of the created figure (x, y)
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values eg. RPx.1, RPx.2,...
            - symbol_a [str, optional] - symbol for a values
            - symbol_b [str, optional] - symbol for b values
            - save_plot [str, optional] - option to save plot in a specified directory
            - size_a [int, optional] - size for a values
            - size_b [int, optional] - size for b values
            - show_nozzle [bool, optional] - option to show/hide nozzle name in title
            - iterations_split [str, optional] - character that differs iterations of data eg. RPx_1, RPx_2,... 
            - tick_size [int, optional] - size of tick characters
            - title [str, optional] - Title of the graph
            - label_size [int, optional] - size of label characters
            - legend_size [int, optional] - size of legend characters
    """
    title = kwargs.get("title", "")
    figure_size = kwargs.get("figure_size", (5, 5))
    show_nozzle = kwargs.get("show_nozzle", False)
    x_names = kwargs.get("x_values", None)
    label_a = kwargs.get("label_a", "values_a")
    label_b = kwargs.get("label_b", "values_b")
    color_a = kwargs.get("color_a", "red")
    color_b = kwargs.get("color_b", "green")
    symbol_a = kwargs.get("symbol_a", ".")
    symbol_b = kwargs.get("symbol_b", "s")
    save_plot = kwargs.get("save_plot", None)
    size_a = kwargs.get("size_a", 50)
    size_b = kwargs.get("size_b", 50)
    sort_x = kwargs.get("sort_x_values", True)
    sort_key = kwargs.get("sort_key", pval.nozzle_sorting_key) #ključ je po defaultu za sortiranje imen nozzlov
    nozzles_split = kwargs.get("nozzles_split", ".")
    iterations_split = kwargs.get("iterations_split", "_")
    both = jdl.compare_lists(list1=list(values_a), list2=list(values_b))
    nozzle_param = kwargs.get("nozzle_parameters", None)
    plot_vline = kwargs.get("plot_vertical_lines", False)
    plot_seperate = kwargs.get("plot_seperate", False)
    v_line_coords = kwargs.get("v_line_coords", None)
    v_line_color = kwargs.get("v_line_color", None)
    v_line_style = kwargs.get("v_line_style", "--")
    tick_size = kwargs.get("tick_size", 5)
    label_size = kwargs.get("label_size", 5)
    legend_size = kwargs.get("legend_size", 5)
    legend_location = kwargs.get("legend_location", "upper left")
    y_tick_axis_spacing = kwargs.get("y_axis_spacing", 5)
    x_tick_axis_spacing = kwargs.get("x_axis_spacing", 5)

    for nozzle in nozzles:
        valuesA_filtr = {}
        valuesB_filtr = {}
        if nozzle_param is not None:
            ratio = nozzle_param.get(nozzle.split(nozzles_split)[0].split(iterations_split)[0])[-1]
        both_in = jdl.filter_names(both, include=[f"{nozzle}-"])
        if sort_x:
            both_in = sorted(both_in, key=sort_key)
        for i in both_in:
            va_filtr = values_a.get(i)
            vb_filtr = values_b.get(i)
            if x_names is not None:
                try:
                    name = x_names.get(i)
                except KeyError:
                    print(f"Element {i} does not have a specified x value. Please corect it")
            else:
                name = i.strip(f"{nozzle}-")
            valuesA_filtr.update({name: va_filtr})
            valuesB_filtr.update({name: vb_filtr})
            
        names = list(valuesA_filtr)
        bp_pos = []
        for i in range(len(names)):
            bp_pos.append(i)
        pointsA = list(valuesA_filtr.values())
        pointsB = list(valuesB_filtr.values())         
        
        fig = plt.figure()
        ax = fig.add_subplot()
        fig.set_figwidth(figure_size[0])
        fig.set_figheight(figure_size[1])
        
        new_names = []
        for name in names:
            name_split = name.split("-")[0:]
            try:
                L = int(name_split[-1].strip("L"))*10**-2
                G = round(int(name_split[0].strip("G")), 2)
            except ValueError:
                continue
            L_cor = str(L)
            new_point= "G" + str(G) + "-" + "L" + L_cor[:2]
            if len(new_point) < 10:
                new_G = str(G)
                new_L = L_cor[:2]
                if len(str(G)) < 2:
                    new_G = "0" + str(G)
                if "." in new_L:
                    new_L = "0" + L_cor[:1]
                
                new_point = "G" + new_G + "-" + "L" + new_L
                new_names.append(new_point)
        try:
            pointsA[0]
            pointsB[0]
        except IndexError:
            continue
        
        if type(pointsA[0]) is not list:
            sc1 = ax.scatter(names, pointsA, c=color_a, label=label_a, marker=symbol_a, s=size_a)
            handles = [sc1]
        else:
            flierprops = dict(marker=symbol_a, markersize=size_a)
            bp1 = plt.boxplot(pointsA, positions = bp_pos, notch=True, flierprops=flierprops, showfliers=False, 
                        showbox=False, medianprops = dict(markerfacecolor = color_a, 
                                                           markeredgecolor = color_a,
                                                           markersize=size_a), 
                        showmeans=True, meanprops = dict(markerfacecolor = color_a,
                                                              markeredgecolor = color_a, 
                                                              marker = symbol_a, markersize=size_a),
                        widths=1.7, patch_artist=True)
            handles = [bp1["means"][0]]
            
        if type(pointsB[0]) is not list:
            sc2 = ax.scatter(names, pointsB, c=color_b, label=label_b, marker=symbol_b, s=size_b)
            handles.append(sc2)
        else:
            flierprops = dict(marker=symbol_b, markersize=size_b)
            bp2 = plt.boxplot(pointsB,  notch=True, flierprops=flierprops, showfliers=False, 
                        showbox=False, medianprops = dict(markerfacecolor = color_b, 
                                                           markeredgecolor = color_b,
                                                           markersize=size_b), 
                        showmeans=False, meanprops = dict(markerfacecolor = color_b,
                                                              markeredgecolor = color_b, 
                                                              marker = symbol_b, markersize=size_b),
                        widths=1.7)
            handles.append(bp2["means"][0])
        #plot vertical lines
        if plot_vline:
            for line in range(len(v_line_coords)):
                ax.vlines(x=v_line_coords[line], colors=v_line_color[line], linestyles=v_line_style[line])
        
        ax.set_xlabel(x_label, fontsize=label_size)
        ax.set_ylabel(y_label, fontsize=label_size)
        if nozzle_param is not None and show_nozzle: 
            ax.set_title(label=f"{title} {nozzle} ; {ratio}")
        if show_nozzle:
            ax.set_title(label=f"{title} {nozzle}")
        else:
            ax.set_title(label=f"{title}")
            
        ax.set_xticks(names)
        ax.set_xticklabels(new_names, rotation=90)
        
        plt.xticks(fontsize=tick_size)
        plt.yticks(fontsize=tick_size)
        plt.tick_params(axis='y', pad=y_tick_axis_spacing)
        plt.tick_params(axis='x', pad=x_tick_axis_spacing)
        ax.legend(handles, [label_a, label_b], fontsize=legend_size, loc=legend_location)
        ax.grid()
        
        #plt.title(f"{nozzle}", fontsize=80)
        
        if save_plot is not None:
            plt.savefig(save_plot, bbox_inches = 'tight')
        else:
            continue
        plt.show()


def plot_multiple_sets_of_data(nozzles, data_list, x_label, y_label, title, **kwargs):
    """
    Plots multiple sets of data for specified nozzles.
    
    Input parameters:
        - nozzles [list] - names of nozzles
        - data_list [list] - list of dictionaries containing data for each nozzle.
                             Each dictionary has nozzle names as keys and corresponding data as values.
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - title [str] - Title of the graph
        - kwargs: Optional keyword arguments:
            - plot_vertical_lines [bool] - Option to plot vertical lines in the plot
            - v_line_coords [list] - List of x coordinates for plotting vertical lines
            - v_line_color [list] - Colors of vertical lines
            - v_line_style [list] - Styles of vertical lines
            - nozzle_parameters [dict] - Geometrical parameters of nozzles. Dictionary is written as {base_name:nozzle_parameter}
            - labels [list] - List of suffixes for legend labels
            - colors [list] - List of colors for plotting data points
            - markers [list] - List of markers for plotting data points
            - sizes [list] - List of sizes for plotting data points
    """
    
    plot_vline = kwargs.get("plot_vertical_lines", False)
    v_line_coords = kwargs.get("v_line_coords", None)
    v_line_color = kwargs.get("v_line_color", None)
    v_line_style = kwargs.get("v_line_style", "--")
    labels = kwargs.get("labels", None)
    colors = kwargs.get("colors", ["red", "green", "blue", "orange", "purple"])
    markers = kwargs.get("markers", ["o", "s", "v", "^", "x"])
    sizes = kwargs.get("sizes", [50, 50, 50, 50, 50])
        
    if labels is None:
        labels=[]
        a = 1
        for i in range(len(data_list)):
            labels.append(f"Points {a}")
            a+=1
    if colors is None:
        colors=[]
        a = 1
        for i in range(len(data_list)):
            colors.append(random_pleasant_color())
            a+=1
    if markers is None:
        markers=[]
        a = 1
        for i in range(len(data_list)):
            markers.append(random_markers())
            a+=1
    if sizes is None:
        sizes=[]
        a = 1
        for i in range(len(data_list)):
            sizes.append(50)
            a+=1
    
    nozzle_all_values = {}
    for i in nozzles:
        filtered_list = []
        for j in data_list:
             filtered_list.append(jdl.filter_dict_keys_by_substring(input_dict=j, substring=i))
        nozzle_all_values.update({i : filtered_list})

    
    for nozzle in nozzles:
        fig, ax = plt.subplots()
        nozzle_val = nozzle_all_values.get(nozzle)
        for i in range(len(nozzle_val)):
            x_values = list(nozzle_val[i].keys())
            for j in range(len(x_values)):
                x_values[j] = x_values[j].strip(f"{nozzle}-")
            y_values = list(nozzle_val[i].values())
            color = colors[i]
            marker = markers[i]
            size = sizes[i]
            label = labels[i]
            
            ax.scatter(x_values, y_values, c=color, label=f'{nozzle} {label}', marker=marker, s=size)
        
        # Plot vertical lines
        if plot_vline:
            for line in range(len(v_line_coords)):
                ax.vlines(x=v_line_coords[line], colors=v_line_color[line], linestyles=v_line_style[line])
        
        ax.set_xlabel(x_label)
        ax.set_ylabel(y_label)
        ax.set_title(title)
        ax.set_xticks(x_values)
        ax.set_xticklabels(x_values, rotation=90)
        ax.legend()
        ax.grid()
        plt.show()
        
    
def boxplot_one_dataset(nozzles, values, x_label, y_label, title, **kwargs):
    """
    Plots a boxplot of sets of data provided.
    
    Input parameters:
        - nozzles [str] - names of nozzles
        - values [dict] - ploted values. Dictionary is writen; {base_name : value}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - title [str] - Title of the graph
        - kwargs:
            - axis_label_pad [int, optional] - option to increase or decrease the space between axis labels and ticks
            - nozzle_parameters [dict] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
            - filter_values [list, optional] - List of base_names that ploted values must contain
            - figure_size [tuple, optiona]  - option to resize a figure (x, y)
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values eg. RPx.1, RPx.2,...
            - iterations_split [str, optional] - character that differs iterations of data eg. RPx_1, RPx_2,... 
            - show_fliers [bool, optional] - option to hide outliers from the plot
            - show_box [bool, optional] - option to hide box of boxplot
            - x_font_style [str, optional] - font style of x axis label
            - y_font_style [str, optional] - font style of y axis label
            - show_title [bool, optional] - option to hide the title of the plot
            - show_mean [bool, optional] - option to show the meanline
            - marker_size [int, optional] - size of the mean line marker
            - mean_marker_style [str, optional] - option to change the style of the mean marker
            - meanline_color [str, optional] - option to change the color of the meanline
            - multiply_values [float, optional] - option to multipy y values
            - show_median [bool, optional] - option to show the median_line
            - median_line_color [str, optional] - option to change the color of the medianline
            - tick_size [int, optional] - size of x and y ticks
            - legend_size [int, optional] - size of text in the legend
            - axis_label_size [int, optional] - size of text in x and y axis labels
    """
    filter_values = kwargs.get("filter_values", None)
    figure_size = kwargs.get("figure_size", (5,5))
    nozzles_split = kwargs.get("nozzles_split", ".")
    iterations_split = kwargs.get("iterations_split", "_")
    nozzle_param = kwargs.get("nozzle_parameters", None)
    sh_fl = kwargs.get("show_fliers", True)
    sh_box = kwargs.get("show_box", True)
    x_font_style = kwargs.get("x_font_style", "italic")
    y_font_style = kwargs.get("y_font_style", "italic")
    show_title = kwargs.get("show_title", True)
    multiply_values = kwargs.get("multiply_values", 1)
    #mean line
    show_mean = kwargs.get("show_mean", False)
    mean_marker = kwargs.get("mean_marker_style", "s")
    mean_color = kwargs.get("meanline_color", "green")
    marker_size = kwargs.get("marker_size", 30)
    #median line
    show_median = kwargs.get("show_median", True)
    median_color = kwargs.get("median_line_color", "orange")
    #font
    tick_size = kwargs.get("tick_size", 5)
    legend_size = kwargs.get("legend_size", 5)
    axis_label_size = kwargs.get("axis_label_size", 5)
    label_pad = kwargs.get("label_pad", 5)
    
    if filter_values is not None: 
        both = jdl.compare_lists(list1=list(values), list2=filter_values)
    else:
        both = list(values)
    if nozzle_param is None:
        no_ratio = True
    else:
        no_ratio = False
    
    #median line
    if show_median is False:
        median_line = ""
    else:
        median_line = "solid"
            
    plt.figure(figsize=figure_size)
    for nozzle in nozzles:
        flierprops = dict(marker='o', markersize=2)
        names = np.array(jdl.filter_names(both, include = [f"{nozzle}-"]))
        names = sorted(names, key=lambda item: (int(item.split('-')[1][1:]), int(item.split('-')[2][1:])))
        if not no_ratio:
            ratio = nozzle_param.get(nozzle.split(nozzles_split)[0].split(iterations_split)[0])[-1]
            if type(ratio) is not float:
                ratio = nozzle_param.get(nozzle.split(".")[0])[-1]
        points = []
        plot_values = []
        for name in names: 
            point = name.strip(f"{nozzle}-")
            name_split = point.split("-")
            L = int(name_split[-1].strip("L"))*10**-2
            G = round(int(name_split[0].strip("G")), 2)
            L_cor = str(L)
            new_point= "G" + str(G) + "-" + "L" + L_cor[:2]
            if len(new_point) < 10:
                new_G = str(G)
                new_L = L_cor[:2]
                if len(str(G)) < 2:
                    new_G = "0" + str(G)
                if "." in new_L:
                    new_L = "0" + L_cor[:1]
                
                new_point = "G" + new_G + "-" + "L" + new_L

            points.append(new_point)
            plot_values.append(np.array(values[name])*multiply_values)

        if type(mean_marker) == str:
            mean_marker_dict = {}
            for i in nozzles:
                mean_marker_dict.update({i : mean_marker})
            mean_marker = mean_marker_dict

        plt.boxplot(plot_values[0:],  notch=True, flierprops=flierprops, showfliers=sh_fl, 
                    showbox=sh_box, medianprops = dict(markerfacecolor = median_color, 
                                                       markeredgecolor = median_color,
                                                       ls=median_line, markersize=marker_size), 
                    showmeans=show_mean, meanprops = dict(markerfacecolor = mean_color,
                                                          markeredgecolor = mean_color, 
                                                          marker = mean_marker[nozzle], markersize=marker_size),
                    widths=1.7)
    
        plt.xlabel(x_label, fontstyle=x_font_style, fontsize=axis_label_size, labelpad=label_pad)
        plt.ylabel(y_label, fontstyle=y_font_style, fontsize=axis_label_size, labelpad=label_pad)
        if show_title:
            if no_ratio:
                plt.title(label=f"{title} {nozzle}")
            else:
                plt.title(label=f"{title} {nozzle} ; {ratio}")
        x_tick_points = [[], points]
        for i in range(len(points)):
            x_tick_points[0].append(i+1)
        plt.xticks(x_tick_points[0], x_tick_points[1], rotation=90, fontsize=tick_size)
        plt.yticks(fontsize=tick_size)
        #plt.legend(fontsize=legend_size)
        #plt.title(f"{ctrl.nozzle_type}{nozzle}")
        plt.grid()
    plt.show()
        
def plot_3D_graphs(nozzles, values_a, values_b, nozzle_parameters, x_label, y_label, z_label, title, **kwargs):
    """
    Plots a 3D graph of values a and b
    
    Input parameters:
        - nozzles [list] - names of nozzles
        - values_a [dict] - calculated values. Dictionary is writen; {base_name : value}
        - values_b [dict] - meassured values. Dictionary is writen; {base_name : value}
        - nozzle_parameters [dict] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - z_label [str] - Label of z axis
        - title [str] - Title of the graph
        - kwargs:
            - label_a [str, optional] - Label of values a (default it is values a)
            - label_b [str, optional] - Label of values b (default it is values b)
            - color_a [str, optional] - color of a points. Default it is red
            - color_b [str, optional] - color of b points. Default it is green
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values eg. RPx.1, RPx.2,...
            - iterations_split [str, optional] - character that differs iterations of data eg. RPx_1, RPx_2,... 
            - marker_size [int, optional] - Size of the markers
            - x_font_style [str, optional] - font style of x axis label
            - y_font_style [str, optional] - font style of y axis label
            - z_font_style [str, optional] - font style of z axis label
            - axis_fontsize [int, optional] - font size of x, y and z axis labels
            - title_fontsize [int, optional] - font size of title
            - title_font_style [str, optional] - font style of title
            - tick_size [int, optional] - font size of x, y and z ticks
    """
    label_a = kwargs.get("label_a", "values_a")
    label_b = kwargs.get("label_b", "values_b")
    color_a = kwargs.get("color_a", "red")
    color_b = kwargs.get("color_b", "green")
    nozzles_split = kwargs.get("nozzles_split", ".")
    iterations_split = kwargs.get("iterations_split", "_")
    marker_size = kwargs.get("marker_size", 5)
    axis_size = kwargs.get("axis_fontsize", 10)
    x_fontstyle = kwargs.get("x_font_style", "italic")
    y_fontstyle = kwargs.get("y_font_style", "italic")
    z_fontstyle = kwargs.get("z_font_style", "italic")
    title_size = kwargs.get("title_fontsize", 20)
    title_style = kwargs.get("title_font_style", "italic")
    tick_size = kwargs.get("tick_size", 3)
    #compare lists of names for jet dictionaries
    both = jdl.compare_lists(list1=list(values_a), list2=list(values_b))
    #plots 3d graph of meassured and calculated jet diameters
    difference = {}
    for nozzle in nozzles:
        d_m = []
        d_c = []
        G = []
        L = []
        ratio = nozzle_parameters.get(nozzle.split(nozzles_split)[0].split(iterations_split)[0])[-1]
        both_in = jdl.filter_names(both, include=[f"{nozzle}-"])
        for i in both_in:
            G.append(int(i.split("-")[1].strip("G")))
            L.append(int(i.split("-")[-1].strip("L")))
            d_m.append(values_b.get(i))
            d_c.append(values_a.get(i))
            diff = (values_a.get(i) - values_b.get(i))
            difference.update({i : diff})
        d_m = np.array(d_m)
        d_c = np.array(d_c)
        G = np.array(G)
        L = np.array(L)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        ax.scatter(L, G, d_m, color_b, label = label_b, s = marker_size)
        ax.scatter(L, G, d_c, color_a, label = label_a, s = marker_size)
        ax.autoscale(enable=True)
        ax.set_xlabel(x_label, fontsize=axis_size, fontstyle=x_fontstyle)
        ax.set_ylabel(y_label, fontsize=axis_size, fontstyle=y_fontstyle)
        ax.set_zlabel(z_label, fontsize=axis_size, fontstyle=z_fontstyle)
        ax.set_title(f"{title} {nozzle} ; {ratio}", fontsize=title_size, fontstyle=title_style)
        ax.legend(fontsize=tick_size)
        ax.grid()
        plt.show()
        
def plot_xy_data(nozzles, values_x, values_y, x_label, y_label, title, **kwargs):
    """
    Plots two sets of data as x and y coordinates. 
    
    Input parameters:
        - nozzles [list] - List of nozzle names
        - values_x [dict] - x values on graph. Dictionary is writen; {base_name : value}
        - values_y [dict] - y values on graph. Dictionary is writen; {base_name : value}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - title [str] - Title of the graph
        - kwargs:
            - axis_fontsize [int, optional] - font size of x and y axis labels
            - axis_scalar_format [bool, optional] - option to switch between scientific and normal tick display
            - colorbar [dict, option] - dict of colorbar values. Keys must mach those of x and y axis
            - colorbar_label [str, optional] - label of the colorbar
            - colorbar_range [list, optional] - option to artificially set colorbar range. By default it is set by colorbar dictionary ranges.
            - color_palett [dict, optional] - dictionary in which colors are defined to certain nozzles
            - connect_points [list, optional] - Option to connect points whose keys share a certain string.
            - connect_points_sort [str, optional] - Option to sort connect points by x and y values. Type "x" or "y" for the order of sorting.
            - connecting_lines_label [bool, optional] - option to display labels of lines connecting points. Labels are keys fo each line
            - connect_line_width [int, optional] - width of the connecting line
            - connect_line_color [dict, optional] - option to specify colors of connecting lines. Keys of the dict are connect line keys, values are colors.
            - curve_color [str or list, optional] - color of the line
            - curve_label [str or list, optional] - label of the line
            - curve_style [str or list, optional] - style of the line
            - curve_width [int or list, optional] - width of the line
            - curve_points [list, optional] - points of the curve. List should be given as [(x_points), (y_points)]. It is possible to plot multiple curves, in that case input multiple nested lists inside main list, format of nested lists should sty the same (first an array of x coordinates than y coordinates).
            - custom_legend [list, optional] - option to replace the legend with a custom one
            - edgecolor_pallet [dict,optional] - option to define the edgecolor of a symbol
            - filter_names [bool, optional] - option to filter points
            - figure_size [tuple, optional] - x and y dimensions of created figure (x, y)
            - fill_between [bool, optional] - option to fill between two or more curves defined with plot_curve.
            - fill_label [str, optional] - option to include fill_between in the legend
            - fill_alpha [float, optional] - Option to addapt the alpha (opacity) value of the fill
            - fill_curves_num [list, optional] - Option to select which list to fill between. Default all curves are used.
            - fill_color [str, optional] - Color of the fill between
            - iterations_split [str, optional] - character that differs iterations of data
            - legend_fontsize [int, optional] - font size of legend, x and y ticks.
            - legend_location [str, optional] - location of the label, eg. upper left.
            - legend_param [tuple, optional] - parameters of the legend box, tuple must be in order (x, y, width, height), width and height are optional.
            - legend_marker_size [int, optional] - option to set a custom size to the markers in a legend
            - marker_size [int, optional] - Size of the markers
            - minor_ticks [bool, optional] - option to plot both major and minor ticks
            - multiply_x [float, optional] - option to multiply all values on x axis by a certain coeficient (useful for converting from different units eg. meters to milimeters)
            - multiply_y [float, optional] - option to multiply all values on y axis by a certain coeficient (useful for converting from different units eg. meters to milimeters)
            - nozzle_parameters [dict] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values
            - n_legend_col [int, optional] - number of collumns the legend will take. By default it is set to 1
            - plot_curve [bool, optional] - mark True if you wish to plot a curve on the plot
            - plot_seperate [bool, optional] - option that by calling the function each time it generates a new plot in its own window
            - plot_font [str, optional] - font of all writing on the plot
            - point_opacity [float, optional] - option to make scatter points transparent
            - save_plot [str, optional] - option to save a plot into a specified directory
            - scale [str, optional] - scale the graph takes (can be log, linear, symlog, ...)
            - show_legend [bool, optional] - option to show/hide legend.
            - symbol_palett [dict, optional] -dictionary from which symbols are defined for nozzles
            - title_fontsize [int, optional] - font size of title
            - tick_fontsize [int, optional] - font size of axis ticks
            - x_font_style [str, optional] - font style of x axis label
            - x_lim [dict, optional] - dimensions of the plot in x-axis. Format is; {right : value, left : value}
            - x_t_minor_interval [float, optional] - interval between x axis minor ticks. If y scale is set to log the parameter defines base value
            - x_tick_interval [float, optional] - interval between x axis major ticks. If x scale is set ot log it defines base value
            - x_tick_axis_spacing [int, optional] - spacing between x axis and its ticks
            - y_font_style [str, optional] - font style of y axis label
            - y_lim [dict, optional] - dimensions of the plot in y-axis. Format is; {up : value, down : value}
            - y_t_minor_interval [float, optional] - interval between y axis minor ticks. If y scale is set to log the parameter defines base value
            - y_tick_interval [float, optional] - interval between y axis major ticks. If y scale is set to log the parameter defines base value
            - y_tick_axis_spacing [int, optional] - spacing between y axis and its ticks
    """
    #unpackign kwargs
    desired_nozzles = nozzles
    both = jdl.compare_lists(list1=list(values_y.keys()), list2=list(values_x.keys()))
    marker_size = kwargs.get("marker_size", 30)
    legend_marker_size = kwargs.get("legend_marker_size", marker_size)
    scale = kwargs.get("scale", "linear")
    scale_x = kwargs.get("scale_x", scale)
    scale_y = kwargs.get("scale_y", scale)
    iterations_split = kwargs.get("iterations_split", "_")
    nozzles_split = kwargs.get("nozzles_split", ".")
    xlim = kwargs.get("x_lim", None)
    ylim = kwargs.get("y_lim", None)
    title_fontsize = kwargs.get("title_fontsize", 50)
    axis_fontsize = kwargs.get("axis_fontsize", 30)
    legend_fontsize = kwargs.get("legend_fontsize", 20)
    legend_loc = kwargs.get("legend_location", "upper left")
    legend_param = kwargs.get("legend_param", (1.0, 1.0))
    plot_curve = kwargs.get("plot_curve", False)
    plot_seperate = kwargs.get("plot_seperate", True)
    point_opacity = kwargs.get("point_opacity", 1)
    curve_points = kwargs.get("curve_points", [(1,-1), (1, -1)])    
    c_label = kwargs.get("curve_label", "Line")
    c_color = kwargs.get("curve_color", "black")
    c_style = kwargs.get("curve_style", "--") 
    c_width = kwargs.get("curve_width", 1)
    edgecolor_pallet = kwargs.get("edgecolor_pallet", None)
    filter_names = kwargs.get("filter_names", True)
    figure_size = kwargs.get("figure_size", None)
    fill_between = kwargs.get("fill_between", False)
    fill_label = kwargs.get("fill_label", None)
    fill_alpha = kwargs.get("fill_alpha", 0.5)
    fill_curves_num = kwargs.get("fill_curves_num", None)
    fill_color = kwargs.get("fill_color", "gray")
    minor_ticks = kwargs.get("minor_ticks", False)
    multy_x = kwargs.get("multiply_x", 1)
    multy_y = kwargs.get("multiply_y", 1)
    nozzle_param = kwargs.get("nozzle_parameters", None)
    n_legend_col = kwargs.get("n_legend_col", 1)
    show_grid = kwargs.get("show_grid", True)
    x_tick_interval = kwargs.get("x_tick_interval")
    y_tick_interval = kwargs.get("y_tick_interval")
    x_t_minor_interval = kwargs.get("x_t_minor_interval")
    y_t_minor_interval = kwargs.get("y_t_minor_interval")
    plot_font = kwargs.get("plot_font", "Times New Roman")
    color_palett = kwargs.get("color_palett", None)
    symbol_palett = kwargs.get("symbol_palett", None)
    show_legend = kwargs.get("show_legend", True)
    save_plot = kwargs.get("save_plot", None)
    axis_scalar_format = kwargs.get("axis_scalar_format", True)
    custom_legend = kwargs.get("custom_legend", None)
    connect_points = kwargs.get("connect_points", None)
    connect_label = kwargs.get("connecting_lines_label", False)
    connect_line_width = kwargs.get("connect_line_width", 1)
    connect_line_color = kwargs.get("connect_line_color", None)
    connect_points_sort = kwargs.get("connect_points_sort", "x")
    tick_fontsize = kwargs.get("tick_fontsize", legend_fontsize)
    x_tick_axis_spacing = kwargs.get("x_tick_axis_spacing", 10)
    y_tick_axis_spacing = kwargs.get("y_tick_axis_spacing", 10)
    colorbar = kwargs.get("colorbar", None)
    colorbar_range = kwargs.get("colorbar_range", None)
    cbar_label = kwargs.get("colorbar_label", "Colorbar")

    #seting the font
    plt.rcParams["font.family"] = plot_font
    #plt.rcParams["mathtext.fontset"] = "cm"
    plt.rcParams['text.usetex'] = True
    plt.rcParams["text.latex.preamble"] = r"\usepackage{mathpazo}"  # Palatino font


                

    #creating a figure
    if plot_seperate:
        if figure_size is not None:
            plt.figure(figsize=figure_size)
        else:
            plt.figure()
    
    if nozzle_param is None:
        no_ratio = True
    else:
        no_ratio = False
        
    if colorbar is not None:
        max_colorbar = max(colorbar.values())
        min_colorbar = min(colorbar.values())
        if colorbar_range is not None:
            max_colorbar = colorbar_range[-1]
            min_colorbar = colorbar_range[0]
        # Create a normalization object for the color range
        norma = Normalize(vmin=min_colorbar, vmax=max_colorbar)  # Set the colorbar range
        
    #Fill between curves. Here, so that it is a part of the background)
    if fill_between and len(curve_points) > 2:
        fil_x = []
        fil_y1 = []
        fil_y2 = []
        curve_points_even = []
        if fill_curves_num is not None:
            for i in fill_curves_num:
                curve_points_even.append(curve_points[i])
        else:
            curve_points_even = curve_points.copy()
            
        if len(curve_points_even) % 2 != 0:
            curve_points_even = curve_points_even[:-1]

        for i in range(int((len(curve_points)-1)/2)):
            fil_x.append(curve_points_even[2*i][0])
            fil_y1.append(curve_points_even[2*i][1])
            fil_y2.append(curve_points_even[2*i+1][1])
            plt.fill_between(x = fil_x[0], y1 = fil_y1[0], y2 = fil_y2[0], color=fill_color, alpha=fill_alpha, edgecolor="none")
    
    #progress_bar initialization
    total_count = len(both)
    progress_bar = tqdm(total = total_count, desc="Processing")
    for i in desired_nozzles:
        if filter_names:
            points = jdl.filter_names(both, include=[f"{i}-"])
        else:
            points = both
        nozzle = i.split("-")[0]
        if not no_ratio:
            ratio = nozzle_param.get(i.split(iterations_split)[0].split(nozzles_split)[0])[-1]
            label = f"{nozzle} ; {ratio}"
        else:
            label = f"{nozzle}"
            
        if color_palett is None:
            color = random_pleasant_color()
        else:
            color = color_palett.get(i)
        if symbol_palett is None:
            symbol = random_markers()
        else:
            symbol = symbol_palett.get(i)
        if edgecolor_pallet is None:
            edge_c = "none"
        else:
            edge_c = edgecolor_pallet.get(i)
            if edge_c is None:
                edge_c = "none"
        c_xpoints = {}
        c_ypoints = {}
        all_c_points = []
        used_labels = []
        used_c_labels = set()
        used_c_colors = []
        
        # Plot main plot (scatter)
        for point in points:
            # Update the progress bar
            progress_bar.update()
            #update color
            if colorbar is not None:
                color = colorbar.get(point)
            if type(values_x.get(point)) != list:
                #ploting scattered points
                y = values_y.get(point)*multy_y
                x = values_x.get(point)*multy_x
                if label not in used_labels and custom_legend is None:
                    if colorbar is None:
                        main_plot = plt.scatter(x, y, label = label, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c)
                        used_labels.append(label)
                    else:
                        main_plot = plt.scatter(x, y, label = label, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c, norm=norma, cmap="jet")
                        used_labels.append(label)
                else:
                    if colorbar is None:
                        main_plot = plt.scatter(x, y, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c)
                    else:
                        main_plot = plt.scatter(x, y, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c, norm=norma, cmap="jet")

            if type(values_y.get(point)) == list:
                for j in range(len(values_x.get(point))):
                    #ploting scattered points
                    y = values_y.get(point)[j]*multy_y
                    x = values_x.get(point)[j]*multy_x
                    if label not in used_labels and custom_legend is None:
                        plt.scatter(x, y, label = label, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c)
                        used_labels.append(label)
                    else:
                         plt.scatter(x, y, c=color, marker = symbol, s=marker_size, alpha=point_opacity, edgecolor=edge_c)
            #getting x and y coordinates for connecting lines (if connect_points is True)
            if connect_points is not None:
                for c_point in connect_points:
                    if c_point in point:
                        if c_point not in all_c_points:
                            all_c_points.append(c_point)
                            c_xpoints.update({c_point : [x]})
                            c_ypoints.update({c_point : [y]})
                        else:
                            c_xpoints[c_point].extend([x])
                            c_ypoints[c_point].extend([y])
        if connect_points is not None:
            #sorting x and y points to connect so that lines are not jumbled
            if connect_points_sort is not None:
                if connect_points_sort == "x":
                    c_xpoints, index_change = jdl.sort_dict_values_by_elements(input_dict=c_xpoints)
                    c_ypoints = jdl.sort_dict_values_by_elements(input_dict=c_ypoints, sort_by_index=index_change)
                else:
                    c_ypoints, index_change = jdl.sort_dict_values_by_elements(input_dict=c_ypoints)
                    c_xpoints = jdl.sort_dict_values_by_elements(input_dict=c_xpoints, sort_by_index=index_change)
                    
            #sort labels by size
            num_val = []
            for p in all_c_points:
                num = jdl.strip_letters(p)
                if num != "":
                    num_val.append(int(num))
            if not num_val:
                index_list = np.argsort(np.argsort(num_val))
                reordered_list = all_c_points.copy()
                # Iterate through the index changes and update the reordered_list
                for i, change in enumerate(index_list):
                    reordered_list[change] = all_c_points[i]
            else:
                reordered_list = all_c_points.copy()
            #ploting connecting lines
            for c_point in reordered_list:
                #chosing line colors
                if connect_line_color is None:
                    lines_color = random_pleasant_color()
                    while lines_color in used_c_colors:
                        lines_color = random_pleasant_color()
                    used_c_colors.append(lines_color)
                else:
                    lines_color = connect_line_color[c_point]
                x = c_xpoints.get(c_point)
                y = c_ypoints.get(c_point)
                if connect_label is True and c_point not in used_c_labels:
                    plt.plot(x, y, label=c_point, c=lines_color, lw=connect_line_width)
                    used_c_labels.add(c_point)
                else:
                    plt.plot(x,y, c=lines_color, lw=connect_line_width)
    if colorbar is not None:
        cbar = plt.colorbar(main_plot)
        cbar.set_label(cbar_label, fontsize = axis_fontsize, labelpad=legend_fontsize/2, math_fontfamily='dejavuserif')
        cbar.ax.tick_params(labelsize=legend_fontsize)

    # Close the progress bar
    progress_bar.close()
    print("Main plot created sucesfully. Compiling optional parameters....")
    
    #___________________Optional parameters________________________
    
    #Extra curves
    if plot_curve:
        if type(curve_points[0]) is list or tuple:
            if type(c_label) is str:
                c_name = c_label[:]
                c_label = []
                for i in range(len(curve_points)):
                    c_label.append(c_name)
            if type(c_color) is str:
                c_c = c_color[:]
                c_color = []
                for i in range(len(curve_points)):
                    c_color.append(c_c)
            if type(c_style) is str:
                c_s = c_style[:]
                c_style = []
                for i in range(len(curve_points)):
                    c_style.append(c_s)
            if type(c_width) is int:
                c_w = copy.deepcopy(c_width)
                c_width = []
                for i in range(len(curve_points)):
                    c_width.append(c_w)
            for l in range(len(curve_points)):
                if c_label is not None:
                    plt.plot(curve_points[l][0], curve_points[l][1], label=c_label[l], color=c_color[l], ls=c_style[l], lw=c_width[l])
                else:
                    plt.plot(curve_points[l][0], curve_points[l][1], color=c_color[l], ls=c_style[l], lw=c_width[l])

        else:
            if c_label is not None:
                plt.plot(curve_points[0], curve_points[1], label=c_label, color=c_color, ls=c_style, lw=c_width)
            else:
                plt.plot(curve_points[0], curve_points[1], color=c_color, ls=c_style, lw=c_width)

    #Plot axis parameters
    
    plt.yscale(scale_y)
    plt.xscale(scale_x)
    ax = plt.gca()
    if axis_scalar_format:
        # Create a ScalarFormatter
        scalar_formatter = ScalarFormatter()
        scalar_formatter.set_scientific(False)  # Disable scientific notation
        # Apply the ScalarFormatter to the axis
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
    plt.title(title, fontsize=title_fontsize)
    if xlim is not None:
        plt.xlim(left = xlim.get("left"), right = xlim.get("right"))
    if ylim is not None:
        plt.ylim(ylim.get("down"), ylim.get("up"))
    plt.xlabel(x_label, fontsize=axis_fontsize, labelpad=axis_fontsize/2)
    plt.ylabel(y_label, fontsize=axis_fontsize, labelpad=axis_fontsize/2)
    if show_legend:
        if custom_legend is not None:
            for i in custom_legend:
                try:
                    edgec = i[4]
                except IndexError:
                    edgec = "none"
                try:
                    alpha_ = i[3]
                except IndexError:
                    alpha_ = 1
                plt.scatter(np.random.randint(0, 100), np.random.randint(0, 100), label = i[0], marker=i[1], color=i[2], alpha=alpha_, s=0, edgecolors=edgec)
        lgnd = plt.legend(fontsize=legend_fontsize, bbox_to_anchor=legend_param, loc=legend_loc, ncol=n_legend_col)
        if custom_legend is not None:
            skip = 0
            if plot_curve:
                if type(c_label) == str:
                    c_label_num = [0]
                if c_label is None:
                    c_label_num = []
                else:
                    c_label_num = range(len(c_label))
                for i in c_label_num:
                    skip +=1
            for i in range(len(custom_legend)+skip):
                lgnd.legendHandles[i]._sizes = [legend_marker_size]
        if fill_label is not None:
            handles, labels = lgnd.legendHandles, [text.get_text() for text in lgnd.get_texts()]
            edge_rgba = (1, 0.5, 0, 1)
            fill_patch = Patch(facecolor="gray", alpha=fill_alpha, label=fill_label, edgecolor = edge_rgba, linestyle="--")
            handles.append(fill_patch)
            labels.append(fill_label)
            plt.legend(handles=handles, labels=labels)

            
    plt.xticks(fontsize=tick_fontsize)
    plt.yticks(fontsize=tick_fontsize)
    plt.tick_params(axis='y', pad=y_tick_axis_spacing)
    plt.tick_params(axis='x', pad=x_tick_axis_spacing)
    # Show the major grid with a denser spacing
    if show_grid:
        plt.grid(True, which='major', linestyle='-', linewidth=0.5, color='black')
    # Set the tick interval for both x and y axes
    try:
        plt.gca().xaxis.set_major_locator(MultipleLocator(x_tick_interval))  # Adjust the interval as needed
        plt.gca().yaxis.set_major_locator(MultipleLocator(y_tick_interval))  # Adjust the interval as needed
    except TypeError:
        pass
    if minor_ticks and show_grid:
        plt.minorticks_on()
        # Show the minor grid with a less dense spacing
        plt.grid(True, which='minor', linestyle=':', linewidth=0.5, color='gray')
        try:
            plt.gca().xaxis.set_minor_locator(MultipleLocator(x_t_minor_interval))  # Adjust the interval as needed
            plt.gca().yaxis.set_minor_locator(MultipleLocator(y_t_minor_interval))  # Adjust the interval as needed
        except TypeError:
            pass
    
    if save_plot is not None:
        plt.savefig(save_plot, bbox_inches = 'tight', transparent = True, dpi=500)
    
    print("Plot created sucesfully.")
    
    plt.show()

def plot_xyz_data(nozzles, values_x, values_y, values_z, x_label, y_label, z_label, title, **kwargs):
    """
    Plots three sets of data as x, y and z coordinates.
    
    Input parameters:
        - nozzles [list] - List of nozzle names
        - values_y [dict] - y values on graph. Dictionary is writen; {base_name : value}
        - values_x [dict] - x values on graph. Dictionary is writen; {base_name : value}
        - values_z [dict] - z values on graph. Dictionary is writen; {base_name : value}
        - x_label [str] - Label of x axis
        - y_label [str] - Label of y axis
        - z_label [str] - Lbael of z axis
        - title [str] - Title of the graph
        - kwargs:
            - axis_fontsize [int, optional] - font size of x and y axis labels
            - axis_scalar_format [bool, optional] - option to switch between scientific and normal tick display
            - color_palett [dict, optional] - dictionary in which colors are defined to certain nozzles
            - iterations_split [str, optional] - character that differs iterations of data
            - legend_fontsize [int, optional] - font size of legend, x and y ticks.
            - surface_color [str, optional] - color of the surface
            - surface_label [str, optional] - label of surface
            - surface_transparency [int, optional] - transparency of the surface
            - surface_points [list, optional] - at least three different points for a trisurface to be ploted on. Positions should be given as [(x_first_point, y_first_point, z_first_point), (x_second_point,....]
            - marker_size [int, optional] - Size of the markers
            - minor_ticks [bool, optional] - option to plot both major and minor ticks
            - nozzles_split [str, optional] - character that differs nozzles with the same geometrical values
            - nozzle_parameters [dict, optional] - geometrical parameters of nozzles. Dictionary is writen - {base_name:nozzle_parameter}
            - plot_line [bool, optional] - mark True if you wish to plot a line on th plot
            - symbol_palett [dict, optional] -dictionary from which symbols are defined for nozzles
            - title_fontsize [int, optional] - font size of title
            - x_scale [str, optional] - scale the x axis takes (can be log, linear, symlog, ...)
            - x_font_style [str, optional] - font style of x axis label
            - x_lim [dict, optional] - dimensions of the plot in x-axis. Format is; {x1 : value, x2 : value}
            - x_t_minor_interval [float, optional] - interval between x axis minor ticks. If y scale is set to log the parameter defines base value
            - x_tick_interval [float, optional] - interval between x axis major ticks. If x scale is set ot log it defines base value
            - y_scale [str, optional] - scale the y axis takes (can be log, linear, symlog, ...)
            - y_font_style [str, optional] - font style of y axis label
            - y_lim [dict, optional] - dimensions of the plot in y-axis. Format is; {y1 : value, y2 : value}
            - y_t_minor_interval [float, optional] - interval between y axis minor ticks. If y scale is set to log the parameter defines base value
            - y_tick_interval [float, optional] - interval between y axis major ticks. If y scale is set to log the parameter defines base value
            - z_scale [str, optional] - scale the z axis takes (can be log, linear, symlog, ...)
            - z_font_style [str, optional] - font style of z axis label
            - z_lim [dict, optional] - dimensions of the plot in x-axis. Format is; {z1 : value, z2 : value}
            - z_t_minor_interval [float, optional] - interval between z axis minor ticks. If z scale is set to log the parameter defines base value
            - z_tick_interval [float, optional] - interval between z axis major ticks. If z scale is set to log the parameter defines base value
    """
    nozzle_parameters = kwargs.get("nozzle_parameters", None)
    marker_size = kwargs.get("marker_size", 30)
    x_scale = kwargs.get("x_scale", "linear")
    y_scale = kwargs.get("y_scale", "linear")
    z_scale = kwargs.get("z_scale", "linear")
    iterations_split = kwargs.get("iterations_split", "_")
    nozzles_split = kwargs.get("nozzles_split", ".")
    xlim = kwargs.get("x_lim", None)
    ylim = kwargs.get("y_lim", None)
    zlim = kwargs.get("z_lim", None)
    title_fontsize = kwargs.get("title_fontsize", 50)
    axis_fontsize = kwargs.get("axis_fontsize", 30)
    legend_fontsize = kwargs.get("legend_fontsize", 20)
    s_points = kwargs.get("surface_points", [(1,1,1), (-1,-1,-1), (1,-1,1)])    
    s_label = kwargs.get("surface_label", "Line")
    s_color = kwargs.get("surface_color", "black")
    s_transp = kwargs.get("surface_transparency", 1)
    minor_ticks = kwargs.get("minor_ticks", True)
    x_tick_interval = kwargs.get("x_tick_interval")
    y_tick_interval = kwargs.get("y_tick_interval")
    z_tick_interval = kwargs.get("z_tick_interval")
    x_t_minor_interval = kwargs.get("x_t_minor_interval")
    y_t_minor_interval = kwargs.get("y_t_minor_interval")
    z_t_minor_interval = kwargs.get("z_t_minor_interval")
    x_font_style = kwargs.get("x_font_style", "italic")
    y_font_style = kwargs.get("y_font_style", "italic")
    z_font_style = kwargs.get("z_font_style", "italic")
    title_font_style = kwargs.get("title_font_style", "italic")   
    color_palett = kwargs.get("color_palett", None)
    symbol_palett = kwargs.get("symbol_palett", None)
    axis_scalar_format = kwargs.get("axis_scalar_format", True)
    all_names = jdl.compare_lists(list1 = list(values_x), list2 = jdl.compare_lists(list(values_y), list(values_z)))
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for i in nozzles:
        points = jdl.filter_names(all_names, include=[f"{i}"])
        if nozzle_parameters is not None:
            ratio = nozzle_parameters.get(i.split(iterations_split)[0].split(nozzles_split)[0])[-1]
            nozzle = i.split("-")[0]
            label = f"{nozzle} ; {ratio}"
        else:
            nozzle = i.split("-")[0]
            label = f"{nozzle}"
        if color_palett is None:
            color = random_pleasant_color()
        else:
            color = color_palett.get(i)
        if symbol_palett is None:
            symbol = random_markers()
        else:
            symbol = symbol_palett.get(i)
        used_labels = []
        for point in points:
            x = values_x.get(point)
            y = values_y.get(point)
            z = values_z.get(point)
            if label not in used_labels:
                ax.scatter(x, y, z, label = label, c=color, marker = symbol, s=marker_size)
                used_labels.append(label)
            else:
                ax.scatter(x, y, z, c=color, marker = symbol, s=marker_size)
    if kwargs.get("plot_surface", False):
        ax.plot_trisurf(s_points[0], s_points[1], s_points[2], label=s_label, alpha=s_transp, color=s_color)
    #Customize plot
    ax.set_yscale(x_scale)
    ax.set_xscale(y_scale)
    ax.set_zscale(z_scale)
    if axis_scalar_format:
        # Create a ScalarFormatter
        scalar_formatter = ScalarFormatter()
        scalar_formatter.set_scientific(False)  # Disable scientific notation
        # Apply the ScalarFormatter to the axis
        for axis in [ax.xaxis, ax.yaxis]:
            axis.set_major_formatter(ScalarFormatter())
    plt.title(title, fontsize=title_fontsize, fontstyle=title_font_style)
    if xlim is not None:
        ax.set_xlim(left = xlim.get("x1"), right = xlim.get("x2"))
    if ylim is not None:
        ax.set_.ylim(ylim.get("y1"), ylim.get("y2"))
    if zlim is not None:
        ax.set_zlim(zlim.get("z1"), zlim.get("z2"))
    ax.set_xlabel(x_label, fontsize=axis_fontsize, fontstyle=x_font_style)
    ax.set_ylabel(y_label, fontsize=axis_fontsize, fontstyle=y_font_style)
    ax.set_zlabel(z_label, fontsize=axis_fontsize, fontstyle=z_font_style)
    #ax.legend(fontsize=legend_fontsize)
    ax.xaxis.label.set_fontsize(legend_fontsize)  # X-axis label
    ax.yaxis.label.set_fontsize(legend_fontsize)  # Y-axis label
    ax.zaxis.label.set_fontsize(legend_fontsize)  # Z-axis label
    # Show the major grid with a denser spacing
    ax.grid(True, which='major', linestyle='-', linewidth=0.5, color='black')
    # Set the tick interval for both x and y axes
    try:
        plt.gca().xaxis.set_major_locator(MultipleLocator(x_tick_interval))  # Adjust the interval as needed
        plt.gca().yaxis.set_major_locator(MultipleLocator(y_tick_interval))  # Adjust the interval as needed
        plt.gca().zaxis.set_major_locator(MultipleLocator(z_tick_interval))  # Adjust the interval as needed
    except TypeError:
        pass
    if minor_ticks:
        plt.minorticks_on()
        # Show the minor grid with a less dense spacing
        plt.grid(True, which='minor', linestyle=':', linewidth=0.5, color='gray')
        try:
            plt.gca().xaxis.set_minor_locator(MultipleLocator(x_t_minor_interval))  # Adjust the interval as needed
            plt.gca().yaxis.set_minor_locator(MultipleLocator(y_t_minor_interval))  # Adjust the interval as needed
            plt.gca().zaxis.set_minor_locator(MultipleLocator(z_t_minor_interval))  # Adjust the interval as needed
        except TypeError:
            pass
    plt.show()
    
    
    
    
    
    
    
    
    