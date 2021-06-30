import networkx as nx

import numpy as np
import pandas as pd

from bokeh.io import show, output_file, export_png
from bokeh.models import (Plot, Range1d, MultiLine, Circle, TapTool, OpenURL, HoverTool, 
                          CustomJS, Slider, Column, StaticLayoutProvider, TapTool, 
                          WheelZoomTool, PanTool, ResetTool, SaveTool, FixedTicker, 
                          LinearColorMapper, LogColorMapper, ColorBar, BasicTicker, LogTicker)
from bokeh.plotting import figure, from_networkx

from matplotlib import cm
from matplotlib.colors import to_hex


# convert a matrix of real values into a matrix of hex colors
def create_color_matrix(color_of_landmarks, my_palette):
    MAX_VALUE = color_of_landmarks.max()
    MIN_VALUE = color_of_landmarks.min()

    color_matrix = []

    for row in color_of_landmarks:
        color_matrix.append([ to_hex(my_palette((value - MIN_VALUE) / (MAX_VALUE - MIN_VALUE))) for value in row ])

    return color_matrix


def sample_size_and_color_matrix(distance_from_landmarks, color_of_landmarks,
                                 my_palette,
                                 num_of_points):
    MAX_VALUE = distance_from_landmarks.max()

    MAX_COLOR = color_of_landmarks.max()
    MIN_COLOR = color_of_landmarks.min()

    size_matrix = []
    color_matrix = []

    for i, row in enumerate(distance_from_landmarks):
        add_size = []
        add_color = []
        for x in np.linspace(0, MAX_VALUE, num=num_of_points):
            index_list = np.where(row<=x)[0]
            add_size.append(len(index_list))

            add_color.append(to_hex(my_palette((color_of_landmarks[i][index_list[-1]] - MIN_COLOR) / (MAX_COLOR - MIN_COLOR))))

        size_matrix.append(add_size)
        color_matrix.append(add_color)


    return size_matrix, color_matrix



def plot_multiscale_BM(d_matrix, 
                       distance_from_landmarks,
                       color_of_landmarks,
                       my_palette = cm.get_cmap(name='jet'),
                       log_scaling=False):
    
    color_matrix = create_color_matrix(color_of_landmarks, my_palette)
    
    # create a graph from the distance matrix
    G = nx.from_numpy_matrix(d_matrix)
    NODE_SIZE = 10
    
    for edge in G.edges:
        G.edges[edge]['alpha'] = 0

    for i, node in enumerate(G.nodes):
        G.nodes[node]['size_scale'] = NODE_SIZE
        G.nodes[node]['size'] = NODE_SIZE
        G.nodes[node]['distance_from_landmark'] = distance_from_landmarks[i].tolist()
        G.nodes[node]['color_list'] = color_matrix[i]
        G.nodes[node]['color'] = color_matrix[i][0]
    
    
    plot = Plot(#plot_width=700, plot_height=700,
                x_range=Range1d(-2, 2), y_range=Range1d(-2, 2),
                sizing_mode="stretch_both",
                toolbar_location = 'right')

    zoom_tool = WheelZoomTool()
    plot.add_tools(PanTool(), zoom_tool,
                        ResetTool(), SaveTool())
    plot.toolbar.active_scroll = zoom_tool

    

    min_epsilon = 0
    max_epsilon = distance_from_landmarks.max()


    graph_renderer = from_networkx(G, nx.kamada_kawai_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size='size', fill_color='color', fill_alpha=0.5)
    graph_renderer.edge_renderer.glyph = MultiLine(line_color = "black", line_alpha = 'alpha', line_width = 2)

    plot.renderers.append(graph_renderer)

    code = """ 
        var current_epsilon = cb_obj.value;

        var node_data = graph_renderer.node_renderer.data_source.data;
        var edge_data = graph_renderer.edge_renderer.data_source.data;

        for (var i = 0; i < node_data['size'].length; i++) {
            var counter = 0;
            for (var j = 0; j < node_data['distance_from_landmark'][i].length; j++) {
                if (node_data['distance_from_landmark'][i][j]  > current_epsilon) { break; }
                counter++;
            }

            if (log_scaling) {
                graph_renderer.node_renderer.data_source.data['size'][i] = Math.max(10,Math.min(Math.log10(counter)*node_data['size_scale'][i], 100));

            }
            else {
                graph_renderer.node_renderer.data_source.data['size'][i] = counter*node_data['size_scale'][i];
            }

            graph_renderer.node_renderer.data_source.data['color'][i] = graph_renderer.node_renderer.data_source.data['color_list'][i][counter-1];
        }

        for (var i = 0; i < edge_data['alpha'].length; i++) {
            graph_renderer.edge_renderer.data_source.data['alpha'][i] = + (edge_data['weight'][i] <= current_epsilon);
        }

        graph_renderer.node_renderer.data_source.change.emit();
        graph_renderer.edge_renderer.data_source.change.emit();
    """



    callback = CustomJS(args = dict(graph_renderer = graph_renderer,
                                    log_scaling = log_scaling),
                        code = code)
    slider = Slider(start=min_epsilon, end=max_epsilon, step=0.1, value=min_epsilon,
                   title='epsilon',
                   )
    slider.js_on_change('value', callback)

    # # continuous colorbar
    num_ticks = 100
    low = color_of_landmarks.min()
    high = color_of_landmarks.max()

    color_mapper = LinearColorMapper(palette=[to_hex(my_palette(color_id)) 
                                              for color_id in np.linspace(0, 1, num_ticks)], 
                                     low=low, high=high)

    color_bar = ColorBar(color_mapper=color_mapper, 
                         major_label_text_font_size='14pt',
                         label_standoff=12,
                        )

    plot.add_layout(color_bar, 'right')

    layout = Column(plot, slider, sizing_mode="scale_both")
    
    return layout




def plot_sampled_multiscale_BM(d_matrix, 
                               distance_from_landmarks,
                               color_of_landmarks,
                               num_of_points = 100,
                               my_palette = cm.get_cmap(name='jet'),
                               log_scaling=False):
    
    sampled_size_matrix, sampled_color_matrix = sample_size_and_color_matrix(distance_from_landmarks, 
                                                                             color_of_landmarks,
                                                                             my_palette,
                                                                             num_of_points)
    G = nx.from_numpy_matrix(d_matrix)
    NODE_SIZE = 10

    for edge in G.edges:
        G.edges[edge]['alpha'] = 0

    for i, node in enumerate(G.nodes):
        G.nodes[node]['size_scale'] = NODE_SIZE
        G.nodes[node]['size'] = NODE_SIZE
        G.nodes[node]['size_list'] = sampled_size_matrix[i]
        G.nodes[node]['filtration_list'] = np.linspace(0, distance_from_landmarks.max(), num=num_of_points)
        G.nodes[node]['color_list'] = sampled_color_matrix[i]
        G.nodes[node]['color'] = sampled_color_matrix[i][0]

    
    plot = Plot(
        #plot_width=700, plot_height=700,
                x_range=Range1d(-2, 2), y_range=Range1d(-2, 2),
                sizing_mode="stretch_both",
                toolbar_location = 'right')

    zoom_tool = WheelZoomTool()
    plot.add_tools(PanTool(), zoom_tool,
                        ResetTool(), SaveTool())
    plot.toolbar.active_scroll = zoom_tool

    min_epsilon = 0
    max_epsilon = distance_from_landmarks.max()


    graph_renderer = from_networkx(G, nx.kamada_kawai_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size='size', fill_color='color', fill_alpha=0.5)
    graph_renderer.edge_renderer.glyph = MultiLine(line_color = "black", line_alpha = 'alpha', line_width = 2)

    plot.renderers.append(graph_renderer)

    code = """ 
        var current_epsilon = cb_obj.value;

        var node_data = graph_renderer.node_renderer.data_source.data;
        var edge_data = graph_renderer.edge_renderer.data_source.data;

        for (var i = 0; i < node_data['size'].length; i++) {
            var counter = 0;
            for (var j = 0; j < node_data['filtration_list'][i].length; j++) {
                if (node_data['filtration_list'][i][j]  > current_epsilon) { break; }
                counter++;
            }

            if (log_scaling) {
                graph_renderer.node_renderer.data_source.data['size'][i] = Math.max(10,Math.min(Math.log10(node_data['size_list'][i][counter])*node_data['size_scale'][i], max_size));
            }
            else {
                graph_renderer.node_renderer.data_source.data['size'][i] = node_data['size_list'][i][counter]*node_data['size_scale'][i];
            }

            graph_renderer.node_renderer.data_source.data['color'][i] = graph_renderer.node_renderer.data_source.data['color_list'][i][counter-1];
        }

        for (var i = 0; i < edge_data['alpha'].length; i++) {
            graph_renderer.edge_renderer.data_source.data['alpha'][i] = + (edge_data['weight'][i] <= current_epsilon);
        }

        graph_renderer.node_renderer.data_source.change.emit();
        graph_renderer.edge_renderer.data_source.change.emit();


    """
    callback = CustomJS(args = dict(graph_renderer = graph_renderer,
                                    log_scaling = True,
                                    max_size = 100),
                        code = code)
    slider = Slider(start=min_epsilon, end=max_epsilon, step=0.1, value=min_epsilon,
                   title='epsilon',
                   )
    slider.js_on_change('value', callback)

    # # continuous colorbar
    num_ticks = 100
    low = color_of_landmarks.min()
    high = color_of_landmarks.max()

    color_mapper = LinearColorMapper(palette=[to_hex(my_palette(color_id)) 
                                              for color_id in np.linspace(0, 1, num_ticks)], 
                                     low=low, high=high)

    color_bar = ColorBar(color_mapper=color_mapper, 
                         major_label_text_font_size='14pt',
                         label_standoff=12,
                        )

    plot.add_layout(color_bar, 'right')

    layout = Column(plot, slider, sizing_mode="scale_both")

    return layout





def save_graph_to_png(epsilon_list, driver, 
                      d_matrix,
                      distance_from_landmarks,
                      color_of_landmarks,
                      plot_width = 500,
                      plot_height = 500,
                      filename = 'img/plot', 
                      add_colorbar = False,
                      num_of_points = 1000,
                      my_palette = cm.get_cmap(name='jet'),
                      log_scaling=False):
    
    sampled_size_matrix, sampled_color_matrix = sample_size_and_color_matrix(distance_from_landmarks, 
                                                                             color_of_landmarks,
                                                                             my_palette,
                                                                             num_of_points)
    G = nx.from_numpy_matrix(d_matrix)
    NODE_SIZE = 10

    for edge in G.edges:
        G.edges[edge]['alpha'] = 0

    for i, node in enumerate(G.nodes):
        G.nodes[node]['size_scale'] = NODE_SIZE
        G.nodes[node]['size'] = NODE_SIZE
        G.nodes[node]['size_list'] = sampled_size_matrix[i]
        G.nodes[node]['filtration_list'] = np.linspace(0, distance_from_landmarks.max(), num=num_of_points)
        G.nodes[node]['color_list'] = sampled_color_matrix[i]
        G.nodes[node]['color'] = sampled_color_matrix[i][0]
        
    plot = Plot(plot_width=plot_width, plot_height=plot_height,
                x_range=Range1d(-2, 2), y_range=Range1d(-2, 2),
                sizing_mode="fixed",
                toolbar_location = None)

    graph_renderer = from_networkx(G, nx.kamada_kawai_layout, scale=1, center=(0, 0))

    graph_renderer.node_renderer.glyph = Circle(size='size', fill_color='color', fill_alpha=0.5)
    graph_renderer.edge_renderer.glyph = MultiLine(line_color = "black", line_alpha = 'alpha', line_width = 2)

    
    plot.renderers.append(graph_renderer)


    # # continuous colorbar
    if add_colorbar:
        num_ticks = 100
        low = color_of_landmarks.min()
        high = color_of_landmarks.max()

        color_mapper = LinearColorMapper(palette=[to_hex(my_palette(color_id)) 
                                                  for color_id in np.linspace(0, 1, num_ticks)], 
                                         low=low, high=high)

        color_bar = ColorBar(color_mapper=color_mapper, 
                             major_label_text_font_size='14pt',
                             label_standoff=12,
                            )

        plot.add_layout(color_bar, 'right')

    
    for current_epsilon in epsilon_list:
        for i, node in enumerate(G.nodes):
            nearest_index = np.where(G.nodes[node]['filtration_list']<=current_epsilon)[0][-1]

            if log_scaling:
                graph_renderer.node_renderer.data_source.data['size'][i] = max(NODE_SIZE, 
                                                                       min(np.log10(G.nodes[node]['size_list'][nearest_index])*G.nodes[node]['size_scale'], 100))
            else:
                graph_renderer.node_renderer.data_source.data['size'][i] = G.nodes[node]['size_list'][nearest_index]*G.nodes[node]['size_scale']
            graph_renderer.node_renderer.data_source.data['color'][i] = G.nodes[node]['color_list'][nearest_index];


        for i, edge in enumerate(G.edges):
            graph_renderer.edge_renderer.data_source.data['alpha'][i] = int(G.edges[edge]['weight'] <= current_epsilon)

        export_png(plot, filename="{}_{}.png".format(filename, current_epsilon),
                  webdriver=driver)
