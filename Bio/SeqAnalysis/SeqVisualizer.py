#!/usr/bin/env python3
from typing import List, Dict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


class SeqVisualizer:
    """
    A class, that takes in two dictionaries from SeqAnalysis.SeqAnalyzer and plots selected values (parameters) in
    chosen type of plot. User can specify if plot should be saved or shown
    Possible plot types:
        - bar - Representing values of chosen parameter from dictionary as bar type plot. Can plot values for all sequences
         together and values for single chosen sequence from dictionary.
            Possible parameters to plot:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)
            - Secondary structure fraction
            - Extinction coefficient
            - Amino acid composition
            - Amino acid percent composition
            - Frequency of C terminus
            - Frequency of N terminus

        - bar_interactive - Representing values of chosen parameter from dictionary as interactive bar type plot. Can plot
         values for all sequences together and values for single chosen sequence from dictionary.
            Possible parameters to plot:
            Same as bar


        - box - Representing values of chosen parameter from dictionary as box type plot. Plotting values only for all
         sequences together from dictionary.
            Possible parameters to plot:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)
            - Secondary structure fraction
            - Extinction coefficient
            - Amino acid composition
            - Amino acid percent composition

        - scatter - Representing values of chosen parameters from dictionary as scatter plot (3 parameters can be passed).
        Plotting values only for all sequences together and values for single chosen sequence from dictionary.
            Possible parameters to plot:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)

    """

    def __init__(self, one_big_dict, dict_of_many_small_dicts):
        """
        Crucial for proper functioning of function is supplying a proper format of dictionaries. Key names in
        one_big_dict should be the same as given below. In dict_of_many_small_dicts

        :param one_big_dict: takes in dictionary with values calculated for all sequences
            For proper functioning of class, specific format of dictionaries is required (names of keys must be same):

                one_big_dict = {
                    'id': ['BX5012', 'PH3210'],
                    'sequence_length': [21, 35],
                    'molecular_weight': [2399.853, 4048.622],
                    'theoretical_pi': [10.461, 9.9],
                    'instability_index': [42.448, 58.869],
                    'flexibility': [1.014, 1.043],
                    'gravy': [-0.752, -1.166],
                    'aromaticity': [0.048, 0.086],
                    'secondary_structure_fraction': [{'Helix': 23.81, 'Turn': 28.571, 'Sheet': 19.048},
                                                     {'Helix': 22.857, 'Turn': 25.714, 'Sheet': 5.714}],
                    'extinction_coefficient': [{'Cysteines_reduced': 1, 'Cysteines_residues': 5},
                                               {'Cysteines_reduced': 2, 'Cysteines_residues': 3}],
                    'amino_acid_comp': [
                        {'A': 1, 'C': 0, 'D': 0, 'E': 1, 'F': 1, 'G': 1, 'H': 0, 'I': 2, 'K': 4, 'L': 1, 'M': 1, 'N': 1,
                         'P': 2, 'Q': 2, 'R': 1, 'S': 2, 'T': 0, 'V': 1, 'W': 0, 'Y': 0},
                        {'A': 0, 'C': 0, 'D': 5, 'E': 1, 'F': 3, 'G': 1, 'H': 0, 'I': 2, 'K': 9, 'L': 0, 'M': 1, 'N': 0,
                         'P': 1, 'Q': 0, 'R': 2, 'S': 7, 'T': 0, 'V': 3, 'W': 0, 'Y': 0}],
                    'amino_acid_percent': [
                        {'A': 4.762, 'C': 0.0, 'D': 0.0, 'E': 4.762, 'F': 4.762, 'G': 4.762, 'H': 0.0, 'I': 9.524,
                         'K': 19.048, 'L': 4.762, 'M': 4.762, 'N': 4.762, 'P': 9.524, 'Q': 9.524, 'R': 4.762,
                         'S': 9.524, 'T': 0.0, 'V': 4.762, 'W': 0.0, 'Y': 0.0},
                        {'A': 0.0, 'C': 0.0, 'D': 14.286, 'E': 2.857, 'F': 8.571, 'G': 2.857, 'H': 0.0, 'I': 5.714,
                         'K': 25.714, 'L': 0.0, 'M': 2.857, 'N': 0.0, 'P': 2.857, 'Q': 0.0, 'R': 5.714, 'S': 20.0,
                         'T': 0.0, 'V': 8.571, 'W': 0.0, 'Y': 0.0}],
                    'c_term_freq': {'M': 0.5, 'I': 0.5},
                    'n_term_freq': {'R': 0.5, 'A': 0.5}
                }

        :param dict_of_many_small_dicts: takes in dictionary of dictionaries similar to one_big_dict, but each
        dictionary stores data just for one sequence:

                dict_of_many_small_dicts = {
                    'BX5012': {
                        'id': ['BX5012'],
                        'sequence_length': [43],
                        'molecular_weight': [4899.524],
                        'theoretical_pi': [6.516],
                        'instability_index': [35.179], 'flexibility': [0.996], 'gravy': [-0.749],
                        'aromaticity': [0.093],
                        'secondary_structure_fraction': [{'Helix': 18.605, 'Turn': 16.279, 'Sheet': 25.581}],
                        'extinction_coefficient': [{'Cysteines_reduced': 8480, 'Cysteines_residues': 8730}],
                        'amino_acid_comp': [{'A': 4, 'C': 4, 'D': 3, 'E': 3, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4,
                            'L': 2, 'M': 2, 'N': 2, 'P': 1, 'Q': 2, 'R': 2, 'S': 1, 'T': 3, 'V': 2, 'W': 1, 'Y': 2}],
                        'amino_acid_percent': [
                               {'A': 9.302, 'C': 9.302, 'D': 6.977, 'E': 6.977, 'F': 2.326, 'G': 6.977, 'H': 2.326,
                                'I': 0.0, 'K': 9.302, 'L': 4.651, 'M': 4.651, 'N': 4.651, 'P': 2.326, 'Q': 4.651,
                                'R': 4.651, 'S': 2.326, 'T': 6.977, 'V': 4.651, 'W': 2.326, 'Y': 4.651}],
                        'c_term_freq': {'G': 1},
                        'n_term_freq': {'M': 1}
                        },
                    'PH3210': {
                        'id': ['PH3210'],
                        'sequence_length': [35],
                        'molecular_weight': [4048.622],
                        'theoretical_pi': [9.9],
                        'instability_index': [58.869],
                        'flexibility': [1.043],
                        'gravy': [-1.166],
                        'aromaticity': [0.086],
                        'secondary_structure_fraction': [{'Helix': 22.857, 'Turn': 25.714, 'Sheet': 5.714}],
                        'extinction_coefficient': [{'Cysteines_reduced': 0, 'Cysteines_residues': 0}],
                        'amino_acid_comp': [
                               {'A': 0, 'C': 0, 'D': 5, 'E': 1, 'F': 3, 'G': 1, 'H': 0, 'I': 2, 'K': 9, 'L': 0,
                                'M': 1, 'N': 0, 'P': 1, 'Q': 0, 'R': 2, 'S': 7, 'T': 0, 'V': 3, 'W': 0, 'Y': 0}],
                        'amino_acid_percent': [{'A': 0.0, 'C': 0.0, 'D': 14.286, 'E': 2.857, 'F': 8.571, 'G': 2.857,
                            'H': 0.0, 'I': 5.714, 'K': 25.714, 'L': 0.0, 'M': 2.857, 'N': 0.0, 'P': 2.857, 'Q': 0.0,
                            'R': 5.714, 'S': 20.0, 'T': 0.0, 'V': 8.571, 'W': 0.0, 'Y': 0.0}],
                        'c_term_freq': {'R': 1},
                        'n_term_freq': {'M': 1}
                        }
                }
        """

        self.main_dict = one_big_dict  # input dictionary
        self.sub_dict = dict_of_many_small_dicts
        self.ids = one_big_dict['id']  # list of sequences ids
        self.lists_args = ['sequence_length', 'molecular_weight', 'theoretical_pi', 'aromaticity', 'instability_index',
                           'flexibility', 'gravy']
        self.dicts_args = ['c_term_freq', 'n_term_freq']
        self.list_of_dicts_args = ['secondary_structure_fraction', 'extinction_coefficient', 'amino_acid_comp',
                                   'amino_acid_percent']
        self.all_args = self.lists_args + self.dicts_args + self.list_of_dicts_args
        self.box_args = self.lists_args + self.list_of_dicts_args
        self.dict_for_ax_units = {
            'id': 'ID',
            'sequence_length': 'Sequence length [Quantity]',
            'molecular_weight': 'Molecular weight [Da]',
            'theoretical_pi': 'Theoretical isoelectric point [pH]',
            'aromaticity': 'Aromaticity [Phe+Trp+Tyr/all aa]',
            'instability_index': 'Instability index',
            'flexibility': 'Flexibility',
            'gravy': 'Average hydrophobicity and hydrophilicity',
            'c_term_freq': 'Frequency of C terminus [amino acid]',
            'n_term_freq': 'Frequency of N terminus [amino acid]',
            'secondary_structure_fraction': 'Secondary structure fraction [%]',
            'extinction_coefficient': 'Extinction coefficient [m^3/(mol*cm)]',
            'amino_acid_comp': 'Amino acid composition [Quantity]',
            'amino_acid_percent': 'Amino acid percent composition [%]'
        }
        self.dict_for_title_bar = {
            'id': 'ID',
            'sequence_length': 'Comparision of lengths of sequences',
            'molecular_weight': 'Comparision of molecular weights of sequences',
            'theoretical_pi': 'Comparision of theoretical isolectric points of sequences',
            'aromaticity': 'Comparision of aromaticities of sequences',
            'instability_index': 'Comparision of instability indexes of sequences',
            'flexibility': 'Comparision of flexibility of sequences',
            'gravy': 'Comparision of hydrophobicity and hydrophilicity of sequences',
            'c_term_freq': 'Frequency of amino acids on C terminus',
            'n_term_freq': 'Frequency of amino acids on N terminus',
            'secondary_structure_fraction': 'Comparision of secondary structure values of sequences',
            'extinction_coefficient': 'Comparision of extinction coefficient values of sequences',
            'amino_acid_comp': 'Comparision of amino acid composition of sequences',
            'amino_acid_percent': 'Comparision of percent amino acid composition of sequences'
        }
        self.dict_for_title_bar = {
            'id': 'ID',
            'sequence_length': 'Comparision of lengths of sequences',
            'molecular_weight': 'Comparision of molecular weights of sequences',
            'theoretical_pi': 'Comparision of theoretical isolectric points of sequences',
            'aromaticity': 'Comparision of aromaticities of sequences',
            'instability_index': 'Comparision of instability indexes of sequences',
            'flexibility': 'Comparision of flexibility of sequences',
            'gravy': 'Comparision of hydrophobicity and hydrophilicity of sequences',
            'c_term_freq': 'Frequency of amino acids on C terminus',
            'n_term_freq': 'Frequency of amino acids on N terminus',
            'secondary_structure_fraction': 'Comparision of secondary structure values of sequences',
            'extinction_coefficient': 'Comparision of extinction coefficient values of sequences',
            'amino_acid_comp': 'Comparision of amino acid composition of sequences',
            'amino_acid_percent': 'Comparision of percent amino acid composition of sequences'
        }
        self.dict_for_x_ax_bar_box = {
            'id': 'ID',
            'sequence_length': 'Sequences',
            'molecular_weight': 'Sequences',
            'theoretical_pi': 'Sequences',
            'aromaticity': 'Sequences',
            'instability_index': 'Sequences',
            'flexibility': 'Sequences',
            'gravy': 'Sequences',
            'c_term_freq': 'Amino acids',
            'n_term_freq': 'Amino acids',
            'secondary_structure_fraction': 'Secondary structures',
            'extinction_coefficient': 'Cystein groups',
            'amino_acid_comp': 'Amino acids',
            'amino_acid_percent': 'Amino acids'
        }
        self.dict_for_title_box = {
            'id': 'ID',
            'sequence_length': 'Distribution of lengths of sequences',
            'molecular_weight': 'Distribution of molecular weights of sequences',
            'theoretical_pi': 'Distribution of theoretical isolectric points of sequences',
            'aromaticity': 'Distribution of aromaticities of sequences',
            'instability_index': 'Distribution of instability indexes of sequences',
            'flexibility': 'Distribution of flexibility of sequences',
            'gravy': 'Distribution of hydrophobicity and hydrophilicity of sequences',
            'secondary_structure_fraction': 'Distribution of secondary structure values of sequences',
            'extinction_coefficient': 'Distribution of extinction coefficient values of sequences',
            'amino_acid_comp': 'Distribution of amino acid composition of sequences',
            'amino_acid_percent': 'Distribution of percent amino acid composition of sequences'
        }
        self.dict_for_title_scatter = {
            'id': 'ID',
            'sequence_length': 'Sequence length',
            'molecular_weight': 'Molecular weight',
            'theoretical_pi': 'Theoretical isoelectric point',
            'aromaticity': 'Aromaticity',
            'instability_index': 'Instability index',
            'flexibility': 'Flexibility',
            'gravy': 'Average hydrophobicity and hydrophilicity',
            'c_term_freq': 'Frequency of C terminus',
            'n_term_freq': 'Frequency of N terminus',
            'secondary_structure_fraction': 'Secondary structure fraction',
            'extinction_coefficient': 'Extinction coefficient',
            'amino_acid_comp': 'Amino acid composition',
            'amino_acid_percent': 'Amino acid percent composition'
        }

    def get_colors(self, for_dict=None) -> Dict:
        """
        Function creates dictionary of color vectors based on length of list of sequences that will be plotted.
        If dictionary is passed, function will create dictionary of color vectors based on length of dictionary.

        :param for_dict: dict; if dictionary is passed, colors are based on it's length

        :return: dictionary with color vectors in format:
            {id1:[color_vector1], id2:[color_vector2]}
            or if dictionary passed:
            {dict_key1:[color_vector2], dict_key2:[color_vector2])

        :Example:

        >>> main_dict = {'id': ['BX5012', 'PH3210'], 'sequence_length': [21, 35], 'molecular_weight': [2399.853, 4048.622], 'theoretical_pi': [10.461, 9.9], 'instability_index': [42.448, 58.869], 'flexibility': [1.014, 1.043], 'gravy': [-0.752, -1.166], 'aromaticity': [0.048, 0.086], 'secondary_structure_fraction': [{'Helix': 23.81, 'Turn': 28.571, 'Sheet': 19.048},                              {'Helix': 22.857, 'Turn': 25.714, 'Sheet': 5.714}], 'extinction_coefficient': [{'Cysteines_reduced': 1, 'Cysteines_residues': 5},                        {'Cysteines_reduced': 2, 'Cysteines_residues': 3}], 'amino_acid_comp': [ {'A': 1, 'C': 0, 'D': 0, 'E': 1, 'F': 1, 'G': 1, 'H': 0, 'I': 2, 'K': 4, 'L': 1, 'M': 1, 'N': 1,  'P': 2, 'Q': 2, 'R': 1, 'S': 2, 'T': 0, 'V': 1, 'W': 0, 'Y': 0}, {'A': 0, 'C': 0, 'D': 5, 'E': 1, 'F': 3, 'G': 1, 'H': 0, 'I': 2, 'K': 9, 'L': 0, 'M': 1, 'N': 0,  'P': 1, 'Q': 0, 'R': 2, 'S': 7, 'T': 0, 'V': 3, 'W': 0, 'Y': 0}], 'amino_acid_percent': [ {'A': 4.762, 'C': 0.0, 'D': 0.0, 'E': 4.762, 'F': 4.762, 'G': 4.762, 'H': 0.0, 'I': 9.524,  'K': 19.048, 'L': 4.762, 'M': 4.762, 'N': 4.762, 'P': 9.524, 'Q': 9.524, 'R': 4.762,  'S': 9.524, 'T': 0.0, 'V': 4.762, 'W': 0.0, 'Y': 0.0}, {'A': 0.0, 'C': 0.0, 'D': 14.286, 'E': 2.857, 'F': 8.571, 'G': 2.857, 'H': 0.0, 'I': 5.714,  'K': 25.714, 'L': 0.0, 'M': 2.857, 'N': 0.0, 'P': 2.857, 'Q': 0.0, 'R': 5.714, 'S': 20.0,  'T': 0.0, 'V': 8.571, 'W': 0.0, 'Y': 0.0}], 'c_term_freq': {'M': 0.5, 'I': 0.5}, 'n_term_freq': {'R': 0.5, 'A': 0.5} }
        >>> small_dict = {'BX5012': { 'id': ['BX5012'], 'sequence_length': [43], 'molecular_weight': [4899.524], 'theoretical_pi': [6.516], 'instability_index': [35.179], 'flexibility': [0.996], 'gravy': [-0.749], 'aromaticity': [0.093], 'secondary_structure_fraction': [{'Helix': 18.605, 'Turn': 16.279, 'Sheet': 25.581}], 'extinction_coefficient': [{'Cysteines_reduced': 8480, 'Cysteines_residues': 8730}], 'amino_acid_comp': [{'A': 4, 'C': 4, 'D': 3, 'E': 3, 'F': 1, 'G': 3, 'H': 1, 'I': 0, 'K': 4, 'L': 2, 'M': 2, 'N': 2, 'P': 1, 'Q': 2, 'R': 2, 'S': 1, 'T': 3, 'V': 2, 'W': 1, 'Y': 2}], 'amino_acid_percent': [ {'A': 9.302, 'C': 9.302, 'D': 6.977, 'E': 6.977, 'F': 2.326, 'G': 6.977, 'H': 2.326, 'I': 0.0, 'K': 9.302, 'L': 4.651, 'M': 4.651, 'N': 4.651, 'P': 2.326, 'Q': 4.651, 'R': 4.651, 'S': 2.326, 'T': 6.977, 'V': 4.651, 'W': 2.326, 'Y': 4.651}], 'c_term_freq': {'G': 1}, 'n_term_freq': {'M': 1} }, 'PH3210': { 'id': ['PH3210'], 'sequence_length': [35], 'molecular_weight': [4048.622], 'theoretical_pi': [9.9], 'instability_index': [58.869], 'flexibility': [1.043], 'gravy': [-1.166], 'aromaticity': [0.086], 'secondary_structure_fraction': [{'Helix': 22.857, 'Turn': 25.714, 'Sheet': 5.714}], 'extinction_coefficient': [{'Cysteines_reduced': 0, 'Cysteines_residues': 0}], 'amino_acid_comp': [ {'A': 0, 'C': 0, 'D': 5, 'E': 1, 'F': 3, 'G': 1, 'H': 0, 'I': 2, 'K': 9, 'L': 0, 'M': 1, 'N': 0, 'P': 1, 'Q': 0, 'R': 2, 'S': 7, 'T': 0, 'V': 3, 'W': 0, 'Y': 0}], 'amino_acid_percent': [{'A': 0.0, 'C': 0.0, 'D': 14.286, 'E': 2.857, 'F': 8.571, 'G': 2.857, 'H': 0.0, 'I': 5.714, 'K': 25.714, 'L': 0.0, 'M': 2.857, 'N': 0.0, 'P': 2.857, 'Q': 0.0, 'R': 5.714, 'S': 20.0, 'T': 0.0, 'V': 8.571, 'W': 0.0, 'Y': 0.0}], 'c_term_freq': {'R': 1}, 'n_term_freq': {'M': 1} }}
        >>> v = SeqVisualizer(one_big_dict=main_dict, dict_of_many_small_dicts=small_dict)
        >>> dictionary_of_colours = v.get_colors()
        >>> print(dictionary_of_colours)
        {'BX5012': array([0.8  , 0.   , 0.128]), 'PH3210': array([0.8, 0. , 0.6])}

        """
        dict_of_colors = {}  # create dicitonary for colors
        if for_dict is None:
            list_for_colors = self.ids
        else:
            list_for_colors = for_dict
        colors = cm.gist_rainbow(np.linspace(0, 1, len(list_for_colors)))  # create array of colors
        colors = np.delete(colors, -1, axis=1)  # remove alpha column
        colors[:, 0:3] *= 0.8  # darken colors

        for seq_id, color in zip(list_for_colors, colors):  # fill final dictionary
            dict_of_colors[seq_id] = color

        return dict_of_colors

    def matplotlib_legend(self, plt, dict_for_leg, leg_size):
        """
        Usable with plots created with matplotlib library. If called, method adds a legend to the plot object.

        :param plt: plot object
        :param dict_for_leg: dict; dictionary of colors created by get.colors() function
        :param leg_size: int; size of legend

        :return: returned plot object with added legend
        """
        labels = dict_for_leg.keys()  # I don't know if I should apply here pycharms tip about a global variable
        handles = [plt.Rectangle((0, 0), 1, 1, color=dict_for_leg[label]) for label in labels]
        plt.legend(handles, labels, prop={'size': leg_size}, loc='center left',
                   bbox_to_anchor=(1, 0.5))  # adjusting legend size
        return plt

    def ax_and_sup_title_bar_box(self, plt, in_parameter, sup_fontsize, x_fontsize, y_fontsize, mode=None):
        """
        Method adds to the plot object title of plot and titles of axes. Works with two modes:
            1. Bar plots
            2. Box plots

        :param plt: plot object
        :param in_parameter: str; parameter chosen by user to plot
        :param sup_fontsize: int; size of font of plot title
        :param x_fontsize: int; size of font of x axis title
        :param y_fontsize: int; size of font of y axis title
        :param mode: 'bar' for bar plots or 'box' for box plots

        :return: plot object with all titles added
        """
        if mode == 'bar':
            plt.title(self.dict_for_title_bar[in_parameter], fontsize=sup_fontsize)
        elif mode == 'box':
            plt.title(self.dict_for_title_box[in_parameter], fontsize=sup_fontsize)
        else:
            print('Choose proper mode of ax_and_sup_title_bar_box function.')
        plt.xlabel(self.dict_for_x_ax_bar_box[in_parameter], fontsize=x_fontsize)
        plt.ylabel(self.dict_for_ax_units[in_parameter], fontsize=y_fontsize)
        return plt

    def ax_and_sup_title_scattter(self, plt, x_parameter, y_parameter, sup_fontsize, x_fontsize, y_fontsize,
                                  z_parameter=None):
        """
        Method adds to the plot object title of plot and titles of axes. Works with scatter plots

        :param plt: plot object
        :param x_parameter: str; parameter chosen by user to plot on x axis
        :param y_parameter: str; parameter chosen by user to plot on y axis
        :param sup_fontsize: str; size of font of plot title
        :param x_fontsize: int; size of font of x axis title
        :param y_fontsize: int; size of font of y axis title
        :param z_parameter: str; parameter chosen by user to plot as size of bubble

        :return: plot object with all titles added
        """
        if z_parameter is None or z_parameter is 100:
            suptitle = 'Comparision of ' + self.dict_for_title_scatter[x_parameter].lower() + 'and' + \
                       self.dict_for_title_scatter[y_parameter].lower() + ' of sequences'
        else:
            suptitle = 'Comparision of ' + self.dict_for_title_scatter[x_parameter].lower() + 'and' + \
                       self.dict_for_title_scatter[y_parameter].lower() + '(' + self.dict_for_title_scatter[
                           z_parameter].lower() + ' as bubble size) of seuqences'

        plt.title(suptitle, fontsize=sup_fontsize)
        plt.xlabel(self.dict_for_ax_units[x_parameter], fontsize=x_fontsize)
        plt.ylabel(self.dict_for_ax_units[x_parameter], fontsize=y_fontsize)
        return plt

    def show_or_save_plot(self, save, out_name):
        """
        Method used to show (by default) or save plots.

        :param save: if passed, plot will be saved in .png format instead of being shown
        :param out_name: name for plot (if not passed, default name will be applied)

        :return: plot showed interactively or saved as .png
        """

        if not save:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig(f'{out_name}.png')
        self.ids = self.main_dict['id']

    def show_or_save_plot_interactive(self, save, figure, out_name):
        """
        Method used to show (by default) or save plots.

        :param save: if passed, plot will be saved in .html format instead of being shown
        :param figure: figure object to save
        :param out_name: name for plot (if not passed, default name will be applied)

        :return: plot showed interactively or saved as .html
        """
        if not save:
            figure.show()
        else:
            figure.write_html(f'{out_name}.html')
        self.ids = self.main_dict['id']

    def bar(self, single_seq=None, parameter=None, save=False, size=(15, 8), x_ax_title_fontsize=14,
            y_ax_title_fontsize=14, suptitle_fontsize=13, add_legend=True, leg_size=13, out=None):
        """
        Method creating bar plot for chosen input parameter. Can work in single sequence mode.
            For parameters:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)
            each bar represents sequence

            For parameter:
            - Secondary structure fraction
            each bar represents secondary structure fraction

            For parameter:
            - Extinction coefficient
            each bar represents cystein groups

            For parameters:
            - Amino acid composition
            - Amino acid percent composition
            - Frequency of C terminus
            - Frequency of N terminus
            each bar represents amino acids

        :param single_seq: str; ID of single sequence from input dictionary; if passed, only this sequenced will be plotted
        :param parameter: str; parameter, which values will be represented
        :param save: boolean; if passed, plot will be saved (if not passed, plot is shown by default), by default False
        :param size: tuple; size of plot, passed as tuple eg. (12, 8)
        :param x_ax_title_fontsize: int; size of font of x axis title
        :param y_ax_title_fontsize: int; size of font of y axis title
        :param suptitle_fontsize: int; size of font of plot title
        :param add_legend: boolean; show legend, by default True
        :param leg_size: int; size of legend
        :param out: str; name of output file

        :return: plot shown or saved
        """
        if parameter is None or parameter not in self.all_args:
            print('Pass a parameter to plot properly.')

        #  Checking if single sequence is passed
        if single_seq in list(self.sub_dict.keys()):
            data = self.sub_dict[single_seq]
            self.ids = data['id']  # overwriting self.ids value, because of considering single-sequence mode
            data = data[parameter]
        else:
            data = self.main_dict[parameter]

        if parameter in self.lists_args or parameter in self.dicts_args:
            fig, ax = plt.subplots(figsize=size)

            if parameter in self.lists_args:
                dict_of_colors_for_legend = self.get_colors()
                ax.bar(x=self.ids, height=data, color=list(self.get_colors().values()))

            else:
                dict_of_colors_for_legend = self.get_colors(for_dict=data)
                ax.bar(x=list(data.keys()), height=list(data.values()), color=self.get_colors().values())

            self.ax_and_sup_title_bar_box(plt, parameter, sup_fontsize=suptitle_fontsize,
                                          x_fontsize=x_ax_title_fontsize,
                                          y_fontsize=y_ax_title_fontsize, mode='bar')

            #  Legend is turned on by default
            if add_legend:
                self.matplotlib_legend(plt, dict_of_colors_for_legend, leg_size)

        else:
            df = pd.DataFrame(data, index=self.ids).T  # Prepare dataframe directly for plotting
            df = df.rename_axis(self.dict_for_x_ax_bar_box[parameter]).reset_index()  # add 1st column with ids names
            df.plot(x=self.dict_for_x_ax_bar_box[parameter], y=self.ids, kind="bar", legend=add_legend,
                    rot=0)  # plotting

            self.ax_and_sup_title_bar_box(plt, parameter, sup_fontsize=suptitle_fontsize,
                                          x_fontsize=x_ax_title_fontsize,
                                          y_fontsize=y_ax_title_fontsize, mode='bar')

            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

        if out is None:
            out = f'bar_{parameter}'

        self.show_or_save_plot(save, out_name=out)

    def bar_interactive(self, single_seq=None, parameter=None, save=False, size=500, out=None):
        """
        Method creating interactive bar plot for chosen input parameter (for parameters check 'bar' documentation).
        Can work in single sequence mode.

        :param single_seq: str; ID of single sequence from input dictionary; if passed, only this sequenced will be plotted
        :param parameter: str; parameter, which values will be represented as height of bars
        :param save: boolean; if passed, plot will be saved (if not passed, plot is shown by default), by default False
        :param size: int; size of plot, by default 500
        :param out: str; name of output file

        :return: plot shown or saved
        """

        if parameter is None or parameter not in self.all_args:
            print('Pass a parameter to plot properly.')

        #  Checking if single sequence is passed
        if single_seq in list(self.sub_dict.keys()):
            data = self.sub_dict[single_seq]
            self.ids = data['id']  # overwriting self.ids value, because of considering single-sequence mode
            data = data[parameter]
        else:
            data = self.main_dict[parameter]

        if parameter in self.lists_args or parameter in self.dicts_args:
            if parameter in self.lists_args:
                df = pd.DataFrame({
                    self.dict_for_x_ax_bar_box[parameter]: self.ids,
                    self.dict_for_ax_units[parameter]: data
                }).reset_index()
            elif parameter in self.dicts_args:
                df = pd.DataFrame({
                    self.dict_for_x_ax_bar_box[parameter]: list(data.keys()),
                    self.dict_for_ax_units[parameter]: list(data.values())
                }).reset_index()
            fig = px.bar(df, x=self.dict_for_x_ax_bar_box[parameter], y=self.dict_for_ax_units[parameter],
                         color=self.dict_for_x_ax_bar_box[parameter], barmode='group', height=size,
                         title=self.dict_for_title_bar[parameter])

        else:
            df = pd.DataFrame(data, index=self.ids).T  # Prepare dataframe directly for plotting
            df = df.rename_axis(self.dict_for_x_ax_bar_box[parameter]).reset_index()  # add 1st column with ids names
            df = df.melt(id_vars=[self.dict_for_x_ax_bar_box[parameter]], value_vars=self.ids, var_name='Sequences',
                         value_name=self.dict_for_ax_units[parameter])

            fig = px.bar(df, x=self.dict_for_x_ax_bar_box[parameter], y=self.dict_for_ax_units[parameter],
                         color='Sequences', barmode='group', height=size, title=self.dict_for_title_bar[parameter])

        if out is None:
            out = f'bar_inter_{parameter}'

        self.show_or_save_plot_interactive(save, figure=fig, out_name=out)

    def box(self, parameter=None, save=False, size=(15, 8), x_ax_title_fontsize=14,
            y_ax_title_fontsize=14, suptitle_fontsize=13, add_legend=True, leg_size=13, out=None):
        """
        Method creating box plot for chosen input parameter.
            For parameters:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)
            box represents distribution of values for all sequences

            For parameter:
            - Secondary structure fraction
            each box represents distribution of values of secondary structure fraction for all sequences

            For parameter:
            - Extinction coefficient
            each box represents distribution of values of cystein groups for all sequences

            For parameters:
            - Amino acid composition
            - Amino acid percent composition
            each box represents distribution of values of amino acids for all sequences

        :param parameter: str; parameter, which values will be represented as height of bars
        :param save: boolean; if passed, plot will be saved (if not passed, plot is shown by default), by default False
        :param size: tuple; size of plot, passed as tuple eg. (12, 8)
        :param x_ax_title_fontsize: int; size of font of x axis title
        :param y_ax_title_fontsize: int; size of font of y axis title
        :param suptitle_fontsize: int; size of font of plot title
        :param add_legend: boolean; show legend, by default True
        :param leg_size: int; size of legend
        :param out: str; name of output file

        :return: plot shown or saved
        """

        if parameter is None or parameter not in self.box_args:
            print('Pass a parameter to plot properly.')

        fig, ax = plt.subplots(figsize=size)  # creating subclasses for plot

        if parameter in self.lists_args:
            ax.boxplot(x=self.main_dict[parameter])
            ax.set_xticklabels('')

        else:
            dictionary_of_lists = {key: [dictionary[key] for dictionary in self.main_dict[parameter]] for key in
                                   self.main_dict[parameter][0]}

            box_positions = np.arange(len(dictionary_of_lists.keys())) + 1
            ax.boxplot(x=list(dictionary_of_lists.values()), positions=box_positions)
            ax.set_xticklabels(list(dictionary_of_lists.keys()))

        self.ax_and_sup_title_bar_box(plt, parameter, sup_fontsize=suptitle_fontsize, x_fontsize=x_ax_title_fontsize,
                                      y_fontsize=y_ax_title_fontsize, mode='box')
        if out is None:
            out = f'box_{parameter}'

        self.show_or_save_plot(save, out_name=out)

    def scatter(self, single_seq=None, x_parameter=None, y_parameter=None, z_parameter=100, save=False, size=(15, 8),
                x_ax_title_fontsize=14, y_ax_title_fontsize=14, suptitle_fontsize=13, point_annotation=True,
                add_legend=True, leg_size=13, out=None):
        """
        Method used to create plot for chosen parameters (min. 2, max3). Can work in single sequence mode.
        Plotting values only for parameters:
            Possible parameters to plot:
            - Sequence length
            - Molecular weight
            - Theoretical isoelectric point
            - Aromaticity
            - Flexibility
            - Average hydrophobicity and hydrophilicity (GRAVY)


        :param single_seq: str; ID of single sequence from input dictionary; if passed, only this sequenced will be plotted
        :param x_parameter: str; parameter, which values will be represented on x axis
        :param y_parameter: str; parameter, which values will be represented on y axis
        :param z_parameter: str; parameter, which values will be represented as point size; optional
        :param save: boolean; if passed, plot will be saved (if not passed, plot is shown by default), by default False
        :param size: tuple; size of plot, passed as tuple eg. (12, 8)
        :param x_ax_title_fontsize: int; size of font of x axis title
        :param y_ax_title_fontsize: int; size of font of y axis title
        :param suptitle_fontsize: int; size of font of plot title
        :param point_annotation: boolean; showing annotation of points, by default True
        :param add_legend: boolean; show legend, by default True
        :param leg_size: int; size of legend
        :param out: str; name of output file

        :return: plot shown or saved
        """

        if x_parameter is None or \
                y_parameter is None or \
                x_parameter not in self.all_args or \
                y_parameter not in self.all_args:
            print('Pass a parameter to plot properly.')

        if x_parameter in self.lists_args and y_parameter in self.lists_args:
            fig, ax = plt.subplots(figsize=size)  # creating subclasses for plot

            if single_seq in list(self.sub_dict.keys()):
                data = self.sub_dict[single_seq]
                self.ids = data['id']  # overwriting self.ids value, because of considering single-sequence mode
            else:
                data = self.main_dict

            x_data = data[x_parameter]
            y_data = data[y_parameter]

            if z_parameter in self.lists_args:
                sc = ax.scatter(x=x_data, y=y_data, s=data[z_parameter], c=list(self.get_colors().values()))
            else:
                m_sizes = []
                for element in range(len(x_data)):
                    m_sizes.append(z_parameter)
                sc = ax.scatter(x=x_data, y=y_data, s=m_sizes, c=list(self.get_colors().values()))

            if point_annotation:
                for txt, size, xi, yi in zip(self.ids, sc.get_sizes(), x_data, y_data):
                    ax.annotate(txt, xy=(xi, yi), xytext=(np.sqrt(size) / 2 + 2, 0),
                                # xytext is used for proper positioning of labels
                                textcoords="offset points",
                                ha='left', va='center', )

            if add_legend:
                self.matplotlib_legend(plt, self.get_colors(), leg_size)

            self.ax_and_sup_title_scattter(plt, x_parameter=x_parameter, y_parameter=y_parameter,
                                           z_parameter=z_parameter, sup_fontsize=suptitle_fontsize,
                                           x_fontsize=x_ax_title_fontsize,
                                           y_fontsize=y_ax_title_fontsize)
            if out is None:
                out = f'scatter_{x_parameter}_{y_parameter}'

            self.show_or_save_plot(save, out_name=out)

        else:
            print('Parameter argument has been passed improperly.')


if __name__ == "__main__":
    import doctest

    doctest.testmod()
