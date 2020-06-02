#!/usr/bin/env python3
from typing import List, Dict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


class SeqVisualizer:
    def __init__(self, one_big_dict, dict_of_many_small_dicts):
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

        :param for_dict: if dictionary is passed, colors are based on it's length

        :return: dictionary with color vectors in format:
        {id1:[color_vector1], id2:[color_vector2]}
        or if dictionary passed:
        {dict_key1:[color_vector2], dict_key2:[color_vector2])

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
        labels = dict_for_leg.keys()  # I don't know if I should apply here pycharms tip about a global variable
        handles = [plt.Rectangle((0, 0), 1, 1, color=dict_for_leg[label]) for label in labels]
        plt.legend(handles, labels, prop={'size': leg_size}, loc='center left',
                   bbox_to_anchor=(1, 0.5))  # adjusting legend size
        return plt

    def ax_and_sup_title_bar_box(self, plt, in_parameter, sup_fontsize, x_fontsize, y_fontsize, mode=None):
        if mode == 'bar':
            plt.title(self.dict_for_title_bar[in_parameter], fontsize=sup_fontsize)
        elif mode == 'box':
            plt.title(self.dict_for_title_box[in_parameter], fontsize=sup_fontsize)
        else:
            print('Choose proper mode of ax_and_sup_title_bar_box funciton.')
        plt.xlabel(self.dict_for_x_ax_bar_box[in_parameter], fontsize=x_fontsize)
        plt.ylabel(self.dict_for_ax_units[in_parameter], fontsize=y_fontsize)
        return plt

    def ax_and_sup_title_box(self, plt, in_parameter, sup_fontsize, x_fontsize, y_fontsize):
        plt.title(self.dict_for_title_box[in_parameter], fontsize=sup_fontsize)
        plt.xlabel(self.dict_for_x_ax_bar_box[in_parameter], fontsize=x_fontsize)
        plt.ylabel(self.dict_for_ax_units[in_parameter], fontsize=y_fontsize)
        return plt

    def ax_and_sup_title_scattter(self, plt, x_parameter, y_parameter, sup_fontsize, x_fontsize, y_fontsize,
                                  z_parameter=None):
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

    def show_or_save_plot(self, save):
        if not save:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig('lol.png')
        self.ids = self.main_dict['id']

    def show_or_save_plot_interactive(self, save, figure):
        if not save:
            figure.show()
        else:
            figure.write_html('lol.png')
        self.ids = self.main_dict['id']

    def bar(self, single_seq=None, parameter=None, save=False, size=(15, 8), x_ax_title_fontsize=14,
            y_ax_title_fontsize=14, suptitle_fontsize=13, add_legend=True, leg_size=13):
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
                dict_for_legend = self.get_colors()
                ax.bar(x=self.ids, height=data, color=list(self.get_colors().values()))

            else:
                dict_for_legend = self.get_colors(for_dict=data)
                ax.bar(x=list(data.keys()), height=list(data.values()), color=self.get_colors().values())

            self.ax_and_sup_title_bar_box(plt, parameter, sup_fontsize=suptitle_fontsize,
                                          x_fontsize=x_ax_title_fontsize,
                                          y_fontsize=y_ax_title_fontsize, mode='bar')

            #  Legend is turned on by default
            if add_legend:
                self.matplotlib_legend(plt, dict_for_legend, leg_size)

            self.show_or_save_plot(save)

        elif parameter in self.list_of_dicts_args:
            df = pd.DataFrame(data, index=self.ids).T  # Prepare dataframe directly for plotting
            df = df.rename_axis(self.dict_for_x_ax_bar_box[parameter]).reset_index()  # add 1st column with ids names
            df.plot(x=self.dict_for_x_ax_bar_box[parameter], y=self.ids, kind="bar", legend=add_legend,
                    rot=0)  # plotting

            self.ax_and_sup_title_bar_box(plt, parameter, sup_fontsize=suptitle_fontsize,
                                          x_fontsize=x_ax_title_fontsize,
                                          y_fontsize=y_ax_title_fontsize, mode='bar')

            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            self.show_or_save_plot(save)

    def bar_interactive(self, single_seq=None, parameter=None, save=False, size=500):

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

            self.show_or_save_plot_interactive(save, figure=fig)

        elif parameter in self.list_of_dicts_args:
            df = pd.DataFrame(data, index=self.ids).T  # Prepare dataframe directly for plotting
            df = df.rename_axis(self.dict_for_x_ax_bar_box[parameter]).reset_index()  # add 1st column with ids names
            df = df.melt(id_vars=[self.dict_for_x_ax_bar_box[parameter]], value_vars=self.ids, var_name='Sequences',
                         value_name=self.dict_for_ax_units[parameter])

            fig = px.bar(df, x=self.dict_for_x_ax_bar_box[parameter], y=self.dict_for_ax_units[parameter],
                         color='Sequences', barmode='group', height=size, title=self.dict_for_title_bar[parameter])

            self.show_or_save_plot_interactive(save, figure=fig)

    def box(self, parameter=None, save=False, size=(15, 8), x_ax_title_fontsize=14,
            y_ax_title_fontsize=14, suptitle_fontsize=13, add_legend=True, leg_size=13):

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

        self.show_or_save_plot(save)

    def scatter(self, single_seq=None, x_parameter=None, y_parameter=None, z_parameter=100, save=False, size=(15, 8),
                x_ax_title_fontsize=14, y_ax_title_fontsize=14, suptitle_fontsize=13, point_annotation=True,
                add_legend=True, leg_size=13):

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
                                textcoords="offset points",
                                ha='left', va='center', )

            if add_legend:
                self.matplotlib_legend(plt, self.get_colors(), leg_size)

            self.ax_and_sup_title_scattter(plt, x_parameter=x_parameter, y_parameter=y_parameter,
                                           z_parameter=z_parameter, sup_fontsize=suptitle_fontsize,
                                           x_fontsize=x_ax_title_fontsize,
                                           y_fontsize=y_ax_title_fontsize)

            self.show_or_save_plot(save)

        else:
            print('Parameter argument has been passed improperly.')


if __name__ == "__main__":
    import doctest

    doctest.testmod()
