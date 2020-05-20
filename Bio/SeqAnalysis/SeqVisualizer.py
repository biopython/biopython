#!/usr/bin/env python3
from typing import List, Dict

import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import pandas as pd


class SeqVisualizer:
    def __init__(self, one_big_dict, dict_of_many_small_dicts):
        self.main_dict = one_big_dict  # input dictionary
        self.sub_dict = dict_of_many_small_dicts
        self.ids = one_big_dict['id']  # list of sequences ids
        self.nr_seqs = np.arange(len(self.ids))  # array out of ids, for X axis
        self.lists_args = ['sequence_length', 'molecular_weight', 'theoretical_pi', 'aromaticity', 'instability_index',
                           'gravy']
        self.dicts_args = ['c_term_freq', 'n_term_freq']
        self.list_of_dicts_args = ['secondary_structure_fraction', 'extinction_coefficient', 'amino_acid_comp',
                                   'amino_acid_percent']
        self.dict_for_x_ax = {
            'id': 'ID',
            'sequence_length': 'Sequence length',
            'molecular_weight': 'Molecular weight',
            'theoretical_pi': 'Theoretical isoelectric point',
            'aromaticity': 'Aromaticity',
            'instability_index': 'Instability index',
            'gravy': 'Average hydrophobicity and hydrophilicity',
            'c_term_freq': 'Frequency of C terminus',
            'n_term_freq': 'Frequency of N terminus',
            'secondary_structure_fraction': 'Secondary structure fraction',
            'extinction_coefficient': 'Extinction coefficient',
            'amino_acid_comp': 'Amino acid composition',
            'amino_acid_percent': 'Amino acid percent composition'
        }
        self.dict_for_y_ax = {
            'id': 'ID',
            'sequence_length': 'Sequence length [Quantity]',
            'molecular_weight': 'Molecular weight [Da]',
            'theoretical_pi': 'Theoretical isoelectric point [pH]',
            'aromaticity': 'Aromaticity [Phe+Trp+Tyr/all aa]',
            'instability_index': 'Instability index',
            'gravy': 'Average hydrophobicity and hydrophilicity',
            'c_term_freq': 'Frequency of C terminus [amino acid]',
            'n_term_freq': 'Frequency of N terminus [amino acid]',
            'extinction_coefficient': 'Extinction coefficient [m^3/(mol*cm)]',
            'secondary_structure_fraction': 'Secondary structure fraction [%]',
            'amino_acid_comp': 'Amino acid composition [Quantity]',
            'amino_acid_percent': 'Amino acid percent composition [%]'
        }
        self.bar_titles = {
            'id': 'ID',
            'sequence_length': 'Comparision of lengths of sequences',
            'molecular_weight': 'Comparision of molecular weights of sequences',
            'theoretical_pi': 'Comparision of theoretical isolectric points of sequences',
            'aromaticity': 'Comparision of aromaticities of sequences',
            'instability_index': 'Comparision of instability indexes of sequences',
            'gravy': 'Comparision of hydrophobicity and hydrophilicity of sequences',
            'extinction_coeffcient': 'Comparision of extinction coefficients of sequences',
            'c_term_freq': 'Frequency of amino acids on C terminus',
            'n_term_freq': 'Frequency of amino acids on N terminus',
            'extinction_coefficient': 'Comparision of extinction coefficient values of sequences',
            'secondary_structure_fraction': 'Comparision of secondary structure values of sequences',
            'amino_acid_comp': 'Comparision of amino acid composition of sequences',
            'amino_acid_percent': 'Comparision of percent amino acid composition of sequences'
        }

    def get_colors(self, for_dict=None) -> Dict:
        """
        Function creates colormap based on number of sequences ploted (or number of keys in passed dictionary)

        :param for_dict: if dictionary an input, then colors will be generated based on number of keys in input
        :return:
        """
        dict_of_colors = {}  # create dicitonary for colors
        if for_dict is None:
            list_for_colors = self.ids
        else:
            list_for_colors = for_dict
        colors = cm.gist_rainbow(np.linspace(0, 1, len(list_for_colors)))  # create array of colors
        colors = np.delete(colors, -1, axis=1)  # remove alpha column
        colors[:, 0:3] *= 0.8  # darken colors

        for seq_id, color in zip(list_for_colors, colors):
            dict_of_colors[seq_id] = color

        return dict_of_colors

    def matplotlib_legend(self, plt, dict_for_leg, leg_size):
        labels = dict_for_leg.keys()  # I don't know if I should apply here pycharms tip about a global variable
        handles = [plt.Rectangle((0, 0), 1, 1, color=dict_for_leg[label]) for label in labels]
        plt.legend(handles, labels, prop={'size': leg_size}, loc='center left',
                   bbox_to_anchor=(1, 0.5))  # adjusting legend size
        return plt

    def ax_and_sup_title(self, plt, in_parameter, sup_fontsize, x_fontsize, y_fontsize):
        plt.title(self.bar_titles[in_parameter], fontsize=sup_fontsize)
        plt.xlabel(self.dict_for_x_ax[in_parameter], fontsize=x_fontsize)
        plt.ylabel(self.dict_for_y_ax[in_parameter], fontsize=y_fontsize)
        return plt

    def show_or_save_plot(self, save):
        if not save:
            plt.show()
        else:
            plt.tight_layout()
            plt.savefig('lol.png')
        self.ids = self.main_dict['id']

    def bar(self, single_seq=None, parameter=None, save=False, size=(15, 8), x_ax_title_fontsize=14,
            y_ax_title_fontsize=14, suptitle_fontsize=13, add_legend=True, leg_size=13):
        """

        :param single_seq:
        :param parameter:
        :param save:
        :param size:
        :param x_ax_title_fontsize:
        :param y_ax_title_fontsize:
        :param suptitle_fontsize:
        :param add_legend:
        :param leg_size:
        :return:
        """

        #  Checking if single sequence is passed
        if single_seq in list(self.sub_dict.keys()):
            data = self.sub_dict[single_seq]
            self.ids = data['id']  # overwriting self.ids value, because of considering single-sequence mode
            data = data[parameter]
        else:
            data = self.main_dict[parameter]

        if parameter is None:
            print('Pass a parameter to plot.')

        if parameter in self.lists_args or parameter in self.dicts_args:
            fig, ax = plt.subplots(figsize=size)

            if parameter in self.lists_args:
                dict_for_legend = self.get_colors()
                ax.bar(x=self.ids, height=data, color=self.get_colors().values())

            else:
                dict_for_legend = self.get_colors(for_dict=data)
                ax.bar(x=list(data.keys()), height=list(data.values()), color=self.get_colors().values())

            self.ax_and_sup_title(plt, parameter, sup_fontsize=suptitle_fontsize, x_fontsize=x_ax_title_fontsize,
                                  y_fontsize=y_ax_title_fontsize)

            #  Legend is turned on by default
            if add_legend:
                self.matplotlib_legend(plt, dict_for_legend, leg_size)

            self.show_or_save_plot(save)

        elif parameter in self.list_of_dicts_args:
            df = pd.DataFrame(data, index=self.ids).T  # Prepare dataframe directly for plotting
            df = df.rename_axis(parameter).reset_index()  # add 1st column with ids names
            df.plot(x=parameter, y=self.ids, kind="bar", legend=add_legend, rot=0)  # plotting

            self.ax_and_sup_title(plt, parameter, sup_fontsize=suptitle_fontsize, x_fontsize=x_ax_title_fontsize,
                                  y_fontsize=y_ax_title_fontsize)

            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))

            self.show_or_save_plot(save)

        else:
            print('Parameter argument has been passed improperly.')

    def scatter(self, x_parameter=None, y_parameter=None, z_parameter=100, size=(12, 6), point_annotation=True):

        if x_parameter in self.lists_args and y_parameter in self.lists_args:
            fig, ax = plt.subplots(figsize=size)  # creating subclasses for plot
            x_data = self.main_dict[x_parameter]
            y_data = self.main_dict[y_parameter]

            if z_parameter in self.lists_args:
                sc = ax.scatter(x=x_data, y=y_data, s=self.main_dict[z_parameter])
            else:
                m_sizes = []
                for element in range(len(x_data)):
                    m_sizes.append(z_parameter)
                sc = ax.scatter(x=x_data, y=y_data, s=m_sizes)

            if point_annotation:
                for txt, size, xi, yi in zip(self.ids, sc.get_sizes(), x_data, y_data):
                    ax.annotate(txt, xy=(xi, yi), xytext=(np.sqrt(size) / 2 + 2, 0),
                                textcoords="offset points",
                                ha='left', va='center', )

            ax.set_xlabel(self.dict_for_x_ax[x_parameter])  # setting x axis title
            ax.set_ylabel(self.dict_for_x_ax[y_parameter])  # setting y axis title
            plt.show()

        else:
            print('Parameter argument has been passed improperly.')
