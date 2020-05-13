#!/usr/bin/env python3
from typing import List

import matplotlib.pyplot as plt
import numpy as np


class SeqVisualizer:
    def __init__(self, df):
        self.df = df  # input dataframe
        self.ids = df['id'] # list of sequences ids
        self.nr_seqs = np.arange(len(self.ids))  # array out of ids, for X axis
        self.colors_per_seq = self.get_colors()  # list of colors depending on sequences amount
        self.lists_args = ['sequence_length', 'molecular_weight', 'theoretical_pl', 'extinction_coefficient', 'instability_index', 'gravy', 'est_half_life']
        self.dicts_args = ['c_term_freq', 'n_term_freq']
        self.list_of_dicts_args = ['amino_acid_comp', 'atomic_comp']
        self.dict_for_axes = {
            'id': 'ID',
            'sequence_length': 'Sequence length',
            'molecular_weight': 'Molecular weight',
            'theoretical_pl': 'Theoretical isoelectric point',
            'extinction_coefficient': 'Extinction coefficient',
            'instability_index': 'Instability index',
            'gravy': 'Average hydrophobicity and hydrophilicity',
            'est_half_life': 'Estimated half-life',
            'c_term_freq': 'Frequency of C terminus',
            'n_term_freq': 'Frequency of N terminus',
            'amino_acid_comp': 'Amino acid composition',
            'atomic_comp': 'Atomic composition'
        }

    def get_colors(self) -> List:
        """
        Function generates list of colors depending on amount of given sequences
        :return:
        """
        colors = plt.cm.get_cmap('hsv', len(self.ids))
        return colors

    def histogram(self, parameter=None, size=(15, 8)):
        """

        :param size:
        :param parameter:
        :return:
        """

        # if podana_sekwencja:
        #     df=slownik[sekwencja]

        if parameter is None:
            print('xd')

        elif parameter in self.lists_args:
            fig, ax = plt.subplots(figsize=size)  # creating subclasses for plot
            data = self.df[parameter]

            ax.bar(x=self.ids, height=data)
            ax.set_xlabel(self.dict_for_axes[parameter])
            plt.show()

        elif parameter in self.dicts_args:
            pass

        elif parameter in self.list_of_dicts_args:
            pass
        else:
            print('Parameter argument has been passed improperly.')

    def scatter(self, x_parameter=None, y_parameter=None, z_parameter=100, size=(12, 6), point_annotation=True):

        if x_parameter in self.lists_args and y_parameter in self.lists_args:
            fig, ax = plt.subplots(figsize=size)  # creating subclasses for plot
            x_data = self.df[x_parameter]
            y_data = self.df[y_parameter]

            if z_parameter in self.lists_args:
                sc = ax.scatter(x=x_data, y=y_data, s=self.df[z_parameter])
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

            ax.set_xlabel(self.dict_for_axes[x_parameter])  # setting x axis title
            ax.set_ylabel(self.dict_for_axes[y_parameter])  # setting y axis title
            plt.show()

        else:
            print('Parameter argument has been passed improperly.')
