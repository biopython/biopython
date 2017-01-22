# -*- coding: utf-8 -*-

# Copyright 2017 by Joanna Zbijewska, Agata Gruszczyńska, Michał Karlicki.
# All rights reserved.
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#
# Note that:
# 1. RMSD scripts and license is available separately.
#    We added it in file: calculate_rmsd and README2.
#
# 2. RNAstructure license is available separately.
#    Please consult rna.urmc.rochester.edu .

"""
Calculate relaxed base pair score or simply base pair score between two structures
saved in bpseq format.
by: Joanna Zbijewska <asia.zbijewska@gmail.com>
On the basis of:
Agius P., Bennett KP., Zucker M.,"Comparing RNA secondary structures using a relaxed base-pair score";RNA,2010 May,16(5):865–878
"""

import re
import math

class BPSEQ():
    """This class serves to parse the bpseq file format of
    RNA secondary structure to get the base indices and base pairs"""

    def __init__(self, filename):
        self.filename = filename

    def bpseq_input(self):
        """Prepares the list of lines saved in bpseq file"""
        data = []
        with open("{}.bpseq".format(self.filename),"r") as sec_struct:
            sec_structure = sec_struct.readlines()
            for line in sec_structure:
                data.append(line)
        sec_struct.close()
        new_data=[]
        for element in data:
            new_data.append(element.strip('\n'))
        list_of_lists = []
        for element in new_data:
            test = '([0-9])\w+'
            if re.match(test,element):
                list_of_lists.append(element.split(' '))
        return(list_of_lists)

    def bpseq_indices(self):
        """Prepares base indices for further calculations"""
        list_of_lists = self.bpseq_input()
        base_indices = {}
        for i in range(len(list_of_lists)):
            alist = list_of_lists[i]
            base_indices[alist[0]] = alist[1]
        return(base_indices)

    def bpseq_pairs(self):
        """Prepares a dictionary of indices of base pairs"""
        list_of_lists = self.bpseq_input()
        pair_indices = {}
        for i in range(len(list_of_lists)):
            alist = list_of_lists[i]
            pair_indices[alist[0]] = alist[2]
        return(pair_indices)

class struct_comparison():
    """This class takes two filenames of structures in bpseq format
    to compare."""

    def __init__(self, filename1, filename2, t):
        """Initializes arguments for struct_comparison class
        filename1 - name of file with first structure to compare;
        filename2 - name of file with second structure to compare;
        t - value of the relaxation parameter, must be bigger than zero."""
        self.name1 = filename1
        self.name2 = filename2
        self.t = t

    def initialize(self):
        """This function prepares structures to compare using the BPSEQ class"""
        one = BPSEQ(self.name1)
        two = BPSEQ(self.name2)
        one = one.bpseq_pairs()
        two = two.bpseq_pairs()
        return([one,two])

    def bp_distances(self):
        """Calculates the distances between base pairs in both analyzed structures"""
        pairs_from_two = self.initialize()
        one_1 = []
        one_2 = []
        two_1 = []
        two_2 = []
        for key in pairs_from_two[0]:
            one_1.append(key)
        for value in pairs_from_two[0]:
            one_2.append(value)
        for key in pairs_from_two[1]:
            two_1.append(key)
        for value in pairs_from_two[1]:
            two_2.append(value)
        list_of_distances_1 = []
        list_of_distances_2 = []
        for n in range(len(one_1)):
            temp = []
            for i in range(len(two_1)):
                temp.append(abs(max((int(one_1[n])-int(two_1[i])),(int(one_2[n])-int(two_2[i])))))
            list_of_distances_1.append(temp)
        for n in range(len(two_1)):
            temp = []
            for i in range(len(one_1)):
                temp.append(abs(max((int(two_1[n])-int(one_1[i])),(int(two_2[n])-int(one_2[i])))))
            list_of_distances_2.append(temp)
        return([list_of_distances_1,list_of_distances_2])

    def bp_score(self):
        """Returns a list of distances between baise pairs (for nonzero distances) in two given structures"""
        list_bp_distances = self.bp_distances()
        score = []
        for part in list_bp_distances:
            for alist in part:
                score.append(min(alist))
        score = [x for x in score if x != 0]
        score.sort(reverse=True)
        return(score)

    def return_bp_score(self):
        """Returns the value of the base pair score - the old method of comparison
        for secondary RNA structures"""
        scores_list = self.bp_score()
        bp_score = len(scores_list)
        print('Your calculated base pare score is {}'.format(bp_score))

    def rbp_score(self):
        """After setting the relaxation parameter (t), this function calculates
        the relaxed base pair score"""
        t = self.t
        bp_score = self.bp_score()
        bp_value = len(bp_score)
        if bp_value == 0:
            print("Structures are identical!")
        else:
            t = float(t)
            if t < 0:
                print("Wrong number. Relaxation parameter must be bigger than zero.")
                self.rbp_score()
            elif t == 0:
                self.return_bp_score()
            else:
                if t > float(max(bp_score)):
                    rbp_score = 0
                    print('RBP score for relaxation parameter bigger than any base pair distance is {}.'.format(rbp_score))
                else:
                    possible_m_vals = []
                    for ind in range(len(bp_score)):
                        min_m = bp_score[ind]/t
                        min_m = math.floor(min_m)
                        if min_m < (ind+1):
                            possible_m_vals.append(min_m)
                        possible_m_vals = [m for m in possible_m_vals if m > 0]
                    if len(possible_m_vals)==0:
                        print('No RBP value for this relaxation parameter. Try again.')
                        self.rbp_score()
                    else:
                        rbp_score = min(possible_m_vals)
                        print('Your calculated relaxed base pare score for relaxation parameter t = {} is {}'.format(t,rbp_score))
                        return(rbp_score)

#a = struct_comparison('TMR_00273_structure_test','TMR_00200_structure_test', 1)
#a.rbp_score()
