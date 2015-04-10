# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Generic functions which are useful for working with HMMs.

This just collects general functions which you might like to use in
dealing with HMMs.
"""

from __future__ import print_function


def pretty_print_prediction(emissions, real_state, predicted_state,
                            emission_title="Emissions",
                            real_title="Real State",
                            predicted_title="Predicted State",
                            line_width=75):
    """Print out a state sequence prediction in a nice manner.

    Arguments:

    o emissions -- The sequence of emissions of the sequence you are
    dealing with.

    o real_state -- The actual state path that generated the emissions.

    o predicted_state -- A state path predicted by some kind of HMM model.
    """
    # calculate the length of the titles and sequences
    title_length = max(len(emission_title), len(real_title),
                       len(predicted_title)) + 1
    seq_length = line_width - title_length

    # set up the titles so they'll print right
    emission_title = emission_title.ljust(title_length)
    real_title = real_title.ljust(title_length)
    predicted_title = predicted_title.ljust(title_length)

    cur_position = 0
    # while we still have more than seq_length characters to print
    while True:
        if (cur_position + seq_length) < len(emissions):
            extension = seq_length
        else:
            extension = len(emissions) - cur_position

        print("%s%s" % (emission_title,
                        emissions[cur_position:cur_position + seq_length]))
        print("%s%s" % (real_title,
                        real_state[cur_position:cur_position + seq_length]))
        print("%s%s\n" % (predicted_title,
                          predicted_state[cur_position:
                                          cur_position + seq_length]))

        if (len(emissions) < (cur_position + seq_length)):
            break

        cur_position += seq_length
