# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Classes to help deal with stopping training a neural network.

One of the key issues with training a neural network is knowning when to
stop the training of the network. This is tricky since you want to keep
training until the neural network has 'learned' the data, but want to
stop before starting to learn the noise in the data.

This module contains classes and functions which are different ways to
know when to stop training. Remember that the neural network classifier
takes a function to call to know when to stop training, so the classes
in this module should be instaniated, and then the stop_training function
of the classes passed to the network.
"""

from __future__ import print_function


class ValidationIncreaseStop(object):
    """Class to stop training on a network when the validation error increases.

    Normally, during training of a network, the error will always decrease
    on the set of data used in the training. However, if an independent
    set of data is used for validation, the error will decrease to a point,
    and then start to increase. This increase normally occurs due to the
    fact that the network is starting to learn noise in the training data
    set. This stopping criterion function will stop when the validation
    error increases.
    """
    def __init__(self, max_iterations=None, min_iterations=0,
                 verbose=0):
        """Initialize the stopping criterion class.

        Arguments:

        o max_iterations - The maximum number of iterations that
        should be performed, regardless of error.

        o min_iterations - The minimum number of iterations to perform,
        to prevent premature stoppage of training.

        o verbose - Whether or not the error should be printed during
        training.
        """
        self.verbose = verbose
        self.max_iterations = max_iterations
        self.min_iterations = min_iterations

        self.last_error = None

    def stopping_criteria(self, num_iterations, training_error,
                          validation_error):
        """Define when to stop iterating.
        """
        if num_iterations % 10 == 0:
            if self.verbose:
                print("%s; Training Error:%s; Validation Error:%s"
                      % (num_iterations, training_error, validation_error))

        if num_iterations > self.min_iterations:
            if self.last_error is not None:
                if validation_error > self.last_error:
                    if self.verbose:
                        print("Validation Error increasing -- Stop")
                    return 1

        if self.max_iterations is not None:
            if num_iterations > self.max_iterations:
                if self.verbose:
                    print("Reached maximum number of iterations -- Stop")
                return 1

        self.last_error = validation_error
        return 0
