# Copyright 2001 Brad Chapman.  All rights reserved.
#
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.


"""Provide classes for dealing with Training Neural Networks."""

# standard modules
import random


class TrainingExample(object):
    """Hold inputs and outputs of a training example."""

    def __init__(self, inputs, outputs, name=""):
        """Initialize the class."""
        self.name = name
        self.inputs = inputs
        self.outputs = outputs


class ExampleManager(object):
    """Manage a grouping of Training Examples.

    This is meant to make it easy to split a bunch of training examples
    into three types of data:

    - Training Data -- These are the data used to do the actual training
      of the network.
    - Validation Data -- These data are used to validate the network
      while training. They provide an independent method to evaluate how
      the network is doing, and make sure the network gets trained
      independent of noise in the training data set.
    - Testing Data -- The data which are used to verify how well a network
      works. They should not be used at all in the training process, so they
      provide a completely independent method of testing how well a network
      performs.

    """

    def __init__(self, training_percent=.4, validation_percent=.4):
        """Initialize the manager with the training examples.

        Arguments:
         - training_percent - The percentage of the training examples that
           should be used for training the network.
         - validation_percent - Percent of training examples for validating
           a network during training.

        Attributes:
         - train_examples - A randomly chosen set of examples for training
           purposes.
         - valdiation_examples - Randomly chosesn set of examples for
           use in validation of a network during training.
         - test_examples - Examples for training purposes.

        """
        if training_percent + validation_percent > 1.0:
            raise ValueError(
                'The sum of training and validation percentages > 100')

        self.train_examples = []
        self.validation_examples = []
        self.test_examples = []

        self.training_percent = training_percent
        self.validation_percent = validation_percent

    def add_examples(self, training_examples):
        """Add a set of training examples to the manager.

        Arguments:
         - training_examples - A list of TrainingExamples to manage.

        """
        placement_rand = random.Random()

        # assign exact example randomly to the example types
        for example in training_examples:
            chance_num = placement_rand.random()
            # assign with the specified percentage
            if chance_num <= self.training_percent:
                self.train_examples.append(example)
            elif chance_num <= (self.training_percent +
                                self.validation_percent):
                self.validation_examples.append(example)
            else:
                self.test_examples.append(example)
