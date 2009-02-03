#!/usr/bin/env python
"""General test for the NeuralNetwork libraries.

This exercises various elements of the BackPropagation NeuralNetwork
libraries.
"""
# standard library
import random
import unittest

# local stuff
from Bio.NeuralNetwork.Training import TrainingExample, ExampleManager
from Bio.NeuralNetwork.StopTraining import ValidationIncreaseStop


class StopTrainingTest(unittest.TestCase):
    """Test functionality for stopping training networks.
    """
    def test_validation_increase_stop(self):
        """Stop training when the ValidationExamples increase.
        """
        stopper = ValidationIncreaseStop(max_iterations = 20,
                                         min_iterations = 2)

        stopper.last_error = 1.0
        do_stop = stopper.stopping_criteria(5, 1.0, 1.5)
        assert do_stop == 1, \
               "Did not tell us to stop when validation error increased."

        stopper.last_error = 1.0
        do_stop = stopper.stopping_criteria(1, 1.0, 1.5)
        assert do_stop == 0, \
               "Told us to stop before we reached the minimum iterations."

        stopper.last_error = 1.0
        do_stop = stopper.stopping_criteria(25, 1.0, 0.5)
        assert do_stop == 1, \
               "Did not tell us to stop when reaching maximum iterations."


class ExampleManagerTest(unittest.TestCase):
    """Tests to make sure the example manager is working properly.
    """
    def setUp(self):
        self.num_examples = 500
        self.examples = []
        for make_example in range(self.num_examples):
            inputs = []
            for input_make in range(3):
                inputs.append(random.randrange(1, 7))
            outputs = [random.randrange(1, 7)]
            self.examples.append(TrainingExample(inputs, outputs))

    def test_adding_examples(self):
        """Make sure test examples are added properly.
        """
        manager = ExampleManager()

        # figure out the expected number of examples in each category
        expected_train = manager.training_percent * self.num_examples
        expected_validation = manager.validation_percent * self.num_examples
        expected_test = self.num_examples - expected_train \
                        - expected_validation

        manager.add_examples(self.examples)

        for expect, actual in [(expected_train, len(manager.train_examples)),
                               (expected_validation,
                                len(manager.validation_examples)),
                               (expected_test, len(manager.test_examples))]:
            
            wrong_percent = abs(expect - actual) / self.num_examples
            assert wrong_percent < .1, \
                   "Deviation in how examples were added, expect %s, got %s" \
                   % (expect, actual)
        
    def test_partioning_examples(self):
        """Test that we can change how to partition the test examples.
        """
        manager = ExampleManager(0, 0)
        manager.add_examples(self.examples)
        assert len(manager.test_examples) == self.num_examples, \
               "Did not partion correctly to test_examples."

        manager = ExampleManager(1.0, 0)
        manager.add_examples(self.examples)
        assert len(manager.train_examples) == self.num_examples, \
               "Did not partition correctly to train_examples."

        manager = ExampleManager(0, 1.0)
        manager.add_examples(self.examples)
        assert len(manager.validation_examples) == self.num_examples, \
               "Did not partition correctly to validation_examples."

if __name__ == "__main__":
    runner = unittest.TextTestRunner(verbosity = 2)
    unittest.main(testRunner=runner)
