"""Represent Neural Networks.

This module contains classes to represent Generic Neural Networks that
can be trained.

Many of the ideas in this and other modules were taken from
Neil Schemenauer's bpnn.py, available from:

http://www.enme.ucalgary.ca/~nascheme/python/bpnn.py

My sincerest thanks to him for making this available for me to work from,
and my apologies for anything I mangled.
"""
# standard library
import math

class BasicNetwork(object):
    """Represent a Basic Neural Network with three layers.

    This deals with a Neural Network containing three layers:

    o Input Layer

    o Hidden Layer

    o Output Layer
    """
    def __init__(self, input_layer, hidden_layer, output_layer):
        """Initialize the network with the three layers.
        """
        self._input = input_layer
        self._hidden = hidden_layer
        self._output = output_layer

    def train(self, training_examples, validation_examples,
              stopping_criteria, learning_rate, momentum):
        """Train the neural network to recognize particular examples.

        Arguments:

        o training_examples -- A list of TrainingExample classes that will
        be used to train the network.

        o validation_examples -- A list of TrainingExample classes that
        are used to validate the network as it is trained. These examples
        are not used to train so the provide an independent method of
        checking how the training is doing. Normally, when the error
        from these examples starts to rise, then it's time to stop
        training.

        o stopping_criteria -- A function, that when passed the number of
        iterations, the training error, and the validation error, will
        determine when to stop learning.

        o learning_rate -- The learning rate of the neural network.

        o momentum -- The momentum of the NN, which describes how much
        of the prevoious weight change to use.
        """
        num_iterations = 0
        while 1:
            num_iterations += 1
            training_error = 0.0
            for example in training_examples:
                # update the predicted values for all of the nodes
                # based on the current weights and the inputs
                # This propogates over the entire network from the input.
                self._input.update(example.inputs)

                # calculate the error via back propogation
                self._input.backpropagate(example.outputs,
                                          learning_rate, momentum)
            
                # get the errors in our predictions
                for node in range(len(example.outputs)):
                    training_error += \
                             self._output.get_error(example.outputs[node],
                                                    node + 1)

            # get the current testing error for the validation examples
            validation_error = 0.0
            for example in validation_examples:
                predictions = self.predict(example.inputs)

                for prediction_num in range(len(predictions)):
                    real_value = example.outputs[prediction_num]
                    predicted_value = predictions[prediction_num]
                    validation_error += \
                            0.5 * math.pow((real_value - predicted_value), 2)

            # see if we have gone far enough to stop
            if stopping_criteria(num_iterations, training_error,
                                 validation_error):
                break

    def predict(self, inputs):
        """Predict outputs from the neural network with the given inputs.

        This uses the current neural network to predict outputs, no
        training of the neural network is done here.
        """
        # update the predicted values for these inputs
        self._input.update(inputs)

        output_keys = self._output.values.keys()
        output_keys.sort()

        outputs = []
        for output_key in output_keys:
            outputs.append(self._output.values[output_key])
        return outputs
