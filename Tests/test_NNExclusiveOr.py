#!/usr/bin/env python
"""Test function to teach the neural network an XOR function.

This is a very basic test of Neural Network functionality.
"""
# Neural Network code we'll be using
from Bio.NeuralNetwork.Training import TrainingExample
from Bio.NeuralNetwork.BackPropagation import Layer
from Bio.NeuralNetwork.BackPropagation.Network import BasicNetwork

VERBOSE = 0


def main():
    """Train a neural network, and then test it to see how it does.

    Since we have so few examples, we use all of them for training,
    validation and testing.
    """
    print "Setting up training examples..."
    # set up the training examples
    examples = []
    examples.append(TrainingExample([0, 0], [0]))
    examples.append(TrainingExample([0, 1], [1]))
    examples.append(TrainingExample([1, 0], [1]))
    examples.append(TrainingExample([1, 1], [0]))

    # create the network
    output = Layer.OutputLayer(1)
    hidden = Layer.HiddenLayer(3, output)
    input = Layer.InputLayer(2, hidden)
    
    network = BasicNetwork(input, hidden, output)

    print "Training the network..."
    # train it
    learning_rate = .5
    momentum = .1
    network.train(examples, examples, stopping_criteria, learning_rate,
                  momentum)

    print "Predicting..."
    # try predicting
    for example in examples:
        prediction = network.predict(example.inputs)
        if VERBOSE:
            print "%s;%s=> %s" % (example.inputs, example.outputs, prediction)

def stopping_criteria(num_iterations, validation_error, training_error):
    """Define when to stop iterating.
    """
    if num_iterations % 100 == 0:
        if VERBOSE:
            print "error:", validation_error
    if num_iterations >= 2000:
        return 1

    return 0

main()
