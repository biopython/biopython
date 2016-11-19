# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.
#

"""Model a single layer in a nueral network.

These classes deal with a layers in the neural network (ie. the input layer,
hidden layers and the output layer).
"""
# standard library
import math
import random

from Bio._py3k import range


def logistic_function(value):
    """Transform the value with the logistic function.

    XXX This is in the wrong place -- I need to find a place to put it
    that makes sense.
    """
    return 1.0 / (1.0 + math.exp(-value))


class AbstractLayer(object):
    """Abstract base class for all layers.
    """
    def __init__(self, num_nodes, has_bias_node):
        """Initialize the layer.

        Arguments:

        o num_nodes -- The number of nodes that are contained in this layer.

        o has_bias_node -- Specify whether or not this node has a bias
        node. This node is not included in the number of nodes in the network,
        but is used in constructing and dealing with the network.
        """
        # specify all of the nodes in the network
        if has_bias_node:
            lower_range = 0
        else:
            lower_range = 1

        self.nodes = list(range(lower_range, num_nodes + 1))

        self.weights = {}

    def __str__(self):
        """Debugging output.
        """
        return "weights: %s" % self.weights

    def set_weight(self, this_node, next_node, value):
        """Set a weight value from one node to the next.

        If weights are not explicitly set, they will be initialized to
        random values to start with.
        """
        if (this_node, next_node) not in self.weights:
            raise ValueError("Invalid node values passed.")

        self.weights[(this_node, next_node)] = value


class InputLayer(AbstractLayer):
    def __init__(self, num_nodes, next_layer):
        """Initialize the input layer.

        Arguments:

        o num_nodes -- The number of nodes in the input layer.

        o next_layer -- The next layer in the neural network this is
        connected to.
        """
        AbstractLayer.__init__(self, num_nodes, 1)

        self._next_layer = next_layer

        # set up the weights
        self.weights = {}
        for own_node in self.nodes:
            for other_node in self._next_layer.nodes:
                self.weights[(own_node, other_node)] = \
                                        random.randrange(-2.0, 2.0)

        # set up the weight changes
        self.weight_changes = {}
        for own_node in self.nodes:
            for other_node in self._next_layer.nodes:
                self.weight_changes[(own_node, other_node)] = 0.0

        # set up the calculated values for each node -- these will
        # actually just be set from inputs into the network.
        self.values = {}
        for node in self.nodes:
            # set the bias node -- always has a value of 1
            if node == 0:
                self.values[0] = 1
            else:
                self.values[node] = 0

    def update(self, inputs):
        """Update the values of the nodes using given inputs.

        Arguments:

        o inputs -- A list of inputs into the network -- this must be
        equal to the number of nodes in the layer.
        """
        if len(inputs) != len(self.values) - 1:
            raise ValueError("Inputs do not match input layer nodes.")

        # set the node values from the inputs
        for input_num in range(len(inputs)):
            self.values[input_num + 1] = inputs[input_num]

        # propagate the update to the next layer
        self._next_layer.update(self)

    def backpropagate(self, outputs, learning_rate, momentum):
        """Recalculate all weights based on the last round of prediction.

        Arguments:

        o learning_rate -- The learning rate of the network

        o momentum - The amount of weight to place on the previous weight
        change.

        o outputs - The output info we are using to calculate error.
        """
        # first backpropagate to the next layers
        next_errors = self._next_layer.backpropagate(outputs, learning_rate,
                                                     momentum)

        for this_node in self.nodes:
            for next_node in self._next_layer.nodes:
                error_deriv = (next_errors[next_node] *
                               self.values[this_node])

                delta = (learning_rate * error_deriv +
                        momentum * self.weight_changes[(this_node, next_node)])

                # apply the change to the weight
                self.weights[(this_node, next_node)] += delta

                # remember the weight change for next time
                self.weight_changes[(this_node, next_node)] = delta


class HiddenLayer(AbstractLayer):
    def __init__(self, num_nodes, next_layer, activation=logistic_function):
        """Initialize a hidden layer.

        Arguments:

        o num_nodes -- The number of nodes in this hidden layer.

        o next_layer -- The next layer in the neural network that this
        is connected to.

        o activation -- The transformation function used to transform
        predicted values.
        """
        AbstractLayer.__init__(self, num_nodes, 1)

        self._next_layer = next_layer
        self._activation = activation

        # set up the weights
        self.weights = {}
        for own_node in self.nodes:
            for other_node in self._next_layer.nodes:
                self.weights[(own_node, other_node)] = \
                                        random.randrange(-2.0, 2.0)

        # set up the weight changes
        self.weight_changes = {}
        for own_node in self.nodes:
            for other_node in self._next_layer.nodes:
                self.weight_changes[(own_node, other_node)] = 0.0

        # set up the calculated values for each node
        self.values = {}
        for node in self.nodes:
            # bias node
            if node == 0:
                self.values[node] = 1
            else:
                self.values[node] = 0

    def update(self, previous_layer):
        """Update the values of nodes from the previous layer info.

        Arguments:

        o previous_layer -- The previous layer in the network.
        """
        # update each node in this network
        for update_node in self.nodes[1:]:
            # sum up the weighted inputs from the previous network
            sum = 0.0
            for node in previous_layer.nodes:
                sum += (previous_layer.values[node] *
                        previous_layer.weights[(node, update_node)])

            self.values[update_node] = self._activation(sum)

        # propagate the update to the next layer
        self._next_layer.update(self)

    def backpropagate(self, outputs, learning_rate, momentum):
        """Recalculate all weights based on the last round of prediction.

        Arguments:

        o learning_rate -- The learning rate of the network

        o momentum - The amount of weight to place on the previous weight
        change.

        o outputs - The output values we are using to see how good our
        network is at predicting things.
        """
        # first backpropagate to the next layers
        next_errors = self._next_layer.backpropagate(outputs, learning_rate,
                                                     momentum)

        # --- update the weights
        for this_node in self.nodes:
            for next_node in self._next_layer.nodes:
                error_deriv = (next_errors[next_node] *
                               self.values[this_node])

                delta = (learning_rate * error_deriv +
                        momentum * self.weight_changes[(this_node, next_node)])

                # apply the change to the weight
                self.weights[(this_node, next_node)] += delta

                # remember the weight change for next time
                self.weight_changes[(this_node, next_node)] = delta

        # --- calculate error terms
        errors = {}
        for error_node in self.nodes:
            # get the error info propagated from the next layer
            previous_error = 0.0
            for next_node in self._next_layer.nodes:
                previous_error += (next_errors[next_node] *
                                   self.weights[(error_node, next_node)])

            # get the correction factor
            corr_factor = (self.values[error_node] *
                           (1 - self.values[error_node]))

            # calculate the error
            errors[error_node] = previous_error * corr_factor

        return errors


class OutputLayer(AbstractLayer):
    def __init__(self, num_nodes, activation=logistic_function):
        """Initialize the Output Layer.

        Arguments:

        o num_nodes -- The number of nodes in this layer. This corresponds
        to the number of outputs in the neural network.

        o activation -- The transformation function used to transform
        predicted values.
        """
        AbstractLayer.__init__(self, num_nodes, 0)

        self._activation = activation

        self.values = {}
        for node in self.nodes:
            self.values[node] = 0

    def update(self, previous_layer):
        """Update the value of output nodes from the previous layers.

        Arguments:

        o previous_layer -- The hidden layer preceding this.
        """
        # update all of the nodes in this layer
        for update_node in self.nodes:
            # sum up the contribution from all of the previous inputs
            sum = 0.0
            for node in previous_layer.nodes:
                sum += (previous_layer.values[node] *
                        previous_layer.weights[(node, update_node)])

            self.values[update_node] = self._activation(sum)

    def backpropagate(self, outputs, learning_rate, momentum):
        """Calculate the backpropagation error at a given node.

        This calculates the error term using the formula:

        p = (z - t) z (1 - z)

        where z is the calculated value for the node, and t is the
        real value.

        Arguments:

        o outputs - The list of output values we use to calculate the
        errors in our predictions.
        """
        errors = {}
        for node in self.nodes:
            calculated_value = self.values[node]
            real_value = outputs[node - 1]

            errors[node] = ((real_value - calculated_value) *
                            calculated_value *
                            (1 - calculated_value))

        return errors

    def get_error(self, real_value, node_number):
        """Return the error value at a particular node.
        """
        predicted_value = self.values[node_number]
        return 0.5 * math.pow((real_value - predicted_value), 2)

    def set_weight(self, this_node, next_node, value):
        raise NotImplementedError("Can't set weights for the output layer")
