# Copyright 2001 by Jeffrey Chang.  All rights reserved.
#
# This file is part of the Biopython distribution and governed by your
# choice of the "Biopython License Agreement" or the "BSD 3-Clause License".
# Please see the LICENSE file that should have been included as part of this
# package.
"""Maximum Entropy code.

Uses Improved Iterative Scaling.
"""
# TODO Define terminology

from functools import reduce

try:
    import numpy
except ImportError:
    from Bio import MissingPythonDependencyError

    raise MissingPythonDependencyError(
        "Install NumPy if you want to use Bio.MaxEntropy."
    )


class MaxEntropy:
    """Hold information for a Maximum Entropy classifier.

    Members:
    classes      List of the possible classes of data.
    alphas       List of the weights for each feature.
    feature_fns  List of the feature functions.

    Car data from example Naive Bayes Classifier example by Eric Meisner November 22, 2003
    http://www.inf.u-szeged.hu/~ormandi/teaching

    >>> from Bio.MaxEntropy import train, classify
    >>> xcar = [
    ...     ['Red', 'Sports', 'Domestic'],
    ...     ['Red', 'Sports', 'Domestic'],
    ...     ['Red', 'Sports', 'Domestic'],
    ...     ['Yellow', 'Sports', 'Domestic'],
    ...     ['Yellow', 'Sports', 'Imported'],
    ...     ['Yellow', 'SUV', 'Imported'],
    ...     ['Yellow', 'SUV', 'Imported'],
    ...     ['Yellow', 'SUV', 'Domestic'],
    ...     ['Red', 'SUV', 'Imported'],
    ...     ['Red', 'Sports', 'Imported']]
    >>> ycar = ['Yes','No','Yes','No','Yes','No','Yes','No','No','Yes']

    Requires some rules or features

    >>> def udf1(ts, cl):
    ...     return ts[0] != 'Red'
    ...
    >>> def udf2(ts, cl):
    ...     return ts[1] != 'Sports'
    ...
    >>> def udf3(ts, cl):
    ...     return ts[2] != 'Domestic'
    ...
    >>> user_functions = [udf1, udf2, udf3]  # must be an iterable type
    >>> xe = train(xcar, ycar, user_functions)
    >>> for xv, yv in zip(xcar, ycar):
    ...     xc = classify(xe, xv)
    ...     print('Pred: %s gives %s y is %s' % (xv, xc, yv))
    ...
    Pred: ['Red', 'Sports', 'Domestic'] gives No y is Yes
    Pred: ['Red', 'Sports', 'Domestic'] gives No y is No
    Pred: ['Red', 'Sports', 'Domestic'] gives No y is Yes
    Pred: ['Yellow', 'Sports', 'Domestic'] gives No y is No
    Pred: ['Yellow', 'Sports', 'Imported'] gives No y is Yes
    Pred: ['Yellow', 'SUV', 'Imported'] gives No y is No
    Pred: ['Yellow', 'SUV', 'Imported'] gives No y is Yes
    Pred: ['Yellow', 'SUV', 'Domestic'] gives No y is No
    Pred: ['Red', 'SUV', 'Imported'] gives No y is No
    Pred: ['Red', 'Sports', 'Imported'] gives No y is Yes
    """

    def __init__(self):
        """Initialize the class."""
        self.classes = []
        self.alphas = []
        self.feature_fns = []


def calculate(me, observation):
    """Calculate the log of the probability for each class.

    me is a MaxEntropy object that has been trained.  observation is a vector
    representing the observed data.  The return value is a list of
    unnormalized log probabilities for each class.
    """
    scores = []
    assert len(me.feature_fns) == len(me.alphas)
    for klass in me.classes:
        lprob = 0.0
        for fn, alpha in zip(me.feature_fns, me.alphas):
            lprob += fn(observation, klass) * alpha
        scores.append(lprob)
    return scores


def classify(me, observation):
    """Classify an observation into a class."""
    scores = calculate(me, observation)
    max_score, klass = scores[0], me.classes[0]
    for i in range(1, len(scores)):
        if scores[i] > max_score:
            max_score, klass = scores[i], me.classes[i]
    return klass


def _eval_feature_fn(fn, xs, classes):
    """Evaluate a feature function on every instance of the training set and class (PRIVATE).

    fn is a callback function that takes two parameters: a
    training instance and a class.  Return a dictionary of (training
    set index, class index) -> non-zero value.  Values of 0 are not
    stored in the dictionary.
    """
    values = {}
    for i in range(len(xs)):
        for j in range(len(classes)):
            f = fn(xs[i], classes[j])
            if f != 0:
                values[(i, j)] = f
    return values


def _calc_empirical_expects(xs, ys, classes, features):
    """Calculate the expectation of each function from the data (PRIVATE).

    This is the constraint for the maximum entropy distribution. Return a
    list of expectations, parallel to the list of features.
    """
    # E[f_i] = SUM_x,y P(x, y) f(x, y)
    #        = 1/N f(x, y)
    class2index = {}
    for index, key in enumerate(classes):
        class2index[key] = index
    ys_i = [class2index[y] for y in ys]

    expect = []
    N = len(xs)
    for feature in features:
        s = 0
        for i in range(N):
            s += feature.get((i, ys_i[i]), 0)
        expect.append(float(s) / N)
    return expect


def _calc_model_expects(xs, classes, features, alphas):
    """Calculate the expectation of each feature from the model (PRIVATE).

    This is not used in maximum entropy training, but provides a good function
    for debugging.
    """
    # SUM_X P(x) SUM_Y P(Y|X) F(X, Y)
    # = 1/N SUM_X SUM_Y P(Y|X) F(X, Y)
    p_yx = _calc_p_class_given_x(xs, classes, features, alphas)

    expects = []
    for feature in features:
        sum = 0.0
        for (i, j), f in feature.items():
            sum += p_yx[i][j] * f
        expects.append(sum / len(xs))
    return expects


def _calc_p_class_given_x(xs, classes, features, alphas):
    """Calculate conditional probability P(y|x) (PRIVATE).

    y is the class and x is an instance from the training set.
    Return a XSxCLASSES matrix of probabilities.
    """
    prob_yx = numpy.zeros((len(xs), len(classes)))

    # Calculate log P(y, x).
    assert len(features) == len(alphas)
    for feature, alpha in zip(features, alphas):
        for (x, y), f in feature.items():
            prob_yx[x][y] += alpha * f
    # Take an exponent to get P(y, x)
    prob_yx = numpy.exp(prob_yx)
    # Divide out the probability over each class, so we get P(y|x).
    for i in range(len(xs)):
        z = sum(prob_yx[i])
        prob_yx[i] = prob_yx[i] / z
    return prob_yx


def _calc_f_sharp(N, nclasses, features):
    """Calculate a matrix of f sharp values (PRIVATE)."""
    # f#(x, y) = SUM_i feature(x, y)
    f_sharp = numpy.zeros((N, nclasses))
    for feature in features:
        for (i, j), f in feature.items():
            f_sharp[i][j] += f
    return f_sharp


def _iis_solve_delta(
    N, feature, f_sharp, empirical, prob_yx, max_newton_iterations, newton_converge
):
    """Solve delta using Newton's method (PRIVATE)."""
    # SUM_x P(x) * SUM_c P(c|x) f_i(x, c) e^[delta_i * f#(x, c)] = 0
    delta = 0.0
    iters = 0
    while iters < max_newton_iterations:  # iterate for Newton's method
        f_newton = df_newton = 0.0  # evaluate the function and derivative
        for (i, j), f in feature.items():
            prod = prob_yx[i][j] * f * numpy.exp(delta * f_sharp[i][j])
            f_newton += prod
            df_newton += prod * f_sharp[i][j]
        f_newton, df_newton = empirical - f_newton / N, -df_newton / N

        ratio = f_newton / df_newton
        delta -= ratio
        if numpy.fabs(ratio) < newton_converge:  # converged
            break
        iters = iters + 1
    else:
        raise RuntimeError("Newton's method did not converge")
    return delta


def _train_iis(
    xs,
    classes,
    features,
    f_sharp,
    alphas,
    e_empirical,
    max_newton_iterations,
    newton_converge,
):
    """Do one iteration of hill climbing to find better alphas (PRIVATE)."""
    # This is a good function to parallelize.

    # Pre-calculate P(y|x)
    p_yx = _calc_p_class_given_x(xs, classes, features, alphas)

    N = len(xs)
    newalphas = alphas[:]
    for i in range(len(alphas)):
        delta = _iis_solve_delta(
            N,
            features[i],
            f_sharp,
            e_empirical[i],
            p_yx,
            max_newton_iterations,
            newton_converge,
        )
        newalphas[i] += delta
    return newalphas


def train(
    training_set,
    results,
    feature_fns,
    update_fn=None,
    max_iis_iterations=10000,
    iis_converge=1.0e-5,
    max_newton_iterations=100,
    newton_converge=1.0e-10,
):
    """Train a maximum entropy classifier, returns MaxEntropy object.

    Train a maximum entropy classifier on a training set.
    training_set is a list of observations.  results is a list of the
    class assignments for each observation.  feature_fns is a list of
    the features.  These are callback functions that take an
    observation and class and return a 1 or 0.  update_fn is a
    callback function that is called at each training iteration.  It is
    passed a MaxEntropy object that encapsulates the current state of
    the training.

    The maximum number of iterations and the convergence criterion for IIS
    are given by max_iis_iterations and iis_converge, respectively, while
    max_newton_iterations and newton_converge are the maximum number
    of iterations and the convergence criterion for Newton's method.
    """
    if not training_set:
        raise ValueError("No data in the training set.")
    if len(training_set) != len(results):
        raise ValueError("training_set and results should be parallel lists.")

    # Rename variables for convenience.
    xs, ys = training_set, results

    # Get a list of all the classes that need to be trained.
    classes = sorted(set(results))

    # Cache values for all features.
    features = [_eval_feature_fn(fn, training_set, classes) for fn in feature_fns]
    # Cache values for f#.
    f_sharp = _calc_f_sharp(len(training_set), len(classes), features)

    # Pre-calculate the empirical expectations of the features.
    e_empirical = _calc_empirical_expects(xs, ys, classes, features)

    # Now train the alpha parameters to weigh each feature.
    alphas = [0.0] * len(features)
    iters = 0
    while iters < max_iis_iterations:
        nalphas = _train_iis(
            xs,
            classes,
            features,
            f_sharp,
            alphas,
            e_empirical,
            max_newton_iterations,
            newton_converge,
        )
        diff = [numpy.fabs(x - y) for x, y in zip(alphas, nalphas)]
        diff = reduce(numpy.add, diff, 0)
        alphas = nalphas

        me = MaxEntropy()
        me.alphas, me.classes, me.feature_fns = alphas, classes, feature_fns
        if update_fn is not None:
            update_fn(me)

        if diff < iis_converge:  # converged
            break
    else:
        raise RuntimeError("IIS did not converge")

    return me


if __name__ == "__main__":
    from Bio._utils import run_doctest

    run_doctest(verbose=0)
