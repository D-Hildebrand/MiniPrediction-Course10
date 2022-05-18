# https://medium.com/maheshkkumar/implementing-a-binary-classifier-in-python-b69d08d8da21

import numpy as np


class Perceptron(object):
    """ Perceptron Classifier

    Parameters
    ------------
    rate : float
        Learning rate (ranging from 0.0 to 1.0)
    number_of_iteration : int
        Number of iterations over the input dataset.

    Attributes:
    ------------

    weight_matrix : 1d-array
        Weights after fitting.

    error_matrix : list
        Number of misclassification in every epoch(one full training cycle on the training set)

    """

    def __init__(self, rate=0.01, number_of_iterations=100):
        self.rate = rate
        self.number_of_iterations = number_of_iterations

    def fit(self, X, y):
        """ Fit training data

        Parameters:
        ------------
        X : array-like, shape = [number_of_samples, number_of_features]
            Training vectors.
        y : array-like, shape = [number_of_samples]
            Target values.

        Returns
        ------------
        self : object

        """

        self.weight_matrix = np.zeros(1 + X.shape[1])
        self.errors_list = []

        for _ in range(self.number_of_iterations):
            errors = 0
            for xi, target in zip(X, y):
                update = self.rate * (target - self.predict(xi))
                self.weight_matrix[1:] += update * xi
                self.weight_matrix[0] += update
                errors += int(update != 0.0)
            self.errors_list.append(errors)
        return self

    def dot_product(self, X):
        """ Calculate the dot product """
        return (np.dot(X, self.weight_matrix[1:]) + self.weight_matrix[0])

    def predict(self, X):
        """ Predicting the label for the input data """
        return np.where(self.dot_product(X) >= 0.0, 1, 0)


if __name__ == '__main__':
    X = np.array(
        [[0, 0, 0], [0, 0, 1], [0, 1, 0], [0, 1, 1], [1, 0, 0], [1, 0, 1],
         [1, 1, 0]])
    y = np.array([0, 1, 1, 1, 1, 1, 1])
    p = Perceptron()
    p.fit(X, y)
    print(
        "Predicting the output of [1, 1, 1] = {}".format(p.predict([1, 1, 1])))
