import numpy as np


class Discretization:
    def __init__(self, DynamicalEquation, mode="fSSE"):
        self.Dynamics = DynamicalEquation
        self.DiscretizationMode = mode

print("a",np.array([
    [1,2],[2,3]
]))