'''
Author: CTC 2801320287@qq.com
Date: 2023-04-24 20:09:42
LastEditors: CTC 2801320287@qq.com
LastEditTime: 2023-05-24 00:16:23
Description: 

Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
'''
import numpy as np


class Discretization:
    def __init__(self, DynamicalEquation, mode="fSSE"):
        self.Dynamics = DynamicalEquation
        self.DiscretizationMode = mode

print("a")