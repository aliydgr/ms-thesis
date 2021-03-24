import networkx as nx
import numpy as np
import pysindy as ps
from .dynamic import *


def get_graph(number_of_nodes, type_):
    if type_ == 'erdos_renyi':
        return nx.erdos_renyi_graph(number_of_nodes, 0.2)
    if type_ == 'small_world':
        return nx.connected_watts_strogatz_graph(number_of_nodes, 4, 0.5)
    if type_ == 'barabashi':
        return nx.barabasi_albert_graph(number_of_nodes, int(number_of_nodes/5))
    if type_ == 'scale_free':
        return nx.scale_free_graph(number_of_nodes)
    print('unsupported type')


def _bio(data, t):
    F = 1
    B = 1
    R = 1

    n = len(data)
    result = np.empty(n)

    for i in range(0,n):
        sigma = 0
        for j in range (0,n):
            sigma = sigma + R*(adjacency[i,j])*data[i]*data[j]
        result[i] = F - B*data[i] - sigma

    return result


_library_functions_bio = [
    lambda : 1,
    lambda x : x,
    lambda x,y : x*y
]

_library_function_names_bio = [
    lambda : 1,
    lambda x : x,
    lambda x,y : '' + x + ' ' + y
]


def _population(data, t):
    a = 2
    b = 3
    beta = 1

    n = len(data)
    result = np.empty(n)

    for i in range(0,n):
        sigma = 0
        for j in range (0,n):
            sigma = sigma + (adjacency[i,j])*pow(data[j], a)
        result[i] = sigma - beta*pow(data[i], b)

    return result


_library_functions_population = [
    lambda : 1,
    lambda x : x,
    lambda x : x*x,
    lambda x : x*x*x
]
_library_function_names_population = [
    lambda : 1,
    lambda x : x,
    lambda x : x + '.' + x,
    lambda x : x + '.' + x + '.' + x
]


def get_dynamic_function(type_):
    if type_ == 'bio':
        return _bio
    if type_ == 'population':
        return _population
    print('unsupported type')


def get_custom_library(type_):
    if type_ == 'bio':
        return ps.CustomLibrary(library_functions=_library_functions_bio, function_names=_library_function_names_bio)
    if type_ == 'population':
        return ps.CustomLibrary(library_functions=_library_functions_population, function_names=_library_function_names_population)
    print('unsupported type')


def create_dataset(graph, dynamic_function, x0_train, dt):
    global adjacency
    adjacency = nx.to_numpy_matrix(graph).A
    x_train = generate_data_by_function(_bio, x0_train, dt)
    return x_train
