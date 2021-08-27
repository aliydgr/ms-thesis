import math
import numpy as np
from scipy.integrate import odeint
from scipy.optimize import fsolve
from networkx.algorithms.distance_measures import diameter


STEADY_THRESHOLD = 0.0001
TIME_RANGE = 100


def library_parser(data, features):
    n = len(data)
    m = len(features)

    result = [None] * m
    for i in range(0, m):
        r = 1
        for elm in str(features[i]).split(' '):
            value = None
            if 'x' in elm:
                elm = elm.replace('x', '')
                if '^' in elm:
                    var_index_str = elm.split('^')[0]
                    power = int(elm.split('^')[1])
                    if '!' in var_index_str:
                        var_index = int(var_index_str.split('!')[1])
                        value = pow(1-data[var_index], power)
                    else:
                        var_index = int(var_index_str)
                        value = pow(data[var_index], power)
                else:
                    if '!' in elm:
                        var_index = int(elm.split('!')[1])
                        value = 1-data[var_index]
                    else:
                        var_index = int(elm)
                        value = data[var_index]
                r = r * value
        result[i] = r
    return result


def find_steady_state(func, x0):
    return fsolve(func, x0, args=(0,))


def sindy(data, t):
    coefficients = model.coefficients()
    features = model.feature_library.get_feature_names()
    cf = library_parser(data, features)
    result = np.empty(len(data))
    for index in range(0, len(data)):
        sum = 0
        for i in range(0, len(cf)):
            sum += coefficients[index][i] * cf[i]
        result[index] = sum
    return result


def generate_data_by_function(dynamic_func, x0, dt, t_range=TIME_RANGE):
    t_test = np.arange(0, t_range, dt)
    x_test = odeint(dynamic_func, x0, t_test)
    return x_test


def generate_data_by_model(dynamic_model, x0, dt):
    global model
    model = dynamic_model
    return generate_data_by_function(sindy, x0, dt)


def create_perturbed_state(alpha, perturbed_node, steady_state):
    number_of_nodes = len(steady_state)
    perturbation = np.zeros(number_of_nodes)
    perturbation[perturbed_node] = steady_state[perturbed_node] * alpha
    perturbed = np.add(steady_state, perturbation)
    return perturbed


def funcp(data, t, func, p, xp):
    result = func(data, t)
    result[p] = xp
    return result


def calculate_g_by_function(perturbation, steady_state, dynamic_func, dt):
    number_of_nodes = len(steady_state)
    g_matrix = np.empty((number_of_nodes,number_of_nodes))
    t_perturbed = np.linspace(0, dt, TIME_RANGE)

    flow = []
    for i in range(0, number_of_nodes):
        perturbation_i = np.zeros(number_of_nodes)
        perturbation_i[i] = perturbation[i]
        perturbed_i = np.add(steady_state, perturbation_i)

        x_perturbed_i = odeint(funcp, perturbed_i, t_perturbed, args=(dynamic_func, i, perturbed_i[i])).T
        flow.append(x_perturbed_i)
        final_state_i = [row[-1] for row in x_perturbed_i]
        diff_i = np.subtract(final_state_i, steady_state)

        dxi_xi = perturbed_i[i]/steady_state[i]
        for j in range(0, number_of_nodes):
            dxj_xj = diff_i[j]/steady_state[j]
            g_matrix[j,i] = abs(dxj_xj/dxi_xi)

    return g_matrix, flow


def calculate_g_by_model(perturbation, steady_state, dynamic_model, dt):
    global model
    model = dynamic_model
    return calculate_g_by_function(perturbation, steady_state, sindy, dt)


def calculate_f(g_matrix, degree):
    number_of_nodes = len(g_matrix)
    f = np.zeros(number_of_nodes)
    logf = np.zeros(number_of_nodes)
    for i in range(0, number_of_nodes):
        f_i = np.zeros(number_of_nodes)
        for n in range(0, number_of_nodes):
            sum_gmn = 0
            for m in range(0, number_of_nodes):
                sum_gmn += g_matrix[m,n]
            sum_gmi = 0
            for m in range(0, number_of_nodes):
                if i==m:
                    continue
                sum_gmi += g_matrix[m,i]/sum_gmn
            f_i[n] = g_matrix[i,n] * sum_gmi

        for n in range(0, number_of_nodes):
            if i==n:
                continue
            f[i] += f_i[n]
        f[i] /= (number_of_nodes-1)
        logf[i] = math.log(f[i], degree[i])

    return f, logf


def calculate_gamma(g_matrix, graph):
    def k_neighbors(graph, node, k):
        def aux(node_list):
            res = set()
            for n in node_list:
                for i in graph.neighbors(n):
                    res.add(i)
            return list(res)
        seen = [node]
        result = aux(seen)
        for _ in range(0, k-1):
            not_seen = [elm for elm in result if elm not in seen]
            result += [elm for elm in aux(not_seen) if elm not in seen]
            seen += not_seen
        return list(set(result))

    def gamma(g_matrix, graph, l):
        _sum = 0
        number_of_nodes = len(graph.degree)
        for j in range(0, number_of_nodes):
            for i in k_neighbors(graph, j, l):
                _sum += g_matrix[i][j]
        return _sum/number_of_nodes

    return [gamma(g_matrix, graph, l) for l in range(1,diameter(graph)+1)]
