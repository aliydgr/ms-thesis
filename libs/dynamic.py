import math
import numpy as np
from scipy.integrate import odeint


STEADY_THRESHOLD = 0.0001
TIME_RANGE = 100


def library_parser(data, features):
    n = len(data)
    m = len(features)

    result = [None] * m
    for i in range(0, m):
        r = 1
        for elm in features[i].split(' '):
            var_index = None
            power = 1
            if 'x' in elm:
                elm = elm.replace('x', '')
                if '^' in elm:
                    var_index = int(elm.split('^')[0])
                    power = int(elm.split('^')[1])
                else:
                    var_index = int(elm)
                r = r * pow(data[var_index], power)
        result[i] = r
    return result


def find_steady_state(X):
    prev_state = X[0]
    for i in range(1,len(X)):
        state = X[i]
        diff = np.subtract(state, prev_state)
        prev_state = state
        if -STEADY_THRESHOLD < max(np.amax(diff), np.amin(diff), key=abs) < STEADY_THRESHOLD:
            return state


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


def generate_data(dynamic_model, x0_test, dt):
    global model
    model = dynamic_model
    t_test = np.arange(0, TIME_RANGE, dt)
    x_test = odeint(sindy, x0_test, t_test)
    return x_test


def create_perturbed_state(alpha, perturbed_node, steady_state):
    number_of_nodes = len(steady_state)
    perturbation = np.zeros(number_of_nodes)
    perturbation[perturbed_node] = (steady_state[perturbed_node] if steady_state[perturbed_node] > 0 else 1) * alpha
    perturbed = np.add(steady_state, perturbation)


def apply_perturbation(dynamic_model, perturbed, dynamic_func, dt):
    global model
    model = dynamic_model
    t_perturbed = np.arange(0, TIME_RANGE, dt)
    x_perturbed = odeint(dynamic_func, perturbed, t_perturbed)
    return x_perturbed


def flow(X, start=0, stop=float('inf')):
    prev_state = X[0]
    result = [prev_state]
    for i in range(1,len(X)):
        state = X[i]
        diff = np.subtract(state, prev_state)
        if start < i < stop:
            temp = [min(abs(s),100) for s in state]
            result = np.vstack((result, temp))
        prev_state = state
        if -STEADY_THRESHOLD < max(np.amax(diff), np.amin(diff), key=abs) < STEADY_THRESHOLD:
            return result
    print('chaos state')
    return result


def calculate_g(perturbation, steady_state, dynamic_func, dt):
    number_of_nodes = len(steady_state)
    g_matrix = np.empty((number_of_nodes,number_of_nodes))
    t_perturbed = np.arange(0, TIME_RANGE, dt)

    for i in range(0, number_of_nodes):
        perturbation_i = np.zeros(number_of_nodes)
        perturbation_i[i] = perturbation[i]
        perturbed_i = np.add(steady_state, perturbation_i[i])

        x_perturbed_i = odeint(dynamic_func, perturbed_i, t_perturbed)
        final_state_i = x_perturbed_i[-1]
        diff_i = np.subtract(final_state_i, steady_state)

        dxi_xi = diff_i[i]/steady_state[i]
        for j in range(0, number_of_nodes):
            dxj_xj = diff_i[j]/steady_state[j]
            g_matrix[i,j] = abs(dxi_xi/dxj_xj)

    return g_matrix


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
