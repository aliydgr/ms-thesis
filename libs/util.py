import matplotlib.pyplot as plt


def plot_flow(flow, mode='single', size=(16, 9), fnumber=1):
    plt.figure(fnumber)
    plt.rcParams["figure.figsize"] = size
    number_of_nodes = len(flow[0])
    if mode == 'single':
        for i in range(0,number_of_nodes):
            plt.plot(flow[:,i], label='x'+str(i))
        plt.legend()
    elif mode == 'multi':
        fpr = 2
        fig, axs = plt.subplots(int(number_of_nodes/fpr), fpr)
        for i in range(0,number_of_nodes):
            axsi = axs[int(i/fpr), i%fpr]
            axsi.plot(flow[:,i])
            axsi.set_title('x'+str(i))
    plt.show()


def find_matrix_dist(matrix):
    dist = []
    for row in matrix:
        for ij in row:
            dist.append(ij)
    return dist


def flow_prediction_error(predict, base, limit=None):
    number_of_nodes = len(base[0])
    limit = limit if limit else min(len(base), len(predict))
    sum_error = [0] * number_of_nodes
    for i in range(0,limit):
        for n in range(0, number_of_nodes):
            sum_error[n] += abs(predict[i][n] - base[i][n])/base[i][n]
    return [se/limit for se in sum_error]
