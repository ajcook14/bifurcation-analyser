import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

def rectangles(ax, boxes, facecolor):

    rectangles = [Rectangle((i[2], i[0]), i[3] - i[2], i[1] - i[0], linewidth=2) for i in boxes]

    pc = PatchCollection(rectangles, facecolor=facecolor, alpha=0.5, edgecolor='r')

    ax.add_collection(pc)



if __name__ == '__main__':

    a_bounds = (0, 2)
    b_bounds = (-1, 1)

    fig, ax = plt.subplots(figsize=(16, 16))

    ax.set_xlim(b_bounds[0], b_bounds[1])
    ax.set_ylim(a_bounds[0], a_bounds[1])

    f = open('./regular1')
    
    regular_boxes1 = []

    while not (line := f.readline()) == '':

        regular_boxes1.append( list(map(float, line.strip().split(','))) )

    f.close()

    f = open('./regular2')
    
    regular_boxes2 = []

    while not (line := f.readline()) == '':

        regular_boxes2.append( list(map(float, line.strip().split(','))) )

    f.close()

    f = open('./special')
    
    special_boxes = []

    while not (line := f.readline()) == '':

        special_boxes.append( list(map(float, line.strip().split(','))) )

    f.close()

    rectangles(ax, regular_boxes1, 'y')
    rectangles(ax, regular_boxes2, 'b')
    rectangles(ax, special_boxes, 'g')

    """
    a = np.linspace(b_bounds[0], b_bounds[1], 50)
    x = np.linspace(a_bounds[0], a_bounds[1], 50)

    A, X = np.meshgrid(a, x)

    Q = A - X**2
    D = -2*X

    plt.contour(A, X, Q, levels=[0.], colors=['k'], alpha=1.0, linewidths=1)
    plt.contour(A, X, D, levels=[0.], colors=['b'], alpha=1.0, linewidths=1)
    """

    ax.set_xlabel('$b$')
    ax.set_ylabel('$a$')

    #plt.savefig('./tanh.png')

