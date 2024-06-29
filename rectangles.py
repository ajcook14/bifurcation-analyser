import numpy as np

import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

import os

def rectangles(ax, boxes, facecolor):

    rectangles = [Rectangle((i[2], i[0]), i[3] - i[2], i[1] - i[0], linewidth=2) for i in boxes]

    pc = PatchCollection(rectangles, facecolor=facecolor, alpha=0.5, edgecolor='k')

    ax.add_collection(pc)



if __name__ == '__main__':

    a_bounds = (0., 2.)
    b_bounds = (-1., 1.)

    fig, ax = plt.subplots(figsize=(20, 20))

    ax.set_xlim(b_bounds[0], b_bounds[1])
    ax.set_ylim(a_bounds[0], a_bounds[1])

    rng = np.random.default_rng(2)

    cnames = list(os.listdir('./components/'))

    cnames.sort()

    regular_colours = np.linspace(0.1, 1.0, len(cnames))

    for i, component_name in enumerate(cnames):

        f = open(f'./components/{component_name}')
        
        boxes = []

        while not (line := f.readline()) == '':

            boxes.append( list(map(float, line.strip().split(','))) )

        f.close()

        if component_name == "verified":

            rectangles( ax, boxes, 'r' )

        elif component_name == "special":

            rectangles( ax, boxes, 'y' )

        else:

            rectangles( ax, boxes, np.array([0., 0., regular_colours[i]]) )



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

    plt.show()

