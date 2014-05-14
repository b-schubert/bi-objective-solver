from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter
import matplotlib.pyplot as plt
import numpy as np

def plot_pareto_front(out, x, y, x_label="X Values", y_label="Y Values"):
    """
    plots the pareo front for biobjective problems given their objective values

    :param out: output file
    :param x: x values as list
    :param y: y values as list
    """
    fig = Figure(figsize=(6,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.set_title("Pareto Front", fontsize=14)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.grid(True,linestyle='-', color='0.75')
    #ax.invert_yaxis()
    #ax.invert_xaxis()
    ax.plot(x,y,"o-")
    canvas.print_figure(out, dpi=300)


def plot_variable_distribution(out, x, y, x_label="X Values", y_label="Y Values", title="Value Correlation"):
    """
        plots the score value distribution of two correlated variables.
        this function makes only sense for one pareto point at a time!

        :param out: output file
        :param x: x values as list -> these are going to be plotted on the x-axis
        :param y: y values as list -> these values are correlated with x and are plotted on the y-axis
        :invariant: x and y must be sorted!
    """
    #rcParams['xtick.minor.size'] = 40
    fig = Figure(figsize=(6,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111)
    ax.set_title(title, fontsize=14)
    ax.set_xlabel(x_label, fontsize=12)
    ax.set_ylabel(y_label, fontsize=12)
    ax.set_ylim((-0.01, 10))
    #ax.grid(False,linestyle='-', color='0.75')
    ax.set_xlim((-1,len(x)+1))
    ax.fill_between(range(len(x)), y)
    ax.set_xticks(range(len(x)))
    ax.set_xticklabels(x, fontsize=5)
    fig.tight_layout()
    canvas.print_figure(out, dpi=300)


def plot_3d_variable_distribution(out, x, ys, x_label="X Values", y_label="Y Values", title="Value Correlation" ):
    """
        plots a set of variable distributions together in a 3d plot

        :param out: output file
        :param x: x values as list of lists -> these are going to be plotted on the x-axis
        :param y: y values as list of lists -> these values are correlated with x and are plotted on the y-axis
        :invariant: x and y must be sorted!
    """
    colors = plt.cm.rainbow(np.linspace(0, 1, len(x)))
    fig = Figure(figsize=(6,6))
    canvas = FigureCanvas(fig)
    ax = fig.add_subplot(111, projection='3d')
    zs = range(len(x))
    xs = range(len(x[0]))
    for i, t in enumerate(zip(ys,colors)):
        y=t[0]
        c=t[1]
        ax.bar(xs, y, zs=i, zdir='y', color=c, alpha=0.6)
    #vert = [list(zip(xs, y)) for y in ys]
    #poly = PolyCollection(vert, facecolors=colors)
    #poly.set_alpha(0.7)
    #ax.add_collection3d(poly, zs=zs, zdir='y')
    ax.set_xlabel(x_label)
    ax.set_xlim3d(0, len(x[0]))
    ax.set_ylabel("Pareto Point")
    ax.set_ylim3d(-1, len(x))
    ax.set_zlabel(y_label)
    ax.set_zlim3d(0, 8)
    canvas.print_figure(out, dpi=300)