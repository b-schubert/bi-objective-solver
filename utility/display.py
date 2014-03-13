from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure


def plot_pareto_front(out, x, y, x_label="Epitope Count", y_label="Energy Score"):
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

