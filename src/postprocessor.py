# pylint: disable=invalid-name

"Class to post-process data"

from itertools import chain

import csv
from math import sqrt, ceil, log, log10
import numpy as np
import matplotlib.pyplot as plt

class Postprocessor:
    "Postprocessor class"
    def __init__(self, data, output_folder):
        self.data = data
        self.output_folder = output_folder

    def plot_errors(self, show_popup):
        """
        Use ``matplotlib`` to plot all errors in a figure.

        Exporting PDFs with:

        .. code-block:: python

            # Write PDF without access to backend, i.e. use PDF backend
            import matplotlib
            matplotlib.use('pdf')
            import matplotlib.pyplot as plt # pylint: disable=C0413
        """

        filename = "convergence_plot_" + self.output_folder  + ".pdf"

        if not show_popup:
            plt.switch_backend('agg')
        plt.figure(num=None, figsize=(16, 9), dpi=100)

        data = self.data
        h = [df["h"] for df in data]

        raw_data = [key for key in data[0] if key != "h"]
        num_fields = len(raw_data)

        plt_sidelength_x = ceil(sqrt(num_fields))
        plt_sidelength_y = plt_sidelength_x
        while plt_sidelength_y * (plt_sidelength_x-1) >= num_fields:
            plt_sidelength_x -= 1

        for i, key in enumerate(raw_data):
            field = [df[key] for df in data]
            plot_i = plt.subplot(plt_sidelength_x, plt_sidelength_y, (i+1))
            for j, etype in enumerate(field[0]):
                values = [df[etype] for df in field]
                plt.loglog(h, values, "-o", label=etype)

                use_top = j == 1
                bot = np.array([
                    [h[-1], 0.5*values[-1]],
                    [h[-2], 0.5*values[-1]],
                    [h[-2], 0.5*values[-2]]
                ])
                top = np.array([
                    [h[-1], 2.0*values[-1]],
                    [h[-2], 2.0*values[-2]],
                    [h[-1], 2.0*values[-2]]
                ])
                tria = top if use_top else bot
                slope_marker = plt.Polygon(
                    tria, color=plot_i.get_lines()[-1].get_color(),
                    alpha=1.5, fill=False
                    )
                plot_i.add_patch(slope_marker)

                conv_rate = (
                    (log(values[-2]) - log(values[-1]))
                    / (log(h[-2]) - log(h[-1]))
                )
                anchor_x = h[-1] if use_top else h[-2]
                anchor_y = (
                    10**(
                        (log10(2.0*values[-2])+log10(2.0*values[-1]))/2
                    ) if use_top
                    else 10**((log10(0.5*values[-2])+log10(0.5*values[-1]))/2)
                )
                h_align = "left" if use_top else "right"
                plot_i.text(
                    anchor_x, anchor_y, str(round(conv_rate, 2)),
                    alpha=1.0, ha=h_align, va="center"
                )

            plt.loglog(h, np.array(2*np.power(h, 1)), "--", label="O(h^1)")
            plt.loglog(h, np.array(0.02*np.power(h, 2)), "--", label="O(h^2)")
            plt.loglog(h, np.array(0.02*np.power(h, 3)), "--", label="O(h^3)")
            plt.xlabel("max(h)")
            plt.ylabel("Error")
            plt.title(key)
            plt.legend(loc='lower right')

        plt.tight_layout()
        plt.savefig(self.output_folder + "/" + filename, dpi=150)
        if show_popup:
            plt.show()

    def write_errors(self):
        "Writes errors to csv file"
        filename = "errors.csv"
        data = self.data

        with open(self.output_folder + "/" + filename, mode='w') as file:
            writer = csv.writer(file, delimiter=',', quotechar='"')

            # header
            writer.writerow(
                ["h"] + list(
                    chain.from_iterable([
                        [first+"_"+second for second in data[0].get(first, "")]
                        for first in [
                            first for first in data[0] if first != "h"
                        ]
                    ])
                )
            )

            # data of all runs
            for run in data:
                writer.writerow(
                    [run["h"]] + list(chain.from_iterable([
                        [run[field]["L_2"], run[field]["l_inf"]]
                        for field in [f for f in run if f != "h"]
                    ]))
                )

def plot_single(data_x_, data_y_, title_, legend_):
    "TODO"
    plt.loglog(data_x_, data_y_, "-o", label=legend_)
    plt.loglog(data_x_, np.array(
        2*np.power(data_x_, 1)), "--", label="1st order")
    plt.loglog(data_x_, np.array(
        0.02*np.power(data_x_, 2)), "--", label="2nd order")
    plt.xlabel("h_max")
    plt.ylabel("norm(theta_i-theta_{i-1})_L2")
    plt.title(title_)
    plt.legend()
    plt.show()
