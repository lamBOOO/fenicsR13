# pylint: disable=invalid-name

"""
Module to store the postprocessor classes.
"""

from itertools import chain

import csv
from math import sqrt, ceil, log, log10
import numpy as np
import matplotlib.pyplot as plt


class Postprocessor:
    """
    The class for postprocessing.

    The constructor expects a list of data dictionaries of the form:

    .. code-block:: python

        [
            {
                'h': 0.9886573325052778,
                'theta': {'L_2': 0.18459230027019588,
                'l_inf': 0.11513941287387819},
                'sx': {'L_2': 0.38232086596749726,
                'l_inf': 0.30555017240112986},
                'sy': {'L_2': 0.51913814853123219,
                'l_inf': 0.2889203116681428}
            },
            {
                'h': 0.6340332990709842,
                'theta': {'L_2': 0.19952286856390053,
                'l_inf': 0.18948213107495288},
                'sx': {'L_2': 0.3528639452958619,
                'l_inf': 0.30118676712697767},
                'sy': {'L_2': 0.34778282823921497,
                'l_inf': 0.30792984640130788}
            }
        ]

    Parameters
    ----------
    data : list of dicts
        Data to plot, see above for shape.
    output_folder : string
        output folder for PDF plot.

    """

    def __init__(self, data, output_folder):
        """
        Construct postprocessor.

        Only sets arguments
        """
        self.data = data
        self.output_folder = output_folder

    def plot_errors(self, show_popup):
        """
        Use ``matplotlib`` to plot all errors in a figure.

        A slope marker is added.

        Exporting PDFs with:

        .. code-block:: python

            # Write PDF without access to backend, i.e. use PDF backend
            import matplotlib
            matplotlib.use('pdf')
            import matplotlib.pyplot as plt # pylint: disable=C0413
        """
        filename = "convergence_plot_" + self.output_folder + ".pdf"

        if not show_popup:
            plt.switch_backend('agg')
        plt.figure(num=None, figsize=(16, 9), dpi=100)

        data = self.data
        h = [df["h"] for df in data]

        raw_data = [key for key in data[0] if key != "h"]
        num_fields = len(raw_data)

        plt_sidelength_x = ceil(sqrt(num_fields))
        plt_sidelength_y = plt_sidelength_x
        while plt_sidelength_y * (plt_sidelength_x - 1) >= num_fields:
            plt_sidelength_x -= 1

        for i, key in enumerate(raw_data):
            field = [df[key] for df in data]
            plot_i = plt.subplot(plt_sidelength_x, plt_sidelength_y, (i + 1))
            for j, etype in enumerate(field[0]):
                values = [df[etype] for df in field]

                # Plot actual data
                plt.loglog(h, values, "-o", label=etype)

                # Add slope marker
                use_top = j == 1
                bot = np.array([
                    [h[-1], 0.5 * values[-1]],
                    [h[-2], 0.5 * values[-1]],
                    [h[-2], 0.5 * values[-2]]
                ])
                top = np.array([
                    [h[-1], 2.0 * values[-1]],
                    [h[-2], 2.0 * values[-2]],
                    [h[-1], 2.0 * values[-2]]
                ])
                tria = top if use_top else bot
                slope_marker = plt.Polygon(
                    tria, color=plot_i.get_lines()[-1].get_color(),
                    alpha=1.5, fill=False
                )
                plot_i.add_patch(slope_marker)

                delta_values = log(values[-2]) - log(values[-1])
                delta_h = log(h[-2]) - log(h[-1])
                conv_rate = delta_values / delta_h
                anchor_x = h[-1] if use_top else h[-2]
                anchor_y = (
                    10**(
                        (log10(2.0 * values[-2]) + log10(2.0 * values[-1])) / 2
                    ) if use_top
                    else 10**(
                        (log10(0.5 * values[-2]) + log10(0.5 * values[-1])
                         ) / 2)
                )
                h_align = "left" if use_top else "right"
                plot_i.text(
                    anchor_x, anchor_y, str(round(conv_rate, 2)),
                    alpha=1.0, ha=h_align, va="center"
                )

            # Add order slopes
            plt.loglog(h, np.array(2 * np.power(h, 1)), "--", label="O(h^1)")
            plt.loglog(h, np.array(0.02 * np.power(h, 2)), "--", label="O(h^2)")
            plt.loglog(h, np.array(0.02 * np.power(h, 3)), "--", label="O(h^3)")

            # Add information
            plt.xlabel("max(h)")
            plt.ylabel("Error")
            plt.title(key)
            plt.legend(loc='lower right')

        plt.tight_layout()

        output_path = self.output_folder + "/" + filename
        print("Write {}".format(output_path))
        plt.savefig(output_path, dpi=150)
        if show_popup:
            plt.show()

    def write_errors(self):
        r"""
        Write errors to a csv file.

        The output file (e.g. ``error.csv``) looks like:

        .. code-block:: text

            h,p_L_2,p_l_inf,ux_L_2,ux_l_inf,uy_L_2,uy_l_inf,sigmaxx_L_2,...
            0.9,0.2,0.2,0.1,0.1,0.2,0.1,0.2,0.5,0.2,0.3,0.2,0.5
            0.6,0.1,0.1,0.0,0.1,0.1,0.1,0.0,0.3,0.0,0.3,0.0,0.3
        """
        filename = "errors.csv"
        output_path = self.output_folder + "/" + filename
        data = self.data

        with open(output_path, mode='w') as file:
            print("Write {}".format(output_path))
            writer = csv.writer(file, delimiter=',', quotechar='"')

            # header
            writer.writerow(
                ["h"] + list(
                    chain.from_iterable([
                        [
                            first + "_" + second
                            for second in data[0].get(first, "")
                        ]
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
