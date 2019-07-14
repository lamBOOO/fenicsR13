# pylint: disable=invalid-name

"Class to post-process data"

from itertools import chain

import csv
import matplotlib.pyplot as plt
import numpy as np

class Postprocessor:
    "Postprocessor class"
    def __init__(self, data):
        self.data = data

    def plot_errors(self):
        "Plot the stored errors"
        data = self.data
        h = [df["h"] for df in data]
        for key in [key for key in data[0] if key != "h"]:
            field = [df[key] for df in data]
            for etype in field[0]:
                values = [df[etype] for df in field]
                plt.loglog(h, values, "-o", label=etype)
            plt.loglog(h, np.array(2*np.power(h, 1)), "--", label="O(h^1)")
            plt.loglog(h, np.array(0.02*np.power(h, 2)), "--", label="O(h^2)")
            plt.loglog(h, np.array(0.02*np.power(h, 3)), "--", label="O(h^3)")
            plt.xlabel("max(h)")
            plt.ylabel("Error")
            plt.title(key)
            plt.legend()
            plt.show()

    def write_errors(self):
        "Writes errors to csv file"
        filename = "errors.csv"
        data = self.data

        with open(filename, mode='w') as file:
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
