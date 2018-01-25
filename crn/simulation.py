from crn import Expression

class Simulation:
    def __init__(self, sim):
        self.sim = sim

    def __getitem__(self, s):
        if type(s) is Expression:
            if not s.is_species():
                raise ValueError(
                        "can't get simulation time series of complex "
                        f"expression {item}")
            s, *_ = s.species.keys()

        return self.sim[s]

    def plot(self, filename=None, title=None):
        if filename:
            import matplotlib as mpl
            mpl.use('Agg')

        import matplotlib.pyplot as plt


        time = self.sim['time']
        for species, series in self.sim.items():
            if species not in ("time", "nothing"):
                plt.plot(time, series, label=f"[{species}]")

        plt.xlabel("time (seconds)")
        plt.ylabel("concentration (M)")
        plt.legend(loc="best")
        if title:
            plt.title(title)

        if filename:
            plt.savefig(filename)
        else:
            plt.show()
        plt.cla()
        plt.clf()

