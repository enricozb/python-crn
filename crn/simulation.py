from crn import Expression

class Simulation:
    """
    A completed simulation of a CRN. Contains the time-series information
    of every species in the CRN throughout the simulation. Allows for quick
    plotting and extraction of simulation data.

    This class probably won't be constructed by a user, thus it's
    implementation is more internal.

    args:
        sim: Dict[str, np.ndarray]
            A dictionary of species name to concentration time series.
            This dictionary also contains other fields such as "time" and
            "nothing" which are used for plotting data. This dictionary
            makes no guarantees that these are the only things it will
            contain, it may have new fields that are added if they are
            needed.
    """
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
        """
        Plots the concentration of all of the species over time.

        args:
            filename: Optional[str]
                if present, save the plot to a file `filename`. Otherwise,
                the plot will show up as a new window.

            title: Optional[str]
                if present, the plot will have a title `title`.
        """
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

