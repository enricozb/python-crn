import numpy as np

from crn import Species, Simulation, utils
from random import random
from scipy.integrate import odeint

with utils.no_output():
    import stochpy

class CRN:
    """ A Chemical Reaction Network (CRN)

    args:
        *system: List[Reaction]
            The chemical reactions that define the system.
    attributes:
        system: List[Reaction]
            The chemical reactions that define the system.
        species: Set[Species]
            The set of species string names that are present in the CRN.
        kwargs: Dict[]
            name: str
                defaults to `id(self)`. Used for naming any files produced
                by the CRN during simulation.
    """

    def __init__(self, *system, **kwargs):
        self.system = system
        self.species = self.get_species()
        self.species_index = self.get_species_index()
        self.reactions_index = self.get_reactions_index()
        self.diffeq_system_func = self.rate_laws()
        self.name = kwargs.get("name", id(self))

    def get_species(self):
        """
        Returns the set of species present in the CRN.
        """
        species = set()
        for s in self.system:
            species |= s.get_species()
        return species

    def get_species_index(self):
        """
        Get a map of int to species, which is unchanged for the lifetime
        of the CRN.

        This is meant for internal use, and won't be very useful for anyone
        using the CRN.
        """
        return dict(enumerate(sorted(self.species)))

    def get_reactions_index(self):
        """
        Get a map of int to reactions, which is unchanged for the lifetime
        of the CRN.

        This is meant for internal use, and won't be very useful for anyone
        using the CRN.
        """
        return dict(enumerate(self.system))

    def rate_law_for_species(self, s):
        """
        Returns the symbolic representation for the rate law of species `s`.

        args:
            s: Union[Species, str]
                the species whose rate law in the CRN `self` is desired.
                `s` can either be an instance of the Species object or
                the string name of the Species.
        """
        if type(s) not in (str, Species):
            raise ValueError("rate_law_for_species called on non-species. "
                             "parameter must be a species name (str) or a "
                             "Species instance.")

        if type(s) is Species:
            s = species.name

        return sum(rxn.net_production(s) * rxn.flux() for rxn in self.system)

    def rate_laws(self):
        """
        Returns a function that takes a list of species concentrations,
        in the same order as specified in `self.species_index`, and returns
        the a vector of the current rate of change of each species, again in
        the same order as specified in `self.species_index`.
        """
        laws = []

        # Has to be iterated this way to laws is in order
        for i in range(len(self.species_index)):
            sp = self.species_index[i]
            law = self.rate_law_for_species(sp)
            for i, s in self.species_index.items():
                sub = 1 if s == "nothing" else f"v[{i}]"
                law = law.subs(s, sub)
            laws.append(str(law))

        func = f"lambda v, t: [{', '.join(laws)}]"
        return eval(func)

    def stoch_simulate(self, amounts, t=20):
        """
        Stochastic discrete simulation of the CRN until time `t` with initial
        molecule count `amounts`. The species that are omitted from the
        dictionary of initial concentrations are assumed to have an initial
        count of 0. The simulation runs until time `t` unless the
        propensities all reach zero. If that's the case, the system has
        reached a steady state, and the simulation stops.

        args:
            amount: Dict[Species, int]
                A map describing each species' initial count.
            t: Union[float, int]
                The upper bound of the time to run the simulation to.
        """

        # Write psc file
        pscfile = utils.datadir(f"{self.name}.psc")
        self.write_pscfile(pscfile, amounts)

        # invoke stochpy
        smod = stochpy.SSA()
        smod.Model(pscfile)
        smod.DoStochSim(mode="time", end=t)

        # return data
        data = {}
        for species in self.species:
            if species != "nothing":
                data[species] = smod.data_stochsim.getSimData(species)[:, 1]
                if "time" not in data:
                    data["time"] = smod.data_stochsim.getSimData(species)[:, 0]

        return Simulation(data, stochastic=True)

    def write_pscfile(self, filename, amounts):
        """
        Write the CRN in PySCeS Model Description Language for stochastic
        discrete simulation to `filename`.

        args:
            filename: str
                name of the .psc file to be written
            amounts: Dict[Species, int]
                A map describing each species' initial count.
        """

        def fmt_expression(expr):
            formatted_expr = []
            for species, c in expr.species.items():
                # special PySCeS Model Description Language Syntax
                if species == "nothing":
                    species = "$pool"

                if c != 1:
                    formatted_expr.append(f"{{{c}}}{species}")
                else:
                    formatted_expr.append(species)

            return " + ".join(formatted_expr)

        with open(filename, "w") as pscfile:
            pscfile.write("# Reactions\n\n")
            for i, rxn in self.reactions_index.items():
                pscfile.write(f"R{i}:\n")
                reac = fmt_expression(rxn.reactants)
                prod = fmt_expression(rxn.products)
                pscfile.write(f"{reac} > {prod}\n")
                pscfile.write(f"k{i}*{rxn.discrete_flux()}\n\n")
            pscfile.write("\n")

            pscfile.write("# Rate constants\n")
            for i, rxn in self.reactions_index.items():
                pscfile.write(f"k{i} = {rxn.coeff}\n")

            pscfile.write("\n# Initial Species Counts\n")

            # The `amounts` dictionary can have keys of type Species, so
            # we convert everything to the string names of the species.
            amounts_copy = {}
            for sp, amount in amounts.items():
                if type(sp) is Species:
                    sp = sp.name
                amounts_copy[sp] = amount

            for sp in self.species:
                if sp not in amounts_copy:
                    amounts_copy[sp] = 0

            if "nothing" in amounts_copy:
                del amounts_copy["nothing"]

            amounts = amounts_copy

            for sp in self.species:
                if sp != "nothing":
                    pscfile.write(f"{sp} = {amounts.get(sp, 0)}\n")


    def simulate(self, conc, t=20, resolution=100):
        """
        Deterministic concentration-continuous simulation of the CRN until
        time t with initial concentrations `conc`.
        The species that are omitted from the dictionary of initial
        concentrations are assumed to have an initial concentration of 0.0.

        args:
            conc: Dict[Species, float]
                A map describing each species' initial concentration.
            t: Union[float, int]
                The upper bound of the time to run the simulation to.
            resolution: int
                How many time steps to simulate between times [0, t).
        """

        t = np.linspace(0, t, resolution)

        conc_temp = {}

        for s, c in conc.items():
            if type(s) not in (str, Species):
                raise ValueError("CRN.simulate: initial concentrations "
                                 f"dictionary got a key with type {type(s)}. "
                                 "Key type should be Species or str.")

            if type(s) is Species:
                s = s.name

            conc_temp[s] = c

        conc = conc_temp

        v0 = [0] * len(self.species)

        for i, s in self.species_index.items():
            v0[i] = conc.get(s, 0)

        sol = odeint(self.diffeq_system_func, v0, t)

        sol_dict = {'time': t}

        for i, s in self.species_index.items():
            sol_dict[s] = sol[:, i]

        return Simulation(sol_dict)

    def validate(self, func, *, input_species, output_species, N=100,
            eps=1e-2, t=500):
        """
        Determine if the CRN `self` actually describes the computation in
        function `func`. This is a probabalistic verification: it runs
        several simulations at several different initial concentrations,
        then computes the desired function output and checks if the CRNs
        simulation is within `eps` of `func`s output.

        args:
            func: Callable[Dict[Expression, float], float]
                This function takes in the initial concentrations of the
                input species to the CRN and outputs a number computed from
                these initial concentrations.
        """
        for i in range(N):
            species = {sp : random() * 10 for sp in input_species}
            theoretical = func(species)
            sim = self.simulate(species, t=t)
            simulated = sim[output_species][-1]

            if abs(theoretical - simulated) > eps:
                return {"success": False, "sim": sim, **species,
                        "theoretical": theoretical, "simulated": simulated}

        return {"success": True}



