import numpy as np

from crn import Expression, Simulation
from random import random
from scipy.integrate import odeint

class CRN:
    def __init__(self, *system):
        self.system = system
        self.species = self.get_species()
        self.index = self.get_species_index()
        self.diffeq_system_func = self.rate_laws()

    def get_species(self):
        species = set()
        for s in self.system:
            species |= s.get_species()
        return species

    def get_species_index(self):
        return dict(enumerate(sorted(self.species)))

    def rate_law_for_species(self, s):
        if type(s) is Expression:
            if not s.is_species():
                raise ValueError(
                        "rate_law_for_species called on complex expression")
            s, *_ = s.species.keys()

        return sum([rxn.net_production(s) * rxn.flux() for rxn in self.system])

    def rate_laws(self):
        laws = []

        # Has to be iterated this way to laws is in order
        for i in range(len(self.index)):
            sp = self.index[i]
            law = self.rate_law_for_species(sp)
            for i, s in self.index.items():
                sub = 1 if s == "nothing" else f"v[{i}]"
                law = law.subs(s, sub)
            laws.append(str(law))

        func = f"lambda v, t: [{', '.join(laws)}]"
        return eval(func)

    def simulate(self, conc, t=20):
        t = np.linspace(0, t, 100)

        conc_temp = {}

        for s, c in conc.items():
            if type(s) is Expression:
                if not s.is_species():
                    raise ValueError("concentrations must only be for "
                            "species, not complex expressions")
                s, *_ = s.species.keys()

            elif type(s) is not str:
                raise ValueError(
                    "concentrations must have key type 'str' or 'Expression'"
                    f" not '{type(s)}'")

            conc_temp[s] = c

        conc = conc_temp

        v0 = [0] * len(self.species)

        for i, s in self.index.items():
            v0[i] = conc.get(s, 0)

        sol = odeint(self.diffeq_system_func, v0, t)

        sol_dict = {'time': t}

        for i, s in self.index.items():
            sol_dict[s] = sol[:, i]

        return Simulation(sol_dict)

    def validate(self, func, *, input_species, output_species, N=100,
            eps=1e-2, t=500):
        for i in range(N):
            species = {sp : random() * 10 for sp in input_species}
            theoretical = func(species)
            sim = self.simulate(species, t=t)
            simulated = sim[output_species][-1]

            if abs(theoretical - simulated) > eps:
                print(abs(theoretical - simulated), eps)
                return {"success": False, "sim": sim, **species,
                        "theoretical": theoretical, "simulated": simulated}

        return {"success": True}



