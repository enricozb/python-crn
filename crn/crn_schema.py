from crn import Species, Simulation
from itertools import product
from numpy import log
from numpy.random import choice
from random import uniform

class CRNSchema:
    def __init__(self, *system, **kwargs):
        self.system = system
        self.name = kwargs.get("name", id(self))

    def stoch_simulate(self, state, t=None, steps=None):
        def reactions(state):
            rxns = []

            state = {k: v for k, v in state.items() if v > 0}

            for rxn in self.system:
                reactants = rxn.reactants
                schema = list(reactants.schemas)
                num_reactants = len(schema)

                for potential_species in product(state, repeat=num_reactants):
                    match = rxn.subs(dict(zip(schema, potential_species)))
                    if match is not None:
                        rxns.append(match)

            return rxns

        def propensities(rxns, state):
            return [r.propensity(state) for r in rxns]

        if t is not None and steps is not None:
            raise ValueError(
                    "CRNSchema.stoch_simulate both `t and `steps` args were "
                    "passed. At most one of the arguments can be passed.")

        if t is None and steps is None:
            steps = 1000
            t = 50000

        if steps is None:
            steps = float('inf')
        if t is None:
            t = float('inf')

        sim = {"time": [0]}

        for sp, count in state.items():
            sim[sp.name] = [count]

        state = state.copy()

        curr_time = 0
        curr_step = 0

        while True:
            if curr_time >= t or curr_step >= steps:
                break

            rxns = reactions(state)
            props = propensities(rxns, state)
            p_tot = sum(props)
            if p_tot == 0:
                print("simulation ended before reaching time or steps")
                break

            dt = -log(uniform(0, 1)) / p_tot
            reaction_to_occur = choice(rxns, p=[p/p_tot for p in props])
            print(f"{curr_time:10.5f}\t{reaction_to_occur}")

            for r, coeff in reaction_to_occur.reactants.species.items():
                sr = Species(r)
                if state.get(sr, 0) - coeff < 0:
                    raise RuntimeError("Impossible reaction chosen")

                state[sr] -= coeff

            for p, coeff in reaction_to_occur.products.species.items():
                sp = Species(p)
                if p not in state:
                    state[sp] = 0
                    sim[p] = [0] * (curr_step + 1)

                state[sp] += coeff

            for sp, count in state.items():
                sim[sp.name].append(count)

            curr_time += dt
            curr_step += 1

            sim["time"].append(curr_time)

        return Simulation(sim)

