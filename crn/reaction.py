import operator

from functools import reduce
from sympy import Symbol

class Reaction:
    def __init__(self, reactants, products, k):
        self.reactants = reactants
        self.products = products
        self.coeff = float(k)

    def __str__(self):
        rcts_str = str(self.reactants)
        return (f"{' ' * len(rcts_str)} {self.coeff:.1f} \n"
                f"{self.reactants} ---> {self.products}")

    def __repr__(self):
        return (f"Reaction({repr(self.reactants)}, {repr(self.products)}, "
                f"{self.coeff})")

    def k(self, coeff):
        self.coeff = coeff
        return self

    def get_species(self):
        return {
            *self.reactants.get_species(),
            *self.products.get_species()
        }

    def net_production(self, species):
        return (self.products.species.get(species, 0) -
                self.reactants.species.get(species, 0))


    def flux(self):
        def flux_part(i):
            s, c = i
            return Symbol(s) ** c

        return self.coeff * reduce(operator.mul,
                map(flux_part, self.reactants.species.items()))

