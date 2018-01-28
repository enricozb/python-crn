import operator

from functools import reduce
from sympy import Symbol

class Reaction:
    """
    Representation of a stoichiometric reaction using a pair of Expressions,
    one for the reactants and one for the products.

    args:
        reactants: Expression
            The left hand side of the stoichiometric equation
        products: Expression
            The right hand side of the stoichiometric equation
        k: float
            The rate constant of the reaction

    properties:
        reactants: Expression
            The left hand side of the stoichiometric equation
        products: Expression
            The right hand side of the stoichiometric equation
        coeff: float
            The rate constant of the reaction
    """
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
        """
        Changes the reaction coefficient to `coeff` and returns `self`.

        args:
            coeff: float
                The new reaction coefficient

        This is useful for including the rate constant during the construction
        of a reaction. For example

            x, y, z = species("X Y Z")
            sys = CRN(
                (x + y >> z).k(2.5),
                (z >> x).k(1.5),
                (z >> y).k(0.5))
            ...
        """
        self.coeff = coeff
        return self

    def get_species(self):
        """
        Returns the set of species present in the products and reactants.
        """
        return {
            *self.reactants.get_species(),
            *self.products.get_species()
        }

    def net_production(self, species):
        """
        Returns the net stoichiometric coefficient of a species in this
        reaction.

        args:
            species: str
                string name of the species
        """
        return (self.products.species.get(species, 0) -
                self.reactants.species.get(species, 0))


    def flux(self):
        """
        Returns a symbolic representation of the reaction rate of this
        reaction.
        """
        def flux_part(i):
            s, c = i
            return Symbol(s) ** c

        return self.coeff * reduce(operator.mul,
                map(flux_part, self.reactants.species.items()))

