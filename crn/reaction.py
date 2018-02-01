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
    def __init__(self, reactants, products, k=1):
        if reactants == 0:
            reactants = Species("nothing")

        if products == 0:
            products = Species("nothing")

        if type(reactants) not in (Species, Expression):
            raise ValueError(
                "Attempted construction of reaction with type of reactants "
                f"as {type(reactants)}. Type of reactants must be Species "
                "or Expression")
        if type(products) not in (Species, Expression):
            raise ValueError(
                "Attempted construction of products with type of products "
                f"as {type(products)}. Type of products must be Species "
                "or Expression")

        if type(reactants) is Species:
            reactants = Expression({reactants.name: 1})
        if type(products) is Species:
            products = Expression({products.name: 1})

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

    def discrete_flux(self):
        """
        Discrete analog of Reaction.flux: Returns a symbolic representation
        of the discrete/stochastic reaction rate of this reaction. Essentially
        the propensity not including the rate constant.
        """
        def flux_part(i):
            s, c = i
            return reduce(operator.mul, (Symbol(s) - i for i in range(c)))

        return reduce(operator.mul,
                map(flux_part, self.reactants.species.items()))

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


class Expression:
    """
    Class used for very basic symbolic manipulation of left/right hand
    side of stoichiometric equations. Not very user friendly; users should
    just use the `species` functions and manipulate those to get their
    reactions.

    args:
        species: Dict[str, int]
            represents species (string names) and their coefficients (ints)
            all added together.

    properties:
        species: Dict[str, int]
            represents species (string names) and their coefficients (ints)
            all added together. The same as the argument passed to the
            constructor
    """
    def __init__(self, species):
        self.species = species

    def __add__(self, other):
        if type(other) is Expression:
            species_copy = self.species.copy()
            for s, c in other.species.items():
                if s not in species_copy:
                    species_copy[s] = 0
                species_copy[s] += c
            return Expression(species_copy)

        return NotImplemented

    def __rmul__(self, coeff):
        if type(coeff) is int:
            species_copy = {}
            for s, c in self.species.items():
                species_copy[s] = c * coeff

            return Expression(species_copy)

        return NotImplemented

    __mul__ = __rmul__

    def __rshift__(self, expr):
        return Reaction(self, expr)

    def __str__(self):
        return ' + '.join(
                map(lambda i: f"{i[1] if i[1] != 1 else ''}{i[0]}",
                    self.species.items()))

    def __repr__(self):
        return ' + '.join(
                map(lambda i: f"{i[1] if i[1] != 1 else ''}{i[0]}",
                    self.species.items()))

    def get_species(self):
        """
        Returns the names of the species in this expression, not their
        coefficients.
        """
        return set(self.species.keys())


class Species:
    def __init__(self, name):
        if name == "time":
            raise ValueError(
                "Failed to create Species 'time' because it is a reserved "
                "Species name. Please choose another name for this Species.")
        self.name = name

    def __add__(self, other):
        if type(other) is Expression:
            return other + Expression({self.name: 1})
        elif type(other) is Species:
            return Expression({self.name: 1, other.name: 1})

        return NotImplemented

    def __radd__(self, other):
        if type(other) is Expression:
            return other + Expression({self.name: 1})
        elif type(other) is Species:
            return Expression({self.name: 1, other.name: 1})

        return NotImplemented

    def __rshift__(self, other):
        return Reaction(self, other)

    def __mul__(self, other):
        if type(other) is int:
            return Expression({self.name: other})

        return NotImplemented

    def __rmul__(self, other):
        if type(other) is int:
            return Expression({self.name: other})

        return NotImplemented

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

def species(species):
    """
    Create a list of Species (Single species Expressions).

    args:
        species: str
            A space-seperated string representing the names of the species
            being created

    This is normally used like this:

        x, y, z = species("X Y Z")
        rxn = x + y >> z
        ...

    The names MUST be valid Python identifiers: "X0" is valid but "0X" is not.

    The name 'nothing' refers to a Species which will always have a
    concentration of 1. This is useful for spontaneous decay or creation of
    species. For example,

        x, y, nothing = species("X Y nothing")

        sys = CRN(
            x + y >> nothing,
            nothing >> x,
            nothing >> y)
    """
    species = species.split()
    if "nothing" in species:
        raise ValueError(
            "Species 'nothing' is reserved and therefore cannot be created \n"
            "using `species` function. Use '0' in your reactions instead. \n"
            "For example,\n\n"
            "    0 >> a\n"
            "    a + b >> 0"
            "\n\nOr create 'nothing' directly using the Species constructor.")

    return map(Species, species)

