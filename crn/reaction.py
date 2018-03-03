import re

from functools import reduce
from itertools import product
from operator import mul
from string import Formatter
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
            reactants = Expression({reactants: 1})
        if type(products) is Species:
            products = Expression({products: 1})

        self.reactants = reactants
        self.products = products
        self.coeff = float(k)
        self.is_schema = self.reactants.is_schema or self.products.is_schema

        if self.is_schema:
            self.schema_reactants = [r for r in self.reactants.species
                    if r.is_schema]

        for r in self.reactants.species:
            r.reactify()

    def __str__(self):
        return f"{self.reactants} -->({self.coeff}) {self.products}"
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

    def possible_reactions(self, state):
        """
        Given the current molecule/species counts in 'state', returns a list
        of non-schema reactions that are possible, but may have propensity
        zero.
        """
        rxns = []
        for possible_sps in product(state, repeat=len(self.schema_reactants)):
            reactants = Expression({})
            products = Expression({})
            groups = {}
            for sp, schema_r in zip(possible_sps, self.schema_reactants):
                match = schema_r.match(sp)
                if match is None:
                    break
                # TODO: debug label. remove this check when not debugging
                for key, val in match.groupdict().items():
                    if key in groups:
                        raise RuntimeError(
                                "Duplicate group name used in the same "
                                "reaction schema.")
                groups.update(match.groupdict())
                reactants += sp * self.reactants.species[schema_r]

            # Enters only if 'break' wasn't reach in the above loop
            else:
                for r, c in self.reactants.species.items():
                    if not r.is_schema:
                        reactants += r * c

                for p, c in self.products.species.items():
                    if p.is_schema:
                        p = Species(p.name.format(**groups))
                    products += p * c

                rxns.append((reactants >> products).k(self.coeff))
        return rxns

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
            return reduce(mul, (Symbol(s.name) - i for i in range(c)))

        return reduce(mul, map(flux_part, self.reactants.species.items()))

    def propensity(self, counts):
        """
        Returns the value of Reaction.discrete_flux given the currently
        present molecules/species in 'counts'.
        """
        def flux_part(i):
            s, c = i
            return reduce(mul, (counts.get(s, 0) - i for i in range(c)))

        return reduce(mul, map(flux_part, self.reactants.species.items()))

    def flux(self):
        """
        Returns a symbolic representation of the reaction rate of this
        reaction.
        """
        def flux_part(i):
            s, c = i
            return Symbol(s.name) ** c

        return self.coeff * reduce(mul,
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
        self.is_schema = any(sp.is_schema for sp in self.species)

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
    def __init__(self, name, is_schema=False, schema_groups=None):
        if name == "time":
            raise ValueError(
                "Failed to create Species 'time' because it is a reserved "
                "Species name. Please choose another name for this Species.")

        self.name = name
        self.is_schema = is_schema

        if self.is_schema:
            if schema_groups is None:
                schema_groups = {}

            self.schema_groups = schema_groups
            self.schema = name
            self.regex_schema = None
        elif schema_groups is not None:
            raise ValueError(
                    "Species constructor passed 'schema_groups' but "
                    "'is_schema' is False.")

    def reactify(self):
        if not self.is_schema:
            return

        format_keys = [i[1] for i in Formatter().parse(self.schema) if i[1]]
        format_dict = {}
        for key in format_keys:
            if key in self.schema_groups:
                format_dict[key] = f"(?P<{key}>{self.schema_groups[key]})"
            else:
                raise RuntimeError(
                        "Schema used in reaction without substituting all "
                        "captured groups")

        self.schema = self.schema.format(**format_dict)
        self.compile_regex()

    def compile_regex(self):
        try:
            self.regex_schema = re.compile(f"^{self.schema}$")
        except:
            print(self.schema)

    def match(self, species):
        if not self.regex_schema:
            raise RuntimeError(
                    "Species.match called before Species.reactify")

        return self.regex_schema.match(species.name)

    def has_groups(self):
        if not self.is_schema:
            return False
        return bool([i[1] for i in Formatter().parse(self.schema) if i[1]])

    def __call__(self, *args):
        args = list(args)
        if not self.is_schema:
            raise RuntimeError("Non-Schema Species cannot be called")

        format_keys = [i[1] for i in Formatter().parse(self.schema) if i[1]]

        if len(args) > len(format_keys):
            raise RuntimeError("Too many arguments called on Schema")

        args += [None] * (len(format_keys) - len(args))

        schema_groups = self.schema_groups.copy()
        format_dict = {}
        for arg, key in zip(args, format_keys):
            if arg is None:
                format_dict[key] = ""
            elif type(arg) is str:
                format_dict[key] = f"{{{arg}}}"
                if key in schema_groups:
                    schema_groups[arg] = schema_groups[key]
                    del schema_groups[key]
            else:
                # Can be used with custom objects and non-strings
                format_dict[key] = f"{arg}"

        ret = Species(self.schema.format(**format_dict), is_schema=True,
                schema_groups=schema_groups)
        ret.compile_regex()

        return ret

    def __add__(self, other):
        if type(other) is Expression:
            return other + Expression({self: 1})
        elif type(other) is Species:
            return Expression({self: 1}) + Expression({other: 1})

        return NotImplemented

    __radd__ = __add__

    def __rshift__(self, other):
        return Reaction(self, other)

    def __rrshift__(self, other):
        return Reaction(other, self)

    def __mul__(self, other):
        if type(other) is int:
            return Expression({self: other})

        return NotImplemented

    def __rmul__(self, other):
        if type(other) is int:
            return Expression({self: other})

        return NotImplemented

    def __repr__(self):
        return self.name

    def __str__(self):
        return self.name

    def __eq__(self, other):
        if type(other) is not Species:
            return NotImplemented

        if self.is_schema and other.is_schema:
            return ((self.name == other.name) and
                    (self.schema_groups == other.schema_groups))
        elif not self.is_schema and not other.is_schema:
            return self.name == other.name
        else:
            False

    def __lt__(self, other):
        if type(other) is Species:
            if not self.is_schema and not other.is_schema:
                return self.name < other.name
            elif self.is_schema and not other.is_schema:
                return False
            elif not self.is_schema and other.is_schema:
                return True
            else:
                if self == other:
                    return False
                return (self.name, id(self)) < (other.name, id(other))

        return NotImplemented

    def __hash__(self):
        return hash((self.name, self.is_schema,
            tuple(sorted(self.schema_groups.items()))
                if self.is_schema else None))

    __req__ = __eq__

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
    """
    species = species.split()
    if len(species) == 1:
        return Species(species[0])

    if "nothing" in species:
        raise ValueError(
            "Species 'nothing' is reserved and therefore cannot be created \n"
            "using `species` function. Use '0' in your reactions instead. \n"
            "For example,\n\n"
            "    0 >> a\n"
            "    a + b >> 0"
            "\n\nOr create 'nothing' directly using the Species constructor.")

    return map(Species, species)

def schemas(schemas, schema_groups=None):
    if schema_groups is None:
        schema_groups = {}
    schema_groups = schema_groups.copy()

    schemas = schemas.split()
    if len(schemas) == 1:
        return Species(schemas[0], is_schema=True, schema_groups=schema_groups)

    return (Species(schema, True, schema_groups) for schema in schemas)

