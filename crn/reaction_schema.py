import re

from crn import Species 
from string import Formatter

# TODO: ExpressionSchema can't hold just names to the schemas as the keys,
# since the SpeciesSchemas have properties.

class ReactionSchema:
    def __init__(self, reactants, products, k=1):
        if type(reactants) is SpeciesSchema:
            reactants = ExpressionSchema({reactants: 1})
        if type(products) is SpeciesSchema:
            products = ExpressionSchema({products: 1})

        self.reactants = reactants.regex_reactants()
        self.products = products
        self.rate_constant = k

    def k(self, rate_constant):
        self.rate_constant = rate_constant
        return self

    def regex_reactants(self, names):
        self.reactants.regex_reactants()

    def subs(self, sub_dict):
        """
        Creates a Reaction object by substituting the schemas for proper
        species
        """
        named_groups = {}
        new_reactants = None
        new_products = None

        for reactant, coeff in self.reactants.schemas.items():
            species = sub_dict[reactant]
            match = reactant.match(species)
            if match is None:
                return None

            for group, val in match.groupdict():
                if group in named_groups and named_groups[group] != val:
                    raise Exception("Duplicate group names")

            named_groups.update(match.groupdict())
            if new_reactants is None:
                new_reactants = coeff * species
            else:
                new_reactants += coeff * species

        for product, coeff in self.products.schemas.items():
            if new_products is None:
                new_products = coeff * product.format(named_groups)
            else:
                new_products += coeff * product.format(named_groups)

        return (new_reactants >> new_products).k(self.rate_constant)


    def __str__(self):
        return f"{self.reactants} ---> {self.products}"

    def __repr__(self):
        return (f"ReactionSchema({repr(self.reactants)}, {repr(self.products)}, "
                f"{self.rate_constant})")


class ExpressionSchema:
    def __init__(self, schemas):
        self.schemas = schemas

    def regex_reactants(self):
        new_schemas = {}
        for s, coeff in self.schemas.items():
            s.regex_reactants()
            new_schemas[s] = coeff

        return ExpressionSchema(new_schemas)

    def __add__(self, other):
        if type(other) is ExpressionSchema:
            schemas_copy = self.schemas.copy()
            for s, c in other.schemas.items():
                if s not in schemas_copy:
                    schemas_copy[s] = 0
                schemas_copy[s] += c

            return ExpressionSchema(schemas_copy)

        return NotImplemented

    def __rmul__(self, coeff):
        if type(coeff) is int:
            schemas_copy = {}
            for s, c in self.schemas.items():
                schemas_copy[s] = c * coeff

            return ExpressionSchema(schemas_copy)

        return NotImplemented

    __mul__ = __rmul__

    def __rshift__(self, expr):
        return ReactionSchema(self, expr)

    def __str__(self):
        return ' + '.join(
                map(lambda i: (f"{i[1]} * " if i[1] != 1 else '') + f"{i[0]}",
                    self.schemas.items()))

    __repr__ = __str__


class SpeciesSchema:
    def __init__(self, schema, regexd=None, names=None):
        if regexd is None:
            regexd = schema
        if names is None:
            names = {}

        self.schema = schema
        self.regexd = regexd
        self.names = names.copy()

    def clones(self, n):
        for i in range(n):
            yield SpeciesSchema(self.schema, self.regexd)

    def format(self, names):
        return Species(self.regexd.format(**names))

    def match(self, species):
        return self.regexd_compiled.match(species.name)

    def regex_reactants(self):
        format_keys = set(
                i[1] for i in Formatter().parse(self.regexd)
                if i[1] is not None)

        named_groups = format_keys - set(self.names)
        named_groups = {name: f"{{{name}}}" for name in named_groups}

        names = {name: f"(?P<{name}>{p})" for name, p in self.names.items()}

        self.regexd = self.regexd.format(**named_groups, **names)
        self.regexd_compiled = re.compile('^' + self.regexd + '$')

    def __getitem__(self, args):
        if type(args) is tuple:
            args = list(args)
        else:
            args = [args]

        format_keys = [
                i[1] for i in Formatter().parse(self.regexd)
                if i[1] is not None]

        if len(format_keys) > len(args):
            args = list(args) + (len(format_keys) - len(args)) * [None]


        new_names = {}
        format_names = {}
        for old, new in zip(format_keys, args):
            if type(new) is str:
                new_names[new] = self.names[old]
                format_names[old] = f"{{{new}}}" if new else ""
            elif type(new) is int:
                format_names[old] = str(new)
            elif new is None:
                format_names[old] = ""
            else:
                raise Exception(
                        "SpeciesSchema.__getitem__: Only int and str are "
                        "supported")

        new_regexd = self.regexd.format(**format_names)

        return SpeciesSchema(self.schema, new_regexd, new_names)

    def __add__(self, other):
        if type(other) is ExpressionSchema:
            return other + ExpressionSchema({self: 1})
        elif type(other) is SpeciesSchema:
            return ExpressionSchema({self: 1, other: 1})

        return NotImplemented

    def __rshift__(self, other):
        return ReactionSchema(self, other)

    def __mul__(self, other):
        if type(other) is int:
            return ExpressionSchema({self: other})

    __radd__ = __add__
    __rmul__ = __mul__

    def __repr__(self):
        return self.regexd

    def __str__(self):
        return self.regexd

def species_schemas(args, names=None):
    args = args.split()
    if len(args) == 1:
        return SpeciesSchema(args[0], names=names)

    return (SpeciesSchema(sp, names=names) for sp in args)

