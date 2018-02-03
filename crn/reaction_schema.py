from string import Formatter

# TODO: ExpressionSchema can't hold just names to the schemas as the keys,
# since the SpeciesSchemas have properties.

class ReactionSchema:
    def __init__(self, reactants, products, k=1):
        if reactants == 0:
            reactants = SpeciesSchema("nothing")

        if products == 0:
            products = SpeciesSchema("nothing")

        self.reactants = reactants
        self.products = products
        self.rate_constant = k

    def k(self, rate_constant):
        self.rate_constant = rate_constant
        return self

    def regexify(self, names):
        self.reactants.regexify(names, "reactants")
        self.products.regexify(names, "products")

class ExpressionSchema:
    def __init__(self, schemas):
        self.schemas = schemas

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
    def __init__(self, schema):
        self.schema = schema
        self.regexd = schema

    def regexify(self, names, prods_or_react):
        if prods_or_react not in {"products", "reactants"}:
            # TODO: Actually write the value error
            raise ValueError("SpeciesSchema.regexify: prods_or_react")

        format_keys = set(
                i[1] for i in Formatter().parse(self.schema)
                if i[1] is not None)

        named_groups = format_keys - set(names)

        if prods_or_react == "products":
            named_groups = {name: f"{{{name}}}" for name in (named_groups)}
        if prods_or_react == "reactants":
            named_groups = {name: f"?P<{name}>" for name in (named_groups)}

        self.regexd = self.schema.format(**named_groups, **names)

    def __add__(self, other):
        if type(other) is ExpressionSchema:
            return other + ExpressionSchema({self.schema: 1})
        elif type(other) is SpeciesSchema:
            return ExpressionSchema({self.schema: 1, other.schema: 1})

        return NotImplemented

    def __rshift__(self, other):
        return ReactionSchema(self, other)

    def __mul__(self, other):
        if type(other) is int:
            return ExpressionSchema({self.schema: other})

    __radd__ = __add__
    __rmul__ = __mul__

    def __repr__(self):
        return self.schema

    def __str__(self):
        return self.schema

def species_schemas(*args):
    yield from map(SpeciesSchema, args)


