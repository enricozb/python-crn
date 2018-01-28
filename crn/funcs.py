from crn import Expression

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
    for s in species.split():
        yield Expression({s:1})

