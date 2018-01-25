from crn import Reaction

class Expression:
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
        return Reaction(self, expr, 1)

    def __str__(self):
        return ' + '.join(
                map(lambda i: f"{i[1] if i[1] != 1 else ''}{i[0]}",
                    self.species.items()))

    def __repr__(self):
        return f"Expression({self.species})"

    def is_species(self):
        return len(self.species) == 1

    def get_species(self):
        return set(self.species.keys())

