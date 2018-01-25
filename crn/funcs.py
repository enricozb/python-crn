from crn import Expression

def species(species):
    for s in species.split():
        yield Expression({s:1})

