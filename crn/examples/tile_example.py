from crn import *

z1 = Tile("1", 1, 1, "nop", "nop")
z2 = Tile("0", 0, 0, "nop", "nop")
z3 = Tile("0", 0, 1, "inc", "inc")
i = Tile("1", 1, 0, "nop", "inc")
bv = Tile("V", "B", "B", "inc", None)
bh = Tile("H", 0, None, "B", "B")
b3 = Tile("B", "B", None, "B", None)

bonds = {0: 1, 1: 1, "B": 2, "inc": 1, "nop": 1}
threshold = 2

sys = TileSystem(z1, z2, z3, i, bv, bh, b3, bonds=bonds, threshold=threshold)

seed = {(0, 0): b3, (-1, 0): bh, (0, 1): bv}
history = sys.simulate(seed, steps=1000)

state = history[-1]

for y in range(60, -1, -1):
    for x in range(-60, 1):
        if (x, y) in state:
            print(state[x, y].name, end='')
        else:
            print('.', end='')
    print()

