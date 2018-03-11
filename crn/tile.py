import numpy as np
import random

from collections import defaultdict
from numpy.random import choice

class TileSystem:
    def __init__(self, *tiles, bonds=None, threshold=None, name=None):
        if name is None:
            name = id(self)

        self.name = name
        self.tiles = set(tiles)
        self.bonds = bonds
        self.threshold = threshold

    def simulate(self, seed, steps):
        def fit(tile, state, loc):
            '''
            returns True if 'tile' fits in 'loc' given 'state'
            '''
            B = 0
            x, y = loc
            for dx, dy, d in {(-1, 0, "w"), (1, 0, "e"),
                    (0, -1, "s"), (0, 1, "n")}:

                if (x + dx, y + dy) not in state:
                    continue

                if tile.matches_tile(state[(x + dx, y + dy)], from_dir=d):
                    B += self.bonds[tile.edges[d]]
                else:
                    return False

            if B >= self.threshold:
                return True
            return False

        if seed is Tile:
            seed = {(0, 0): seed}

        state = seed.copy()
        history = [State(state)]

        empty_locs = set()

        # fill initial empty locs
        for (x, y) in state:
            for dx, dy in {(-1, 0), (1, 0), (0, -1), (0, 1)}:
                loc = (x + dx, y + dy)
                if loc not in empty_locs and loc not in state:
                    empty_locs.add(loc)

        for i in range(0, steps):
            potential_tiles = defaultdict(list)

            for loc in empty_locs:
                for tile in self.tiles:
                    if fit(tile, state, loc):
                        potential_tiles[tile].append(loc)

            if not potential_tiles:
                print(f"simulation quit early: {i}/{steps} steps completed")
                break

            # select random fillable empty location to fill with tile
            pt_lst = list(potential_tiles.items())
            lengths = np.array(list(map(lambda x: len(x[1]), pt_lst)))
            tile = choice(list(map(lambda x: x[0], pt_lst)),
                    p=lengths / sum(lengths))
            selected_loc = (x, y) = random.choice(potential_tiles[tile])

            # update empty locs
            for dx, dy in {(-1, 0), (1, 0), (0, -1), (0, 1)}:
                neighbor = (x + dx, y + dy)
                if neighbor not in empty_locs and neighbor not in state:
                    empty_locs.add(neighbor)
            empty_locs.remove(selected_loc)

            # update state & history
            state[selected_loc] = tile
            history.append(State(state))

        return history

class Tile:
    opposite_dir = {"n": "s", "s": "n", "e": "w", "w": "e"}
    def __init__(self, name, n, s, w, e):
        self.name = name
        self.edges = {"n": n, "s": s, "w": w, "e": e}

    def matches_tile(self, tile, *, from_dir):
        a = self.edges[from_dir]
        b = tile.edges[Tile.opposite_dir[from_dir]]
        if None in (a, b):
            return False
        else:
            return a == b

class State:
    def __init__(self, state):
        self.state = state.copy()

    def __str__(self):
        maxx, maxy = minx, miny = list(self.state.keys())[1]
        for (x, y) in self.state:
            minx = min(minx, x)
            miny = min(miny, y)
            maxx = max(maxx, x)
            maxy = max(maxy, y)

        s = ''
        for y in range(maxy, miny - 1, -1):
            for x in range(minx, maxx + 1):
                s += self.state[x, y].name if (x, y) in self.state else '.'
            s += '\n'
        return s

