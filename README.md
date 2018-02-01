# Python CRN
A CRN simulator based on
[David Soloveichik's Mathematica CRN Simulator](
http://users.ece.utexas.edu/~soloveichik/crnsimulator.html)

## Version & Promises
Currently this is in a pre-alpha state, but already has some reasonable
features. No guarantees are made for backward-compatibility. That is, some
of the API calls might change, some of the reaction literals syntax might
change, etc.

## Dependencies
The package is written for Python 3.6+ and definitely won't work on an older
Python version.

The current list of dependencies are
```
numpy
scipy
sympy
stochpy
```

## Examples
A simple example of creating a crn, simulating it, and plotting it, is given
here. An explanation of this code in detail can be found in
`crn/examples/basic_example.py`.

```python
from crn import *

a, a1, a2, b, c, t, z = species("A A1 A2 B C T Z")

sys = CRN(
    a >> a1 + a2,
    a1 + b >> t,
    c >> z,
    (a2 >> z).k(2.5),
    z + t >> 0)

sys.simulate({a: 2.5, b: 2.0}, t=5).plot("sim.png", title="Example Simulation")
```

More examples will be added soon in the `crn/examples/` folder.


## Installation
Install it by placing the `crn/` folder in your python path. I'll setup an
install script soon, or I'll add it to `pip` or something idk.

