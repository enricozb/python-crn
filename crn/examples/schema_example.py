from crn import *

# states and initial stacks
s1, s2, s3, halt = species("s1 s2 s3 halt")

# schema group patterns
schema_groups = {"rest": "[01]*", "top": "[01]"}

# stacks
Stack1, Stack2 = schemas("Stack1<{rest}{top}> Stack2<{rest}{top}>",
        schema_groups)

sys = CRN(
    s1 + Stack1() >> halt + Stack1(),
    s1 + Stack1("r1", 1) >> s2 + Stack1("r1"),
    s1 + Stack1("r1", 0) >> s3 + Stack1("r1"),
    s2 + Stack2("r2") >> s1 + Stack2("r2", 1),
    s3 + Stack2("r2") >> s1 + Stack2("r2", 0))


sim = sys.schema_simulate({s1: 1, Stack1(101010): 1, Stack2(): 1})
for rxn in sim.reactions:
    print(rxn)

