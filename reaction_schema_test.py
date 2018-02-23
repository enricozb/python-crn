from crn import *

names = {"rest": "[01]*"}


Stack1 = species_schemas("Stack1<{rest}{top}>", names=names)
Stack2 = species_schemas("Stack2<{rest}{top}>", names=names)
Stack3 = species_schemas("Stack3<{rest}{top}>", names=names)

s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12 = species_schemas(
        "s1 s2 s3 s4 s5 s6 s7 s8 s9 s10 s11 s12")

sys = CRNSchema(
        s1 + Stack1["r1", 1] >> s2 + Stack1["r1"],
        s1 + Stack1[None] >> s12 + Stack1[None],
        s2 + Stack2["r2", 0] >> s3 + Stack2["r2"],
        s2 + Stack2["r2", 1] >> s4 + Stack2["r2"],
        s2 + Stack2[None] >> s3 + Stack2[None],
        s3 + Stack2["r2"] >> s1 + Stack2["r2", 1],
        s4 + Stack3["r3"] >> s6 + Stack3["r3", 0],
        s6 + Stack2["r2", 1] >> s4 + Stack2["r2"],
        s6 + Stack2["r2", 0] >> s8 + Stack2["r2"],
        s6 + Stack2[None] >> s8 + Stack2[None],
        s8 + Stack3["r3"] >> s7 + Stack3["r3", 1],
        s7 + Stack2["r2", 0] >> s5 + Stack2["r2"],
        s7 + Stack2["r2", 1] >> s8 + Stack2["r2"],
        s7 + Stack2[None] >> s9 + Stack2[None],
        s5 + Stack3["r3"] >> s7 + Stack3["r3", 0],
        s9 + Stack3["r3", 0] >> s10 + Stack3["r3"],
        s9 + Stack3["r3", 1] >> s11 + Stack3["r3"],
        s9 + Stack3[None] >> s1 + Stack3[None],
        s10 + Stack2["r2"] >> s9 + Stack2["r2", 0],
        s11 + Stack2["r2"] >> s9 + Stack2["r2", 1])


s1_init, s2_init, s3_init, s1 = species(
        "Stack1<1111111111> Stack2<> Stack3<> s1")

sim = sys.stoch_simulate({s1_init: 1, s2_init: 1, s3_init: 1, s1: 1}, steps=1000)

