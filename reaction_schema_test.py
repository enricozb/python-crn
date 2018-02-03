from crn import *

dna, rna = species_schemas(
        "DNA[({left}{B}*){promoter}({middle}{B}*){terminator}({right}{B}*)]",
        "RNA[{middle}]")

names = {"B": "[ATCG]", "promoter": "GAGGAG", "terminator": "TATTAT"}

rxn = 2 * dna >> rna
rxn.regexify(names)

print(rxn)


