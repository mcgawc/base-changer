# base-changer
NEBaseChanger knock off

Primer melting temperature calculator that uses nearest-neighbor thermodynamics (SantaLucia, Hicks. 2004) to determine Tm and takes
into account the buffer conditions of Q5 Hot Start Polymerase. Built to mimic NEBaseChanger.

But...
Let's start by building a working Tm calculator that takes into account mismatches

Then build a primer designer, that can construct potential primer pairs according to the below rules:
    -5' end forward primer has 10 complementary nucleotides above the substitution
    -3' end of forward can be extended as needed to increase Tm
    -5' of reverse primer is anchored to where the 5' of forward ends, i.e. no overlap
    -3' end of reverse can be extended to increase Tm
    -Minimum reverse primer length is ~18 nucleotides I think

Then the program must test the Tm of those primer pairs to find the best solution, meaning the primers
have Tms within 5 degrees of each other, are not too long, etc.

Finally the program should be made to work with bulk input, and output a table containing
the primer sequences, both forward and reverse, and the optimal annealing temperature. 
