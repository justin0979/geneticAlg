# geneticAlg

This program simulates the Genetic Algorithm's ability to progressively match the given target
of 1010101010 from 20 randomly initialized "chromosomes" (represented by integers 0,1 arrays).
The algorithm uses a percentage value (Pco 0.9, 0.7, 0.5, 0.3, 0.0) to determine the number
of "chromosomes" that will be replicated. The remaining "chromosomes" will undergo a simulated
cross-over process, based on the mask 0000011111 (or even split) with random "chromosome" pairings.
"Chromosome" fitness values are calculated by the number of gene (array index) matches to the
target pattern. Random mutations will occur with every Pco value. Without the random mutation of at
least one "gene" in the 20 "chromosome" population, the Pco 0.0 run would never stop, due to 100% of
of the "chromosomes" just being replicated. The algorithm will run until at least one chromosome matches
the target by displaying the fitness value of 10.

Program will display the initial generated population and all iterations up to the target value,
followed by the 2nd and 3rd iterations and then the second to last and last iterations for
Pco's 0.7 and 0.0.

# Execute with Docker

From linux command line: <br />
`docker build -t ga .` <br />
`docker run ga`

# Execute if GCC installed

`gcc -o geneticAlgo geneticalgorithm.cpp -lstdc++` <br/>
`./geneticAlgo`
