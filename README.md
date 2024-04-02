# Basic Numerical Demonstration of the Genetic Algorithm

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

<details>
    <summary><strong>Execute Program with Command Line</strong></summary>

## Execute with Docker

From linux command line: <br />

```sh
cd genetic_algorithm_cpp
docker build -t ga .
docker run ga
```

## Execute if GCC installed

```sh
cd genetic_algorithm_cpp
g++ -o geneticAlgo src/main.cpp
./geneticAlgo
```

</details>

<details>
    <summary><strong>Issues Being Worked On</strong></summary>

#### Design

This program lacks quality architecture; however, this is being worked on. The geneticAlgorithm
function does take on the qualities of the Template design pattern in how each
function is ordered in such a way as to provide accurate data outputs. This program
was designed and constructed prior to me having taken Software Engineering or Software
Design courses. After refactoring one year later, I see odd arrangements with
private and public categorizations; however, the arrangement worked with the
sloppy design style. I am also attempting to implement suggestions from "Clean Code" as well.
So far, I've given variables better names.

#### BUG

In the event that there are 4 or less iterations before the optimal solution is found,
there will be excess generations produced on printout of the "The first generation
after initial population"; "The second generation after the initial population";
"The second to last generation"; "The last generation". E.g. if there are only 4
iterations for any run, "The second to last generation" will contain the optimal
fitness value of 10 chromosome and "The last generation" will contain an excess
generation (which might and might not contain the optimal chromosome since that
chromosome only has 10 spots in the chromosome roulette vector. A possible SOLUTION
is to include the condition of 'if there are less than 5 generations, then fill the
generations after the optimal solution generation with null input'.

</details>

<details>
    <summary><strong>Improvements To Be Made</strong></summary>

Utilize separate files for better code comprehension. Implement
helper classes for better code organization. Removal of unused functions for
final product. Overall better naming practices will help in understanding what each
variable and function does (e.g., I donw remember what `int endSpot = 0;` is doing).

</details>

<details>
    <summary><strong>Minor Note</strong></summary>

Only minor refactors were implemented after having completed this course. The
overall state of this program is left for reference of how I developed this program
during the class. I do plan on refactoring each function to be smaller (i.e., more
legible) and to have better names, which will lead to the removal of the many comments
littering the program.

</details>
