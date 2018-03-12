#  Law of Mass Action Reactions to ODE Converter

At the time when this code was originally written, I was working on a project
that required writing out several dozen ODE's that were generated from an
enzyme-kinetic reaction scheme.  To remove user error (and to work on a project
in Python), I wrote this set of short scripts to convert a plain text file with 
reaction schemes into a set of ODE's contained within a Matlab script.

While I later found better software that was already written (I highly recommend
the perl based [ALC](https://bmcsystbiol.biomedcentral.com/articles/10.1186/1752-0509-2-91)),
I still found this simple set of scripts to be pretty useful.

## Example Usage

Write out your reaction scheme in a comma separated list, one equation per line, with
the desired reaction rate names seperated by commas.  So for the following
scheme:
<p align="center">
  <img width="150" height="150" src="https://latex.codecogs.com/gif.latex?%5Cdpi%7B200%7D%20%5Cbegin%7Bmatrix%7D%20A%20&plus;%202B%20%5Cunderset%7Bk_2%7D%7B%5Cstackrel%7Bk_1%7D%7B%5Crightleftharpoons%7D%7D%20C%20%5C%5C%20B%20&plus;%20D%20%5Cstackrel%7Bk_1%7D%7B%5Crightarrow%7D%20E%20%5C%5C%203E%20&plus;%20A%20%5Cunderset%7Bk_4%7D%7B%5Cstackrel%7Bk_3%7D%7B%5Crightleftharpoons%7D%7D%20F%20%5Cend%7Bmatrix%7D" />
</p>
we might input the following into example.csv

```text
A + 2 * B = C , k_1 , k_3
B + D -> E , k_1
3*E + A <-> F, k_3, k_4
```

Running the command
```bash
make_stoic_matrix.py example.csv example_stoich.csv
```
will generate a [stoichmetric matrix](https://en.wikipedia.org/wiki/Stoichiometry#Stoichiometry_matrix) for the set of reactions, which can be used for a myriad of applications.  Here, we will only be using it to construct ODEs, but feel free to be creative.  The output from the above example would be
```text
0,k_1,k_3,k_1,k_3,k_4
A,-1,1,0,-1,1
B,-2,2,-1,0,0
C,1,-1,0,0,0
D,0,0,-1,0,0
E,0,0,1,-3,3
F,0,0,0,1,-1
```
and gets saved to example_stoich.csv.  The final step needed to construct our system of odes is to issue the command
```bash
make_ode_m_file.py example_stoich.csv example_ode.m
```
The created m-file seperates the ODEs into a driver function and a right-hand side for
easy modification.  Rates and intial conditions can be set within the m-file once it's created.

See the comments in the code for further details.
