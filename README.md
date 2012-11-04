qmc
===

`qmc` is a naive implementation of the [Quine-McCluskey algorithm](https://en.wikipedia.org/wiki/Quine%E2%80%93McCluskey_algorithm). It can be used optimize two level multi-output boolean functions. The algorithm implemented here is described in detail in the book "Reti Logiche" by C. Bolchini et al, published by Apogeo.

Given the ONset (the set of minterms - input for which the output is 1) and the DCset (don't care set - inputs for which it is not relevant the output) of the function to be optimized, the algorithm will output the boolean form of the optimized function.

The bibliographic reference for the algorithm hereby implemented is C.Bolchini et al, "Reti Logiche", Apogeo.

Usage
-----

A tiny GUI written using PyQt is provided (`gui.py`). Simply insert the onset and the dcset of each function to be minimized, as lists of integers separed by a comma and click "GO". For example, the XOR function `F = a'b + ab'` has the following truth table
	    
	    a   b   F(a,b)
	0:  0   0     0
	1:  0   1     1
	2:  1   0     1
	3:  1   1     0

and can be represented as `ON(1,2)` plus an empty DC as the function is defined for all the input variables. The number of variables is inferred by the maximum value in the ONset or in the DCset.

The algorithm can be used directly from within python code. An example of use is the following. To minimize the three-output combinatorial network

    F1 = ON(3,4,6,7) DC(0,2)
    F2 = ON(3,4,6) DC(2,5)
    F3 = ON(0,2) DC(4,6)

you can write the following code
```python
q = QuineMcCluskey([
                   QmcFunction([3,4,6,7], [0,2]),
                   QmcFunction([3,4,6], [2,5]),
                   QmcFunction([0,2], [4,6])
                   ], LITERALS_COST_FUNCTION)
print ' --- [Expansion] Prime implicants: ---'
q.findPrimeImplicants()
q.printImplicantList()
q.printTable()
print 'Solution found:'
print q.simplify()
print 'Cost:', q.sol.getCost()
```
and this is what you'll get on the console:

     --- [Expansion] Prime implicants: ---
      A = a'b  	110 	2 3 
      B = bc'  	111 	2 6 
      C = ab'  	010 	4 5 
      D = ac'  	111 	4 6 
      E = c'  	101 	0 2 4 6 
      F = b  	100 	2 3 6 7 
               3     4     6     7   |   3     4     6   |   0     2   |  C
      A -> |   X                     |   X               |             |  2
      B -> |               X         |               X   |         X   |  2
      C -> |                         |         X         |             |  2
      D -> |         X     X         |         X     X   |             |  2
      E -> |         X     X         |                   |   X     X   |  1
      F -> |   X           X     X   |                   |             |  1
    Solution found:
    F1 = b + D
    F2 = D + a'b
    F3 = c'
    D = ac'

    Cost: 7

By setting `VERBOSE = True`, you can output to the console all the intermediate steps of the algorithm as it is executed.

Known Bugs
----------

* The implicants cost functions is different than the one suggested in the reference book
