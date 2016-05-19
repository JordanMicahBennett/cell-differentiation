cell-differentiation
====================

Cell-Diff - Exact computation of cell differentiation model trees 

Version: 1.0.0
Date: 2016-05-18

Authors:
Joshua L. Phillips <Joshua.Phillips@mtsu.edu>
Department of Computer Science, Middle Tennessee State University

Michael E. Colvin <mcolvin@ucmerced.edu>
School of Natural Sciences, University of California, Merced

Thanks for considering this code for your cell differentiation needs!

The development of Cell-Diff is mainly funded by academic research grants.
To help us fund development, we humbly ask that you cite the poster abstract:

* Analytic parameter fitting in stochastic stem cell models.
  J. L. Phillips, J. E. Manilay, and M. E. Colvin
  Biophysical Journal 98(3), 739a. (2010)
  doi:10.1016/j.bpj.2009.12.4052

This work was supported by NSF Award Number 0960480 and was supported
in part by NIH Grants GM077520 and in part by the U.S. Department of
Energy, Office of Science, Offices of Advanced Scientific Computing
Research, and Biological & Environmental Research through the
U.C. Merced Center for Computational Biology.

*************
GENERAL INFO:
*************

Cell-Diff is free software, distributed under the GNU General Public License. 

If you want to distribute a modified version or use part of Cell-Diff
in your own program, remember that the entire modified code must be licensed 
under GPL, and that it must clearly be labeled as derived work. It should 
not use the name "Cell-Diff", and make sure support questions are
directed to you instead of the Cell-Diff developers.

Cell-diff is a suite of tools for performing analyic parameter fitting
for stochastic cell differentiation models, such as the classic Till Model
[A stochastic model of stem cell proliferation, based on the growth of spleen
colony-forming cells. J.E. Till, E.A. McCulloch, and L. Siminovitch,
PNAS 51:29-36 (1964)]

*******************
Basic Documentation
*******************

The code is designed to aid parameter estimation when fitting cell
differentiation models to experimental cell differentiation data. As
such, a few prerequisites are required to run the code.

1. You will need a set of rules describing all of the possible cell
   differentiation events during each generation of the cell lines.
   Some examples are provided in the files till_rules.txt and tcell_rules.txt.

Here is the till_rules.txt example:
```
S D : 1 - P1
3*S : P1
```

Each line contains one rule, so there are two rules in this model. The symbols
on the left side of the colon (:) on each line indicate the possible states of
the model (in this case, S and D). The symbols to the right of the colon describe
the probability with which the rule occurs in mathematical terms. Let's take a
closer look at the first line:
```
S D : 1 - P1
```

This rule indicates that a of cell type, S, will transition to a cell type of D,
with a probability of 1 - P1, where P1 is a free parameter in range [0-1]. In this
model, the S cell type is a living cell, and the D cell type is a dead cell. Hence,
at each generation, the model assumes that cells die with probability 1 - P1.

Now let's look at the second rule:
```
3*S : P1
```

The 3\*S on the left side is a short-hand for three S symbols in a row (i.e. S S S).
Therefore, the rule is describing the situation in which a cell matching the left-most
symbol (S), will transition into two S cells. This will occur with probability P1. In
this model, the S cell type (a living cell) will divide into two living cells with
probability P1 (same free parameter as above).

You could also write this same rule in a different way:

```
S S S : P1
```

The code doesn't care whether you use the shorthand notation or not. What does matter
is that the left-most symbol is the "matching" part of the rule. Only a cell in the
state described by this symbol will cause the rule to "fire", and the rest of the
symbols (to the left of the colon, but to the right of the left-most symbol) describe
what kinds of cells are generated by the matching cell on the next generation. Overall,
the model above indicates that a cell will divide with probability P1, or will die
with probability 1 - P1. Feel free to make your models as complex as you like and note
that you can use more complex strings, so long as they qualify as sympy variable names.

However, you must ensure that the probabilities in your model are consistent. Take a look
at these two rules:

```
S S S : P1
S D : P2
```

These look very similar to the model above, but the use of -two- free parameters in the
model may lead to problems down the road. For example, once the model trees are built
by the code, it is important that you only consider fitting P1 and P2 such that P1+P2=1.
This is because there is an implicit assumption made by the code that if P1+P2 < 1, then
the remaining probability (1-P1-P2) is captured by the assumed rule (S S : 1-P1-P2) where
the cell simply lives to the next generation but doesn't divide. This is also implicitly
assumed by the D state (D D : 1), but this is what we wanted so this is OK. However, if
you specific P1+P2 > 1 in your final fitting, then the model breaks and the code will
provide -no warnings- about this. Instead, we say it here: you have been warned!

2. You will need an initially assumed population for your model. Examples of this are
   found in the files till_states.txt and tcell_states.txt.

The example file till_states.txt only contains the following information:

```
S
```

This means that we start with a single live cell. However, you may want more
complicated starting configurations. One example would be:

```
3\*S 4\*D
```

This would assume your model will start with a population of 3 live cells and
4 dead cells. Of course, having more cell states for living cells makes more sense
than including dead cells in the initial generation, but you can essentially set
this to whatever you like.

3. You will need to specifiy the number of generations you are assuming will match
   well with you experimental data. Note that the algorithms are exponential in
   the number of generations, so asking for too many will result in an expensive
   final tree compilation (the actual code generation is not too slow, but the
   compilation step described below -is-). Generally speaking then, if 7 generations
   is the maximum number of generations you will be considering, it is trivial to
   generate the trees and code for generations 2-6. Let's just start with 3 generations
   to make it fast:

```
$ mkdir tmp

$ cd tmp

$ ../build_tree.py 3 ../till_rules.txt ../till_states.txt
```

Assuming you started from the cell-differentiation directory which contains the code,
then the above commands will generate the C-code needed to instantiate 3 generations
of the Till Model. However, you still need to build the C-code now in order to
start fitting the probabilities in the model.

4. Utilize the generated Makefile to compile the C-code for the 3-generation Till
   Model tree from the last step:

```
$ make -j N
```

Be careful to set N here to the number of cores on your machine. For large models with
many rules, or for many generations, the code compilation can take a significant
amount of time (much longer than it took to generate the code) so be prepared to wait
awhile to see everything come together in this cases. To get a sense of how bad things
can get, try running for 8 generations instead of the 3 specified above and remember
that generation 9 will take exponentially more time than 8. Even simple models can get
complicated quite quickly, but experiments on cell lines are unlikely to be very long
as well since the cell type counts need to be obtained and are typically exponentially
growing as well. That's just life with cellular proliferation. :)

5. Finally, after you have build the generated C-code, you will see that several binaries
   have been created to help with parameter fitting: rule_probabilities,
   generation_XXX_summary, and generation_XXX_summary_unreduced

Each of these codes can help you fit the model to data, and each takes the same number
of arguments equivalent to the number of free parameters in your model. You can run
any of them without arguments, and they will report the order they are expecting for
these arguments. Plugging in a probability will generate an appropriate report. For
example:

```
$ ./rule_probabilities 0.3

Symbol  Rule    Probability
S       D       0.699999999999999956
S       2*S     0.299999999999999989
S       S       0
D       D       1
```

This reports all of the rules from your rules file, and the additional implicit rules
(if they were not defined explicitly in your rules file). Also, it reports what the
numeric probabilites of each rule are for the model, given your provided parameters.
More complicated models (like in tcell_rules.txt), utilize several free parameters and
this is a nice check to see if the rules and parameter values are consistent. If you
get probabilities outside of the [0-1] range from this utility, then you need to go
back and rethink your model -OR- make sure you are using parameter values that make
sense for the model you have built.

```
$ ./generation_001_summary 0.3

S D
0.699999999999999956 0.299999999999999989
0 0.699999999999999956
0.299999999999999989 0
```

This reports the proportion (probability) of all cell types given all possible
outcomes of the model in summary form. For example, starting with a single S, we
can see that one possible outcome is that 


