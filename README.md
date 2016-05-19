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

The code is designed to aid parameter estimation when fitting cell
differentiation models to experimental cell differentiation data. As
such, a few prerequisites are required to run the code.

1. You will need a set of rules describing all of the possible cell
   differentiation events during each cell generation. Some examples
   are provided in the files <till_rules.txt> and <tcell_rules.txt>

