# INFO-F-308

## Purpose

This git repository intends to summarize the work on the article *Uncovering disease-disease relationships
through the incomplete intneractome* from J. Menche et al.

As a 3rd year student in computer sciences at Belgian university *Université Libre de Bruxelles*, one of the
projects I have been assigned to is the following: work on the previously cited article with Pr. Tom Lenaerts
in order to reproduce the obtained results and upgrade it to the new versions of the used databases.

## Hierarchy

The repository is for now organized as follows:

+ `source/` contains all the python source both provided with the article and produced by myself.
+ `data/` contains the provided data and their newer versions.
+ `report/` contains the report on the produced content and the conclusions of it.

## TODO

### Bioinformatics

+ Have the provided source code to evolve from Python2 to Python3.  [**DONE**]
+ Factorize scripts in order to avoid duplication.  [**DONE**]
+ Run analysis source code on both *old* data (provided by the authors) and newer versions.
+ Get last version of required datasets.
+ Write report.

### Mathematics

+ Check if current proof of $\Lambda_{k,\alpha}^m(V, \cdot)$ isomorphism (and thus cardinalty) is correct, and check
on several examples to see if result can be used in Python3 scripts.
+ If not correct, then correct it, and if correct, then re-arrange report and improve rigorous treatments.
+ (If result is correct, then compare efficiency with curent simulation model).
