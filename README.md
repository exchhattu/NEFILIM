## NEFILIM

NEFILIM (http://www.riken.jp/zhangiru/software.html) is iterative fragment-assembly 
de novo protein structure program. It generates a new set of fragments from the 
lowest energy de novo models. Subsequently, the fragments are further used in prediction. 
NEFILIM performed better in comparison. NEFILIM was also tested in CASP11
(http://predictioncenter.org/casp11/index.cgi).

## Dependencies
* Rosetta3.2 (https://www.rosettacommons.org/software) 
* g++  
* cppunit 

After installing the dependencies, their installed paths should replace ### ROSETTA_INSTALL_PATH
and ### CPPUNIT_PATH in SConstruct file to assign rosetta_home and cppunit_top variables. 

## Example 
$ ./example.sh

## Reference
* R. Shrestha, K. Y. Zhang, Improving fragment quality for de novo structure prediction. 
Proteins: Structure, Function, and Bioinformatics, 2014, 82, 2240-2252  
