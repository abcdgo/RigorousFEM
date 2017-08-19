# RigorousFEM


### Installation

#### CAPD
To use the code from the repository one must install [CAPD library](http://capd.ii.uj.edu.pl/).
While currently 4th version of the library is publicly available the software was written with unpublished early 4th version which in fact is something between v3 and v4. Therefore it is included into repository as *capd3.5.tar.gz* file.

##### CAPD Installation
Please, follow full installation procedure which can be found on [CAPD](http://capd.sourceforge.net/capdDynSys/docs/html/capd_compilation.html) web page.

**Warning:** because of a change in GCC compiler the CAPD provided within this repository does not compile under GCC 4.9. Use GCC 4.8 or earlier. If you are C++ savvy you can use GCC 4.9 and fix compilation errors by yourself.

###### Quick guide
As stated before it is strongly recommended to read and follow [CAPD](http://capd.sourceforge.net/capdDynSys/docs/html/capd_compilation.html) installation manual.

TL;DR version:

Download repository and unarchive the *capd3.5.tar.gz* file.
```
$ git clone -b v1.0 --depth 1 https://github.com/abcdgo/RigorousFEM.git
$ cd RigorousFEM
$ tar -zxvf capd3.5.tar.gz
```
In the same place where the unarchived folder is do:
```
$ mkdir capd_bin
$ cd capd_bin
$ ../capd/configure --prefix $PWD -with-filib=check -with-mpfr=check -without-gui
$ make -j8
$ make install
```
where in place of $PWD you should put the absolute path to the current directory

#### Solver
Except CAPD part, the rest of the software does not need any further installation. Only compiling the code is required. In a *Proof* subfolder there are *candidate.cpp* and *inclusion.cpp* files which do the proofing procedure. To help automate the process *makeProof.sh* script is included (see Usage section). In short, download the *Solver* folder and run the script.

### Usage
*makeProof.sh* - is a bash script which automates the proofing procedure. It compiles the code and executes it, prepares proof summary and saves computation time data. The script takes 4 arguments and all of them are required.
```
$ ./makeProof.sh $M $m $NI $HULL,
```
where M - number of modes, m - number of non dissipative modes, NI - is a viscous parameter of the equation, HULL - size of a hull. Example call:
```
$ ./makeProof.sh 100 13 6 1e-3
```
**Prior first use, lines 32, 46 and 51 of the script must be altered to reflect CAPD installation paths.**

###### Script characteristic
What the script does is:
- compiles *Proof/candidate.cpp* code then puts binary file to *bin* folder and executes it
 - output of the candidate program is directed to *Scripts/candidateToRigorous.py* script and prepared for rigorous computation
  - all output is saved in created *experiments* folder
- executes rigorous computation *Proof/inclusion.cpp* (same behavior as for *candidate.cpp*)
 - checks if the starting disk includes the result disk if not then disk refinement is proceeded and rigorous computation with new disk is executed via *Proof/iter2.cpp* (it is the same file as inclusion but for bash script reasons it was needed to split it that way, bash script should be refactored)
 - the checks does not have stop condition, so it should be human supervised and break if necessary
 - each iteration result is processed by *Scripts/checkProof.py* script

### Sample output
Folder *Sample results* contains few sample proof outputs produced by the *Solver*

### Citation
A BibTeX entry for LaTeX users is:
```
@article{,
	Author = {Pa{\l}ka, Wies{\l}aw},
	Journal = {GitHub repository},
	Publisher = {GitHub},
	Title = {RigorousFEM},
	Url = {https://github.com/abcdgo/RigorousFEM},
	Year = {2016}}
```
### License
For details see the LICENSE file.
