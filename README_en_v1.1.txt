
####################
Progressive BKZ Library
Version 1.1
Released date 2016/06/29

Contact email address:
pbkz-info@ml.nict.go.jp

Contact postal address: 
4-2-1, Nukui-Kitamachi, Koganei, Tokyo, 184-8795, Japan.

The latest version is available at http://www2.nict.go.jp/security/pbkzcode/index.html

####################
    INTRODUCTION
####################

This progressive BKZ library is an implementation of the algorithm proposed by
Y. Aono, Y. Wang, T. Hayashi and T. Takagi, 
in "Improved Progressive BKZ Algorithms and their Precise Cost Estimation by Sharp Simulator",
published in Eurocrypt 2016. [REFERENCE 1] 
The full-version is available at https://eprint.iacr.org/2016/146.
The experimental results given in this paper also can be verified using this progressive BKZ library.

Author and copyright: 
Security Fundamentals Laboratory, Cybersecurity Research Institute in
National Institute of Information and Communications Technology.

The core part of BKZ (pbkzmain.cpp) and ENUM (vectorenumeration.cpp) are 
modified from the BKZ subroutine in NTL library.
We follow the copyright notices in "doc/copying.txt" in Shoup's NTL library 9.8.0:

----------------------------------------------------------------------
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version. 

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
----------------------------------------------------------------------

This is an open-source library distributed under the terms of the GNU General Public License version 3,
without any warranty.


####################
    REQUIREMENTS
####################

Please ensure that you have already installed the following libraries,
which are required to support our project.
We recommend you to use versions no lower than the given ones.
Please see their homepages to get more information of the latest version.

gcc-5.3.0 (https://gcc.gnu.org/)
NTL-9.7.0 (http://www.shoup.net/ntl/)
gmp-6.1.0 (https://gmplib.org/)
gsl-1.16 (http://www.gnu.org/software/gsl/)
boost 1.58.0 (http://www.boost.org/)


####################
    INSTALLATION
####################

###
 A
###

Download and extract the file Ågpbkzlib-xxx.zipÅh.

###
 B
###

Some assistant subroutines are directly picked from open source libraries.
Please copy them to /include/external if there is any absence.

B1. Bases generator from SVP Challenge:
http://www.latticechallenge.org/svp-challenge/download/generator.zip
In the generator folder, Ågtools.hÅh and Åggenerate_random.cppÅh are necessary.

B2. Bases generator from Ideal Lattice Challenge:
http://www.latticechallenge.org/ideallattice-challenge/download/generator.zip
In the generator folder, Ågideal.hÅh and Åggenerate_ideal.cppÅh are necessary.

B3. Subroutines from NTL Library:
http://www.shoup.net/ntl/
Ågsrc/LLL_QP.cÅh and Ågsrc/LLL_RR.cÅh are necessary.

NOTE: Here in B3 we recommend using NTL 9.6.0.

Totally the 6 supplements above are required.

###
 C
###

Configure the bkz.conf, to set a path for storing the cache file of
pruning functions and simulations.

###
 D
###

run %make in the extracted directory.
Note that according to the compilers,
some warning messages may be imported during the compiling process.
Please ignore them if the file Åga.outÅh is output successfully. 



####################
    VERIFICATION
####################

Here are several examples for verifying the facts narrated in [REFERENCE 1].
Some shortening options as follows are commonly used.

if = input file
of = output file
sf = strategy file (for progressing blocksize in BKZ)
lf = logfile
si = LLLStatusInterval (seconds between status reports)
nt = number of threads to be used in BKZ and ENUM

There are also some other options as Åg-genbasisÅh, Åg-lllÅh, Åh-bkzÅh, Åg-enumÅh, Åg-genstrategyÅh, Åg-pfuncÅh,
Åh-randomizebasisÅh, Åg-withenumÅh, Åg-type (svp/ideal/unimodular)Åh and so on.
We will give some usages of these options in the examples below.
Also please reference Ågmain.cppÅh for more details.

Note that the bases generating program works in a single thread.


###
Ex.1 Verify the contents in section 6.2 by Ideal Lattice Challenge.
###

#1.1# Generate an LLL-reduced basis of Ideal Lattice Challenge basis with index=421 and seed=0.
%./a.out -genbasis -type ideal -index 421 -seed 0 -of isvpc421.txt -lll -si 10 -of lll_i421.txt 

#1.2# Generate the blocksize progressing strategy for approximate Ideal Lattice Challenge.
%./a.out -if lll_i421.txt -genstrategy -type shortvec -det 420 -of str421approx

#1.3# Run progressive BKZ to solve approximate Ideal Lattice Challenge.
%./a.out -if lll_i421.txt -bkz -sf str421approx -of bkz_i421.txt

Note that the generating step 1.2 takes some time.
The output file Ågstr421approxÅh includes the BKZ blocksize strategy. 


###
Ex.2 Verify the contents in section 8 by SVP Challenge.
###

#2.1# Generate an LLL-reduced basis of SVP Challenge basis with dim=100 and seed=0.
%./a.out -genbasis -type svp -dim 100 -seed 0 -of svpc100.txt -lll -of lll100.txt

#2.2# Generate BKZ blocksize strategy to get the first vector of output basis shorter than 1.05GH with ENUM.
      Here GH means Gaussian Heuristic and ENUM means Enumeration algorithm using modified extreme pruning technic. 
      For more details please see [REFERENCE 1].
%./a.out -if lll100.txt -genstrategy -type shortvec -a 1.05 -withenum -of str100_1.05gh -nt 6 -stopatfound -lf svp100.log -ll 1

Note that in step 2.2, it outputs the BKZ blocksize strategy file: str100_1.05gh 
and the bash script file: str100_1.05gh.sh

#2.3# Change the file system modes of the output shell script 
	 and run it to find a vector shorter than 1.05GH.
%chmod +x str100_1.05gh.sh
%./str100_1.05gh.sh

You can check the BKZ reduced basis from Ågstr100_1.05gh.bkz.log.**Åh,
get the details of ENUM search process from Ågstr100_1.05gh.enum.log.**Åh
and ENUM search results from Ågsvp100.logÅh.


####################
    MORE USAGES
####################

###
One-by-one progressive BKZ
###
Run progressive BKZ by shifting the blocksize with one-by-one strategy.
%./a.out -if (inputfile) -bkz -eb (endblocksize) -of (outputfile)

###
Generating pruning function
###
Generate an optimized pruning function for extreme pruning strategy used in ENUM subroutine.
%./a.out -if (input basis file) -pfunc -pf (input coefficients file if necessary) -p (probability) -a (radius) -vl (output level) -optimize (optimizing time in seconds) -of (outputfile)

Example: ./a.out -if lll100.txt -pfunc -pf pfunc100.txt -p 0.0005 -a 1 -vl 3 -optimize 30 -of pfunc100.txt

###
Lattice vector enumeration
###
Run ENUM search algorithm to find short vectors with a given radius parameter alpha and success probability.
%./a.out -if (input basis file) -enum -alpha (radius) -prob (probability) -of (output file) -lf (logfile) -nt (number of threads) -optimize (time in seconds to optimize pruning function) -preprocess (time in seconds to preprocess the input basis)

To find close vectors given by a vector file (NTL format), use "-enum -cvp vec.txt" instead of "-enum" option.


####################
   REPORTING BUGS
####################

If you find a concrete bug please report it to pbkz-info@ml.nict.go.jp.

All bug reports should include:

       The version number of progressive BKZ, and where you obtained it.
       The hardware, the operating system.
	  The version numbers of libraries listed in REQUIREMENTS above, if there is any difference.
       A description of the bug behavior.
       A short program which can reproduce the bug.

It is useful if you can check whether the same problem occurs after you update your libraries.
Thank you.

Any errors or omissions in the manual can also be reported to 
the same address as pbkz-info@ml.nict.go.jp.


#############
   HISTORY
#############


2016/06/29 v1.1

 - Fixed bugs in computing pruning functions in higher dimensions (>150)
 - Added new command for simulating Gram-Schmidt Lengths 
 - Added new command for lattice vector enumeration to find close vectors

2016/05/02 v1.0


########################
   MISCELLANEOUS INFO.
########################

Subroutine to generate Unimodular matrices (lattice/gen_uni_mat.cpp) 
is recoded from the open source code Sage-7.1, where the function is Ågrandom_unimodular_matrix()Åh.

This library is not thread-safe.









