# MGPT code


Authors:

Mario Alberto Rodriguez-Meza (ININ, Mexico), 
marioalberto.rodriguez@inin.gob.mx

Alejandro Aviles (Conacyt/ININ, Mexico),
alejandro.aviles.conacyt@inin.gob.mx, avilescervantes@gmail.com 

#


MGPT (Modified Gravity Perturbation Theory) code computes 2-point statistics for LCDM model and Hu-Sawicky-Starobinsky f(R) gravity. It is easily modifiable to other models. 

Specifically, it computes:

- SPT matter power spectrum

- SPT Lagrangian-biased tracers power spectrum

- CLPT matter correlation function

- CLPT Lagrangian-biased tracers correlation function

- A set of Q and R functions (extensions to those of Matsubara 2008b) from which other statistics, as leading bispectrum can be constructed. 

The code units are Mpc/h. 

The power spectrum convention is 
(2pi)^3 delta_D(k+k') P(k) = <delta(k) delta(k')> 


The code is divided in three pieces: I. MGPT-PS; II. MGPT-qfunctions; and III. MGPT-CLPT


## I. MGPT_PS

A C code that computes the matter and biased SPT power spectra. It computes the kernels in LPT by solving the set of differential equations of arXiv:1705.10719 and from them the functions Q(k) and R(k) of arXiv:1808.XXXXX 

DOWNLOAD:

Git clone

```
git clone https://github.com/cosmoinin/MGPT.git
```

or download it from http://www.github.com/cosmoinin/MGPT


Compile:

```
/MGPT/src$ make
```

Run: in the parent directory

```
/MGPT$ ./mgpt
```
For help:

```
/MGPT$ ./mgpt -help
```


In help you can see how to change parameters, for example:

```
/MGPT$ ./mgpt om=0.3 h=0.7 fR0=1.0e-6 
```

computes Hu-Sawicky f_R0 = -10^-6, and background cosmology h=0.7, Omega_m = 0.3.

Option screening=1 is the default with screenings, set to screening=0 if you want no screenings.

Alternatively you can run the code with a parameters file, for example

```
/MGPT$ ./mgpt parameters.in
```


The input of the code is the LCDM linear power spectrum extrapolated to present time with two columns 

| column  | function  |
| ------------: |:---------------| 
| #1            | k   (in h/Mpc)|
| #2            | P_L   (in (h/Mpc)^3)        |  




By default it is located in /MGPT/Input/psLCDM.in


The code gives two files as outputs: 

### output a) SPTPowerSpectrum.dat 

with structure


| column  | function  |
| ------------: |:---------------| 
| #1            | k   |
| #2            | PSL   (linear power spectrum in MG)        |  
| #3            | P22        |  
| #4            | P13        |  
| #5            | a10        |  
| #6            | a01        |  
| #7            | a20        |  
| #8            | a11        | 
| #9            | a02        | 

By default the range of wavenumbers is very large, because is necessary to compute the CLPT correlation function. To reduce the range, run the code as

```
/MGPT$ ./mgpt Nk=100 kmin=0.001 kmax=0.2
```

so the output will contain 100 k-points equally spaced in log intervals, from k=0.001 to kmax=0.2 h/Mpc



The 1-loop matter power spectrum is 

Ploop_m = PSL + P22 + P13 

and the 1-loop Lagrangian-biased tracers (X) power spectrum is 

Ploop_X = PSL + P22 + P13 + b1 * a10 + b2 * a01 + b1^2 * a20 +  b1 * b2 * a11 + b2^2 a02 

for b1 and b2 local Lagrangian biases.


The bias b_{01} from operator (\nabla^2 \delta) can be introduced by adding 

Ploop_X = Ploop_X - 2 (1 + b1)b_{01} k^2 PSL + b_{01}^2 k^4 PSL  

In this notation b_{01} = - b_{\nabla^2}


### output b) MGPT/CLPT/kfunctionsT.dat

kfunctionsT.dat file contains all the Q(k) and R(k) functions. It is the input of the code MGPT-qfunctions.



| column  | function  |
| ------------: |:---------------| 
| #1            | k   |
| #2            | Q1          |  
| #3            | Q2        |  
| #4            | Q3        |  
| #5            | Q5        |  
| #6            | Q7        |  
| #7            | Q8        |  
| #8            | Q9        | 
| #9            | Q11        | 
| #10           | Q12        | 
| #11           | Q13        | 
| #12           | QI        | 
| #13           | R1        | 
| #14           | R2        | 
| #15           | R1plus2        | 
| #16           | RI   | 
| #17           | Dpk  (D+(k): linear growth function as a function of k)         | 
| #18           | PSL  (Linear power spectrum in MG)       | 

## II. MGPT_qfunctions

CLPT/MGPT_qfunctions.nb is a Mathematica notebook that post-process the file kfunctionsT.dat to obtain a set of q-functions
that are the building blocks of the CLPT correlation function. This code is independent of the gravitational or dark energy model, just fed it with the appropriate file structure. 

The input is the file kfunctionsT.dat obtained with MGPT-PS (or by other code of your preference)

The output is the file qfunctionsT.dat, with columns structure






| column  | function  |
| ------------: |:---------------| 
| #1            | q      |
| #2            | X_linear          |  
| #3            | Y_linear        |  
| #4            | X_loop        |  
| #5            | Y_loop        |  
| #6            | V        |  
| #7            | T        |  
| #8            | X10        | 
| #9            | Y10        | 
| #10           | U_linear        | 
| #11           | U_loop        | 
| #12           | U11        | 
| #13           | U20        | 
| #14           | xi_linear    (linear correlation function)     | 
| #15           | nabla2_xi_linear   (The laplacian of the linear correlation function)      | 
| #16           | nabla4_xi_linear   (The \nabla^4 of the linear correlation function)   | 



## III. MGPT_CLPT


CLPT/MGPT_CLPT.nb is Mathematica notebook post-process the ouput of MGPT_qfunctions.nb to compute the Zel'dovich approximation, and CLPT correlation function for both matter and biased tracers. 

Input: qfunctions.dat obtained from MGPT_qfunctions.nb (or by other means)

Output: CorrelationFunction.dat, with columns structure





| column  | function  |
| ------------: |:---------------| 
| #1            | r     |
| #2            |  xi_L     (linear correlation function)          |  
| #3            |  xi_ZA    (Zel'dovich approximation correlation function)        |  
| #4            | xi_A       |  
| #5            | xi_W        |  
| #6            | xi10_linear        |  
| #7            | xi10_loop       |  
| #8            | xi20_linear       | 
| #9            | xi20_loop       | 
| #10           | xi01        | 
| #11           | xi02        | 
| #12           | xi11        | 
| #13           | xi_nabla2        | 
| #14           | xi_nabla4    | 


The matter CLPT correlation function is given by 

xi_CLPTm = xi_ZA + xi_A + xi_W

The biased tracers CLPT correlation function is

xi_CLPT_X = xi_CLPTm + b1 * (xi10_linear + xi10_loop)  + b2 * xi01 + b1^2 * (xi20_linear + xi20_loop) +  b1 * b2 * xi11 + b2^2 xi02 

for b1 and b2 local Lagrangian biases.



The bias b_{01} from operator (\nabla^2 \delta) can be introduced by adding 

xi_CLPT_X = xi_CLPT_X  + 2 (1 + b1)b_{01} xi_nabla2 + b_{01}^2 xi_nabla4.


## References

If you use this code please cite the following papers:

1. Alejandro Aviles and Jorge-Luis Cervantes-Cota [Phys. Rev. D 96, 123526 (2017)] https://arxiv.org/abs/1705.10719

2. Alejandro Aviles, Mario Alberto Rodriguez-Meza, Josue De-Santiago, and Jorge-Luis Cervantes-Cota [arXiv:1809.XXXXX]


For the theory in LCDM  see the following papers:

1. https://arxiv.org/abs/0807.1733
2. https://arxiv.org/abs/1209.0780
3. https://arxiv.org/abs/1506.05264
4. https://arxiv.org/abs/1805.05304


