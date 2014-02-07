# autoCorrelation

## Purpose

Computes the autocorrelation function of a three-dimensional array signal for multiple particles.

## Background

Please refer to manual.pdf for detailed instructions and a discussion of the mathematical background.

## Author

Written by Maximilien Levesque, École Normale Supérieure, Paris, France.

## Important contributors

* Vincent Dahirel, UPMC, PECSA, Paris, France, for rewritting the brute force algorithm in a faster and more readable way, Oct. 2013.
* Marie Jardat, UPMC, PECSA, Paris, France, for discussions, extensive testing and bug reports. Sept. 2013.
* Xudong Zhao, UPMC, PECSA, Paris, France, for providing test cases, for discussions and Apple specific bug reports. Sept. 2013

## Changelog

* Sept. 2013:   Basic version of the code, with bruteforce algorithm only.
* Oct. 2013:    Various improvements of the code and of the bruteforce algorithm.
* Feb. 2014:    Introduction of the Fourier space algorithm, inducing runtimes divided by 100.

## How to make it work

### Download it

The best and easiest way is to use git:
```
    git clone https://github.com/maxlevesque/autoCorrelation  
```
will create a folder called autoCorrelation and download the latest version of the sources automatically.

You may also go to https://github.com/maxlevesque/autoCorrelation/archive/master.zip from your web brower.

Or use `get`, `wget`, ... or whatever seems good to you. I repeat: the best and easiest way to download and stay up-to-date is to use `git`.

### Compile it

The Fourier space algorithm requires the `FFTW3` library. Please visite their [website](www.fftw.org) for more informations.

You need to install `scons`, which is a newer and better and easier equivalent to the much complicated `make`.
Please visit [the scons website](www.scons.org) for download, or you should better use the repositories of your own distribution:  
Under Ubuntu: `sudo apt-get install scons`  
Under Fedora: `sudo yum install scons`  
and so on for other distros.

Once you're in the directory where you downloaded `autoCorrelation`, just type  
```
$ scons
```  
In a nutshell, scons will look at the file called `SConstruct` I built for you and compile everything smartly.

Of course, you may compile the simple fortran file(s) by yourself if you prefer the complicated ways.

## How to use `autoCorrelation`

### Inputs
All you need is a file, whatever its name, containing all your velocities (or forces, or ...) in a format  
``` 
x(1,t) y(1,t) z(1,t)  
x(2,t) y(2,t) z(2,t)  
x(i,t) y(i,t) z(i,t)  
...    
x(Nat,t) y(Nat,t) z(Nat,t)    
x(1,t+1) y(1,t+1) z(1,t+1)  
x(2,t+1) y(2,t+1) z(2,t+1)  
x(i,t+1) y(i,t+1) z(i,t+1)  
...  
x(Nat,t+1) y(Nat,t+1) z(Nat,t+1)  
...  
...  
x(Nat,t+Nstep) y(Nat,t+Nstep) z(Nat,t+Nstep)  
```  

where `Nat` is the total number of atoms you have in your supercell,  
and `Nstep` is the number of timesteps in your trajectory.

In other words, you print the coordinates of all sites for a given timestep, then for the next one, etc.  
You do not print blank lines anywhere. You can add as many spaces you want between the columns.

At the end, you *must* have 3 columns and `Nat x Nstep` rows.

That's all you need.

Note that if you have one-dimensional (two-dimensional) vector at each time step, you'll have to fill the last 2 (1) columns with zero.

## Execution

The executable is waiting for arguments:  
1. `Nat`, defined above  
2. `a`, an integer with value `1` or `2`, for bruteforce or fourierspace algorithm.
2. `filename` of the trajectory in the format discussed above.  
  
So, you have to execute:  
`$ autoCorrelation 800 2 ./analysis/velocities.out`

## Outputs

Execution of `autoCorrelation` will result in a single file: `acf.out`.  
This file has a very simple ASCII format:
``` 
1 1.12931233  
2 1.24034134  
3 1.5908O123  
...  
Nstep 123987.12383  
``` 
where the first column is the timestep, and the second column is the autocorrelation function.  

### Units
Both timesteps and autocorrelation functions are in units of your simulation.  


## Disclaimer

This program has been thoroughly tested and validated.
I would not share a tool that is not production ready. I am very confident in this program.
*Nevertheless*, this program may contain bugs or restrictions that would lead to unexpected results.
Please, read carefully this readme file, and test it thoroughly on data for which you already know the results.

Comments, bug reports or even thanks ;) are welcome!
