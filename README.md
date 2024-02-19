# halo

Control orbit to shadow a halo using RL

## Description

Apply corrections to a given orbit so that it stays close to a halo orbit. 
Corrections are applied to the velocities vx, vy. The problem to find
appropriate corrections is formulated as a Reinforcement Learning (RL) problem.
An agent learns from interacting with the environment (the RTBP equations),
with the goal of maximizing the reward / minimizing the loss (distance to the
halo).

## Getting Started

### Dependencies

* You will need a linux OS with the usual C/C++ developing
environment (compiler, libraries), as well as the python interpreter. Moreover,
you will need to install the GNU Scientific Library (GSL).

### Installing

* Download code to any directory and type `make` to build C programs.
* `make clean` to clean up C objects and binaries.

### Executing programs

* This project consists of several different programs, some written in C and
  some in python. They are run from the command line. They usually accept
  arguments on the command line, and output results to stdout. 
* Python programs: 
```
python halo.py
```
* C programs:
```
./fft orbit.dat > coefs.dat
./reward < coefs.dat
./correction
./shadowing
```

## Help

Please read source code for a detailed description of each program, including
arguments and results. Usually, running the program without arguments will
complain and explain the right usage, e.g.
```
./FT
Num of args incorrect. Usage: ./FT scaled_I phip
```

## Authors

Contributors names and contact info

Pablo Roldan
[@roldanpau](https://www.linkedin.com/in/pauroldan/)

<!---
## Version History

* 0.2
    * Various bug fixes and optimizations
    * See [commit change]() or See [release history]()
* 0.1
    * Initial Release
-->

## License

This project is licensed under the GNU GENERAL PUBLIC LICENSE - see the COPYING file for details

## Acknowledgments

<!---
Inspiration, code snippets, etc.
* [awesome-readme](https://github.com/matiassingers/awesome-readme)
* [PurpleBooth](https://gist.github.com/PurpleBooth/109311bb0361f32d87a2)
* [dbader](https://github.com/dbader/readme-template)
* [zenorocha](https://gist.github.com/zenorocha/4526327)
* [fvcproductions](https://gist.github.com/fvcproductions/1bfc2d4aecb01a834b46)
-->
