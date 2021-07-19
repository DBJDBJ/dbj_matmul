
# DBJ-MATMUL

## Version 0.5

### Benchmarking various Matrix Multiplications.

(c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj


## Building

Both `dbj_omp.c` and `dbj_mamtmul.c` can be copy-pasted and used in Godbolt. Of course with much smaller data sets  vs the local usage on your workstation.

### DEVENV

| Part     | Description                                                                         |
| -------- | ----------------------------------------------------------------------------------- |
| OS       | Windows 10                                                                          |
|          | Linux has not be tried yet, but I see no reason for that not to work                |
| IDE      | I use both Visual Studio 2019 and VS Code. Mainly because of VS debugger            |
| Language | C99 or better                                                                       |
| Compiler | clang 11.0, from inside Visual Studio 2019. Namely: `clang-cl.exe`, 64 bit version. |
VS Code does not recognise clang-cl.exe (as of 2021 July). Only sometimes. 
|          | TDM-GCC-64 is used when CTRL+SHIFT+B compiling from inside the Code. You will be more than capable to fix `.vscode/tasks` to match your devenv 
|          | cl.exe can not be used. That is still not full ISO C compiler

### Open MP 

`dbj_omp.c` purpose is to measure OMP vs NO OMP. With very large datasets OMP show a bit better performance. But in that case one is wondering if and where she will need so large datasets, and why not using GPU in that case. For normal sizes datasets gcc/clang optimised code is faster. With added simplicity of not using OMP. '-fopenmp'

| compiler      | Comment                                                                                                             |
| ------------- | ------------------------------------------------------------------------------------------------------------------- |
| gcc or clang  | Compile with: `-fopenmp -O3 -lm` . As already set in `.vscode/tasks.json`                                           |
|               | that also applies to Godbolt                                                                                        |
| Visual Studio | Go to project properties; under "Configuration Properties --> C/C++ --> Language --> Open MP Support" select "Yes". |

### What is this other stuff?

I use the simple and good [`ubench.h`](https://github.com/sheredom/ubench.h) and [`utest.h`](https://github.com/sheredom/utest.h) . 

Instead of using them as git submodules, I have copied them over to ubench and utest folders in here. Why?

The current situation is I found few little bugs and fixed them but only locally to this repo. I will do the proper PR, ASAP. Simply there was no enough time.

The stuff I useth but discardedth is in the `r_and_d` folder. Simply I have focused on "job at hand" and that was benchmarking matrix multiplications and omp situation. 

Other folders are I hope descriptively named so that one can immediately understand what is in them.

`build_time_stamp.inc` is artefact of my simple VS project timestamping. As [described in here](https://dbj.org/visual-studio-project-simple-time-stamp-for-your-c-c-code/).

## Caveat Emptor

I do use result of these benchmarks myself.

You use the code in here 100% as your own responsibility. I provide no guarantees or warranties of any kind.

Any question or critique please do use the GitHub issues tab.

Many thanks...

