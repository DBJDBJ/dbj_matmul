
# DBJ-MATMUL

Parallel and Sequential versions of Matrix Multiplication of non-square matrices.

Inspiration: https://github.com/ivanbgd/Matrix-Multiplication-MatMul-C
The rest: (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj


## Building

### DEVENV

| Part | Description
|---|----------
OS | Windows 10
IDE | Visual Studio 2019
Language | C99 or better
Compiler | clang 11.0, from inside Visual Studio 2019. Namely: `clang-cl.exe`, 64 bit version.
VS Code does not recognise clang-cl (as of 2021 July) | so one can not just hit F5 and debug. Thus we use build.cmd

### Open MP 

| compiler | comment
|----------|----------
| gcc or clang | Compile with: `-fopenmp -O2 file_name.c -o file_name -lm`
| Visual Studio | Go to project properties; under "Configuration Properties --> C/C++ --> Language --> Open MP Support" select "Yes".

### Third party

| Part | Comment
|------|-----------
|`macro.h` and `macro-fundamental.h` | nicked from `systemd` GitHub repository. Refactored.
