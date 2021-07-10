<h1> DBJ*MATMUL</h1>
<h3> Benchmarking </h3>
<h3>&nbsp;</h3>

#### (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj/

<h3>&nbsp;</h3>

 This is benchmarking of a collection of matrix multiplication algorithms.
 Algorithms are kept as simple as possible. No structs are passed as arguments.
 No "clever" "generic" matrix macros are used.

 Different compliers multiplied with different platforms multiplied by selection
 of data types  yield a complex picture of benchmarking results.

 > Here is strong hint for you: The simplest algorithm is the fastest. 
 Keep in mind compiler has the easiest job optimizing the simplest code.

 Use this file to recompile and re measure whenever selecting
 the right matrix multiplication algorithm

 [Godbolt code is here](https://godbolt.org/z/joMo4Tn3T)

 
## Windows 10 PRO

clang-cl aka clang 10.0.0 release build

```
Microsoft Windows [Version 10.0.19042.1052]
```
```
Testing various matrix multiplication algorithms
All matrices are square, side size is: 99       
[==========] Running 6 benchmarks.
[       OK ] simple (mean 647.468us)
[       OK ] not_faster (mean 752.153us)
[       OK ] snazzy (mean 645.381us)
[       OK ] simple_mat_mul_0 (mean 638.535us)
[       OK ] SSE_simpler_matmul (mean 200.209us)
[       OK ] SSE_better_matmul (mean 192.055us)
[==========] 6 benchmarks ran.
[  PASSED  ] 6 benchmarks.
```

It seems SSE has the role for large(r)matrices. 


## GODBOLT

[Godbolt](https://org/z/oYenYbxh1) 

CLANG 12.0.0 -s -lm -O3

```
Testing various matrix multiplication algorithms
All matrices are square, side size is: 99
[==========] Running 6 benchmarks.
[       OK ] simple (mean 974.103us)
[       OK ] not_faster (mean 1.102ms)
[       OK ] snazzy (mean 906.416us)
[       OK ] simple_mat_mul_0 (mean 950.146us)
[       OK ] SSE_simpler_matmul (mean 209.751us)
[       OK ] SSE_better_matmul (mean 215.459us)
[==========] 6 benchmarks ran.
[  PASSED  ] 6 benchmarks.
```

[Godbolt](https://org/z/joMo4Tn3T)

GCC 11.1 -s -lm -O3

```
Testing various matrix multiplication algorithms
All matrices are square, side size is: 99
[==========] Running 6 benchmarks.
[       OK ] simple (mean 902.017us)
[       OK ] not_faster (mean 959.547us)
[       OK ] snazzy (mean 907.815us)
[       OK ] simple_mat_mul_0 (mean 900.294us)
[       OK ] SSE_simpler_matmul (mean 165.025us)
[       OK ] SSE_better_matmul (mean 162.341us)
[==========] 6 benchmarks ran.
[  PASSED  ] 6 benchmarks.
```