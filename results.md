# DBJ MATMUL BENCHMARKING

## Windows 10 PRO

clang-cl aka clang 10.0.0 release build

```
Microsoft Windows [Version 10.0.19042.1052]
(c) Microsoft Corporation. All rights reserved.

C:\Users\ddbj_matmul_Release.exe atmul\out>dbj

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