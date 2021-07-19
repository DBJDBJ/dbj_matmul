<h1> DBJ*MATMUL</h1>
<h3> Benchmarking </h3>
<h3>&nbsp;</h3>

#### (c) 2021 by dbj at dbj dot org -- https://dbj.org/license_dbj/

<h3>&nbsp;</h3>

 This is benchmarking of a collection of matrix multiplication algorithms.
 Algorithms are kept as simple as possible. No structs are passed as arguments.
 No "clever" and/or "generic" matrix macros are used.

 Different compliers multiplied with different platforms multiplied by selection
 of data types  yield a complex picture of benchmarking results.

 > Here is strong hint for you: In general situation where no special treatments are possible, the simplest algorithm is the fastest. 

 Keep in mind compiler has the easiest job optimizing the simplest code.

 Use `dbj_matmul.c` to recompile and re measure whenever selecting
 the right matrix multiplication algorithm for your project needing one. And that is best done on the hardware where your app will be running. 

 Results of benchmarks are almost always surprising. And it is not good to be surprised after you have released something to the customer.



 