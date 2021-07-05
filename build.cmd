@cls

@if [%1] == [clean] goto clean

:build
clang-cl /openmp /TC /Zi matmul_1d_seq.c /o out/test.exe
@goto exit

:clean
@del out\*.pdb
@del out\*.ilk
@del out\*.exe
@goto exit

:exit