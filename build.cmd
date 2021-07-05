@cls

@if [%1] == [clean] goto clean

:build
@rem obviously you know /Zi makes it a debug build
clang-cl /openmp /TC /Zi main.c /o out/test.exe
@goto exit

:clean
@del out\*.pdb
@del out\*.ilk
@del out\*.exe
@goto exit

:exit