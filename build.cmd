@rem
@rem
@rem primary use of this is using clang-cl 
@rem
@rem
@cls

@if [%1] == [clean] goto clean

:build
@rem obviously you know /Zi makes it a debug build
@rem and you know /MDd is dll runtime lib debug version
@rem clang-cl /TC /MDd /Zi godbolt.c /o out/dbj_matmul_clang_debug.exe
@clang-cl /TC /MD godbolt.c /o out/dbj_matmul_clang_release.exe
@goto exit

:clean
@del out\*.pdb
@del out\*.ilk
@del out\*.exe
@del *.pdb
@del *.ilk
@del *.exe
@del r_and_d\*.pdb
@del r_and_d\*.ilk
@del r_and_d\*.exe
@goto exit

:exit