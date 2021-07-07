@rem
@rem
@rem primary use of this is VCode building and debugging
@rem for people who do not want VStudio
@rem
@rem this build outputs to the same  
@rem folder as VStudio solution does
@rem so if changinf the name out/test.exe 
@rem be carefull not to hit the same
@rem exe names as VStudio solution uses 
@rem
@rem dbj_matmul_Release.exe and dbj_matmul_Debug.exe
@cls

@if [%1] == [clean] goto clean

:build
@rem obviously you know /Zi makes it a debug build
@rem and you know /MDd is dll runtime lib debug version
clang-cl /openmp /TC /MDd /Zi main.c /o out/test.exe
@goto exit

:clean
@del out\*.pdb
@del out\*.ilk
@del out\*.exe
@del r_and_d\*.pdb
@del r_and_d\*.ilk
@del r_and_d\*.exe
@goto exit

:exit