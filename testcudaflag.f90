PROGRAM TESTSFLAGS

#ifdef cublas
  write(*,*) 'defined'
#else
  write(*,*) 'not defined'
#endif

END PROGRAM
