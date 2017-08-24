      PROGRAM FTEST
      
! Tests the Fortran-C Interface

      REAL value(2)
      
      value(1)=2.0
      value(2)=2.0

      WRITE(6,*)
      WRITE(6,*)
     &  '--------------------------------------------------------------'      
      WRITE(6,*) 'This program shows the Fortran/C interface'
      WRITE(6,*) 'required for building the Unified Model.'
      WRITE(6,*) 
     &  'Please refer to Documentation Paper X4 for more details'
      WRITE(6,*)
      
      CALL CTEST(value)

      WRITE(6,*)
      WRITE(6,*)
     &  '--------------------------------------------------------------'
           
      END
