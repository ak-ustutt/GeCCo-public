      FUNCTION XFAC(N)
*
* N !  as double precision real
*
      IMPLICIT REAL*8(A-H,O-Z)
      IF( N .LT. 0 ) THEN
       IFAC = 0
       WRITE(6,*) ' WARNING FACULTY OF NEGATIVE NUMBER SET TO ZERO '
      ELSE
C
       XFACN = 1.0D0
       DO 100 K = 2,N
        XFACN = XFACN * DFLOAT(K)
  100  CONTINUE
       XFAC = XFACN
      END IF
C
      RETURN
      END
