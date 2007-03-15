      FUNCTION IFAC(N)
C
C N !
C
      IF( N .LT. 0 ) THEN
       IFAC = 0
       WRITE(6,*) ' WARNING FACULTY OF NEGATIVE NUMBER SET TO ZERO '
      ELSE
C
       IFACN = 1
       DO 100 K = 2,N
        IFACN = IFACN * K
  100  CONTINUE
       IFAC = IFACN
      END IF
C
      RETURN
      END
