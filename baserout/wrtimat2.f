      SUBROUTINE WRTIMAT2(A,NROW,NCOL,NMROW,NMCOL)
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      include 'stdunit.h'
      INTEGER A
      DIMENSION A(NMROW,NMCOL)
C
      ICOLMX=5
      ICOLL=0
      ICOLH=0
      DO WHILE (ICOLH.NE.NCOL)
        ICOLL = ICOLH+1
        ICOLH = MIN(ICOLL-1+ICOLMX,NCOL)
        WRITE(luout,1000) (J,J=ICOLL,ICOLH)
        write(luout,'(x,4("-"),"+",55("-"))')
        DO I=1,NROW
          WRITE(luout,1010) I,(A(I,J),J=ICOLL,ICOLH)
        END DO
      END DO

      RETURN
 1000 FORMAT(4X,1X,"|",5(1X,2X,I6,2X))
 1010 FORMAT(1X,I3,1X,"|",5(1X,I10))
      END
