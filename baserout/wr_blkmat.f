      SUBROUTINE WR_BLKMAT(A,LROW,LCOL,NBLK,ISYM)
C
C PRINT BLOCKED MATRIX
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'stdunit.h'
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)
C
      IBASE = 1
      WRITE(6,*) ' Blocked matrix '
      WRITE(6,*) '================'
      WRITE(6,*)
      DO 100 IBLK = 1, NBLK
        WRITE(lulog,'(A,I3)') ' Block ... ',IBLK
        IF(ISYM.EQ.0) THEN
          IF(IBLK .NE. 1 ) IBASE = IBASE + LROW(IBLK-1)*LCOL(IBLK-1)
          CALL WRTMAT2(A(IBASE),LROW(IBLK),LCOL(IBLK),
     &         LROW(IBLK),LCOL(IBLK) )
        ELSE
          IF(IBLK .NE. 1 ) 
     &         IBASE = IBASE + LROW(IBLK-1)*(LCOL(IBLK-1)+1)/2
c          CALL PRSYM(A(IBASE),LROW(IBLK))              
          CALL PRTRLT(A(IBASE),LROW(IBLK))              
        END IF

  100 CONTINUE
      RETURN
      END
