      SUBROUTINE WR_BLKMAT2(A,LROW,LCOL,NBLK,ISYM,ITRI)
C
C PRINT BLOCKED MATRIX WITH SYMMETRY (COLUMN SYMMETRY ORDERED)
C
      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
      include 'stdunit.h'
      include 'multd2h.h'
      DIMENSION A(*)
      DIMENSION LROW(NBLK),LCOL(NBLK)
C
      IBASE = 1
      WRITE(6,*) ' Blocked matrix '
      WRITE(6,*) '================'
      WRITE(6,*)
      DO IBLK = 1, NBLK
        JBLK = MULTD2H(IBLK,ISYM)
        IF (ITRI.NE.0 .AND. JBLK.LT.IBLK) CYCLE
        WRITE(luout,'(A,2I3)') ' Block ... ',JBLK,IBLK
        IF(ITRI.EQ.0.OR.IBLK.NE.JBLK) THEN          
          CALL WRTMAT2(A(IBASE),LROW(JBLK),LCOL(IBLK),
     &         LROW(JBLK),LCOL(IBLK) )
          IBASE = IBASE + LROW(JBLK)*LCOL(IBLK)
        ELSE
          CALL PRTRLT(A(IBASE),LROW(IBLK))              
          IBASE = IBASE + LROW(IBLK)*(LCOL(IBLK)+1)/2
        END IF

      END DO
      RETURN
      END
