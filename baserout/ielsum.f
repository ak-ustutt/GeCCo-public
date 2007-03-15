      FUNCTION IELSUM(IVEC,NELMNT)
*
* Sum elements of integer vector IVEC
*
      DIMENSION IVEC(*)
*
      ISUM = 0
      DO 100 IELMNT = 1, NELMNT
        ISUM = ISUM + IVEC(IELMNT)
  100 CONTINUE
*
      IELSUM = ISUM
*
      RETURN
      END

