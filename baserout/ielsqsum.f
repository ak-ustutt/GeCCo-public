
      FUNCTION IELSQSUM(IVEC,NELMNT)
*
* Square-sum of elements of integer vector IVEC
*
      DIMENSION IVEC(*)
*
      ISUM = 0
      DO IELMNT = 1, NELMNT
        ISUM = ISUM + IVEC(IELMNT)*IVEC(IELMNT)
      END DO
*
      IELSQSUM = ISUM
*
      RETURN
      END
