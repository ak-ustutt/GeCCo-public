      subroutine mollab(a,lu,luerr)
C
C  adapted from hjaaj (DALTON)
C
C  Purpose:
C     Search for MOLECULE labels on file LU
C
      character*8 a, b(4), c
      data c/'********'/
    1 read (lu,end=3,err=6) b
      if (b(1).ne.c) go to 1
      if (b(4).ne.a) go to 1
      if (luerr.lt.0) luerr = 0
      return
c
    3 if (luerr.lt.0) then
         luerr = -1
         return
      else
         write(luerr,4)a,lu
         call quit(0,'mollab','molecule label not found on file')
      end if
    4 format(/' MOLECULE label ',A8,
     *        ' not found on unit',I4)
C
    6 if (luerr.lt.0) then
         luerr = -2
         return
      else
         write (luerr,7) lu,a
         call quit(0,'mollab', 'error reading file')
      end if
    7 format(/' error reading unit',I4,
     *       /T22,'when searching for label ',A8)
      
      end
