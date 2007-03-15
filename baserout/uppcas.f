* adapted from Jeppe:
      subroutine uppcas(line)
*
* Convert letters in character string LINE to upper case
*
* very stupid and not vectorized !
*
      implicit none
      include 'stdunit.h'

      character*(*), intent(inout) ::
     &     line

      integer, parameter ::
     &     ntest = 00,
     &     nchar = 43
      character*1 ::
     &     lower(nchar), upper(nchar)
*
      DATA LOWER/'a','b','c','d','e',
     &           'f','g','h','i','j',
     &           'k','l','m','n','o',
     &           'p','q','r','s','t',
     &           'u','v','w','x','y',
     &           'z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9','_',' '/
      DATA UPPER/'A','B','C','D','E',
     &           'F','G','H','I','J',
     &           'K','L','M','N','O',
     &           'P','Q','R','S','T',
     &           'U','V','W','X','Y',
     &           'Z','+','-','<','>',
     &           '=','0','1','2','3',
     &           '4','5','6','7','8',
     &           '9','_',' '/

      integer ::
     &     icha, i

      if (ntest.ge.100)
     &     write(luout,*) 'uppcas: on entry "',trim(line),'"'

      do icha = 1, len_trim(line)
        do i = 1,nchar
          if(line(icha:icha).eq.lower(i))
     &       line(icha:icha) = upper(i)
        end do
      end do

      if (ntest.ge.100)
     &     write(luout,*) 'uppcas: on exit "',trim(line),'"'
*
      return
      end
