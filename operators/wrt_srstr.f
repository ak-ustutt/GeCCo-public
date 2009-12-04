*----------------------------------------------------------------------*
      subroutine wrt_srstr(luout,irestr,ngasc,nspin)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     luout, ngasc, nspin, irestr(2,ngasc,2,nspin)
      integer, external ::
     &     ielsqsum

      if (nspin.eq.1) then
        write(luout,'(x,10(x,2i3))') irestr(1:2,1:ngasc,1,1)
        if (ielsqsum(irestr(1,1,2,1),2*ngasc).eq.0) then
          write(luout,'(x,a)') 'unmasked'
        else
          write(luout,'(x,a)') 'masked by:'
          write(luout,'(x,10(x,2i3))') irestr(1:2,1:ngasc,2,1)
        end if
      else
        write(luout,'(x,"a",10(x,2i3))') irestr(1:2,1:ngasc,1,1)
        write(luout,'(x,"b",10(x,2i3))') irestr(1:2,1:ngasc,1,2)
        if (ielsqsum(irestr(1,1,2,1),2*ngasc).eq.0 .and.
     &      ielsqsum(irestr(1,1,2,2),2*ngasc).eq.0) then
          write(luout,'(x,a)') 'unmasked'
        else
          write(luout,'(x,a)') 'masked by:'
          write(luout,'(x,"a",10(x,2i3))') irestr(1:2,1:ngasc,2,1)
          write(luout,'(x,"b",10(x,2i3))') irestr(1:2,1:ngasc,2,2)
        end if
      end if

      return
      end
