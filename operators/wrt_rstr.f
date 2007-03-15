
*----------------------------------------------------------------------*
      subroutine wrt_rstr(luout,irestr,ngas)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     luout, ngas, irestr(2,ngas,2,2)
      integer, external ::
     &     ielsqsum

      write(luout,'(x,"C",10(x,2i3))') irestr(1:2,1:ngas,1,1)
      write(luout,'(x,"A",10(x,2i3))') irestr(1:2,1:ngas,2,1)
      if (ielsqsum(irestr(1,1,1,2),2*2*ngas).eq.0) then
        write(luout,'(x,a)') 'unmasked'
      else
        write(luout,'(x,a)') 'masked by:'
        write(luout,'(x,"C",10(x,2i3))') irestr(1:2,1:ngas,1,2)
        write(luout,'(x,"A",10(x,2i3))') irestr(1:2,1:ngas,2,2)
      end if

      return
      end
