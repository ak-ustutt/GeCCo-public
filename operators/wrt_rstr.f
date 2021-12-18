
*----------------------------------------------------------------------*
      subroutine wrt_rstr(lulog,irestr,ngas)
*----------------------------------------------------------------------*

      implicit none

      integer, intent(in) ::
     &     lulog, ngas, irestr(2,ngas,2,2)
      integer, external ::
     &     ielsqsum

      write(lulog,'(x,"C",10(x,2i3))') irestr(1:2,1:ngas,1,1)
      write(lulog,'(x,"A",10(x,2i3))') irestr(1:2,1:ngas,2,1)
      if (ielsqsum(irestr(1,1,1,2),2*2*ngas).eq.0) then
        write(lulog,'(x,a)') 'unmasked'
      else
        write(lulog,'(x,a)') 'masked by:'
        write(lulog,'(x,"C",10(x,2i3))') irestr(1:2,1:ngas,1,2)
        write(lulog,'(x,"A",10(x,2i3))') irestr(1:2,1:ngas,2,2)
      end if

      return
      end
