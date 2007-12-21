*----------------------------------------------------------------------*
      subroutine diag_guess(mel_out,xlist,idxlist,nlist,iroot,isign)
*----------------------------------------------------------------------*
*     initialize list with diagonal guess for (iroot)th root
*     using the arrays xlist, idxlist which provide the lowest entries
*     (+ indices) of some diagonal approximation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'

      real(8), parameter ::
     &     thresh = 1d-10

      type(me_list), intent(inout) ::
     &     mel_out
      integer, intent(in) ::
     &     nlist, isign, iroot
      integer, intent(in) ::
     &     idxlist(nlist)
      real(8), intent(in) ::
     &     xlist(nlist)

      integer ::
     &     jroot, idx, nset
      integer ::
     &     idxset(2)
      real(8) ::
     &     valset(2)

      ! a) find out which elements to set
      jroot = 0
      idx = 1
      nset = 0
      do while (idx.lt.nlist)
        if (isign.ne.0.and.abs(xlist(idx)-xlist(idx+1)).lt.thresh) then
          ! set a spin combination
          jroot = jroot+1
          nset = 2
          idxset(1) = idxlist(idx)
          idxset(2) = idxlist(idx+1)
          valset(1) = 1d0/sqrt(2d0)
          valset(2) = dble(isign)/sqrt(2d0)
          idx = idx+2
        else if (isign.eq.0.or.isign.eq.+1) then
          ! set a single element
          jroot = jroot+1
          nset = 1
          idxset(1) = idxlist(idx)
          valset(1) = 1d0
          idx = idx+1
        else
          idx = idx+1
        end if
        if (jroot.eq.iroot) exit
      end do

      if (nset.eq.0)
     &     call quit(1,'diag_guess','nset==0 ??')

      if (iprlvl.ge.5) then
        if (nset.eq.1) then
          write(luout,'(3x,"guess for root ",i3,": ",f12.6,i10)')
     &         iroot,valset(1),idxset(1)
        else
          write(luout,'(3x,"guess for root ",i3,": ",2(f12.6,i10))')
     &         iroot,valset(1),idxset(1),valset(2),idxset(2)
        end if
      end if

      call set_list(mel_out,idxset,valset,nset)

      return
      end
