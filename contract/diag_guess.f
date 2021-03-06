*----------------------------------------------------------------------*
      subroutine diag_guess(mel_out,xlist,idxlist,nlist,iroot,isign,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     initialize list with diagonal guess for (iroot)th root
*     using the arrays xlist, idxlist which provide the lowest entries
*     (+ indices) of some diagonal approximation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'mdef_operator_info.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
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
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     jroot, idx, nset
      integer ::
     &     idxset(2)
      real(8) ::
     &     valset(2)

      integer ::
     &     idxlist_ba(nlist)

      ! if requested, get alpha<->beta interchanged indices
      if (isign.ne.0)
     &     call symidx_ab(idxlist_ba,
     &     idxlist,nlist,mel_out,
     &     op_info,str_info,strmap_info,orb_info)

      ! a) find out which elements to set
      jroot = 0
      idx = 1
      nset = 0
      do idx = 1, nlist
        if (isign.ne.0) then
          if (idxlist(idx).lt.abs(idxlist_ba(idx))) then
            ! set a spin combination
            jroot = jroot+1
            nset = 2
            idxset(1) = abs(idxlist(idx))
            idxset(2) = abs(idxlist_ba(idx))
            valset(1) = 1d0/sqrt(2d0)
            valset(2) = dble(isign)*dble(sign(1,idxlist_ba(idx)))
     &                  /sqrt(2d0)
          else if (idxlist(idx).eq.abs(idxlist_ba(idx))
     &           .and.isign.eq.+1) then
            if (idxlist_ba(idx).lt.0)
     &           call quit(1,'diag_guess','unexpected case')
            ! set a single element
            jroot = jroot+1
            nset = 1
            idxset(1) = abs(idxlist(idx))
            valset(1) = 1d0
          end if
        else
          ! set a single element
          jroot = jroot+1
          nset = 1
          idxset(1) = abs(idxlist(idx))
          valset(1) = 1d0
        end if
        if (jroot.eq.iroot) exit
      end do

      if (nset.eq.0)
     &     call quit(1,'diag_guess','nset==0 ??')

      if (iprlvl.ge.5) then
        if (nset.eq.1) then
          write(lulog,'(3x,"guess for root ",i3,": ",f12.6,i10)')
     &         iroot,valset(1),idxset(1)
        else
          write(lulog,'(3x,"guess for root ",i3,": ",2(f12.6,i10))')
     &         iroot,valset(1),idxset(1),valset(2),idxset(2)
        end if
      end if

      call set_list(mel_out,idxset,valset,nset)

      return
      end
