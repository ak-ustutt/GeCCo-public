
      subroutine davidson_assemble_results(dvdsbsp,
     &     rvals,
     $     me_res, nlists, nroots, xnrm,
     $     xbuf1, xbuf2, nincore, lbuf)
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'



      integer, parameter::
     &     ntest =  00
      character(len=*),parameter::
     &     i_am="davidson_assemble_results"

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp


      integer,intent(in)::
     &     nroots,
     &     nlists,
     &     nincore, lbuf

      type(me_list_array),intent(in):: !lists change, array doesn't
     &     me_res(*)

      real(8),intent(out)::
     &     xnrm(nroots,nlists)

      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*),
     &     rvals(nroots,2)

      real(8),allocatable::
     &     eigenvecs(:,:)

      real(8) ::
     &     eigi(nroots),
     &     lxnrm(nlists),
     &     eigr(nroots)

      integer::
     &     leigenvec,
     &     iroot,
     &     ilist,
     &     maxlist
      integer,external::
     &     dvdsbsp_get_curlen
      real(8), external ::
     &     me_ddot




      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      leigenvec=dvdsbsp_get_curlen(dvdsbsp)
      allocate(eigenvecs(leigenvec,nroots))

      call dvdsbsp_get_eigenvec(dvdsbsp, eigenvecs, eigr, eigi,
     &     nroots, leigenvec)

      rvals(1:nroots,1) = eigr
      rvals(1:nroots,2) = eigi
      !print *,rvals
      maxlist=min(nroots,mel_get_maxrec(me_res(1)%mel)) !TODO remove the hardcoded one
      do iroot=1,maxlist
         do ilist=1,nlists
            call switch_mel_record(me_res(ilist)%mel,iroot)
         end do
         call vecsp_assemble_vector(dvdsbsp%vspace,
     &        me_res, nlists,
     &        eigenvecs(1:leigenvec,iroot) , leigenvec,
     &        0d0,
     &        xbuf1,xbuf2,nincore,lbuf)
         do ilist=1,nlists
            xnrm(iroot,ilist)=
     &           me_ddot(me_res(ilist)%mel,me_res(ilist)%mel,
     &           xbuf1, xbuf2, nincore, lbuf)
         end do
      end do
      deallocate(eigenvecs)
      return
      contains
      pure function mel_get_maxrec(mel)
      integer:: mel_get_maxrec
      type(me_list),intent(in)::mel
      mel_get_maxrec=mel%fhand%active_records(2)

      end function
      end subroutine
