
      subroutine davidson_assemble_residuals(dvdsbsp,
     &     rvals, nnew,
     $     me_res, nlists, nroots, xnrm,
     $     xbuf1, xbuf2, nincore, lbuf)
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'



      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="davidson_update_residuals"
      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      
      integer,intent(in)::
     &     nroots,
     &     nlists,
     &     nincore, lbuf

      integer,intent(inout)::
     &     nnew
      
      type(me_list_array),intent(in):: !lists change, array doesn't
     &     me_res(*)

      real(8),intent(out)::
     &     xnrm(nroots,nlists)
      
      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*),
     &     rvals(nroots)

      real(8),pointer::
     &     eigenvecs(:,:)
      
      real(8) ::
     &     eigi(nroots),
     &     lxnrm2(nlists)

      integer::
     &     leigenvec,
     &     iroot,
     &     ilist, ii
      integer,external::
     &     dvdsbsp_get_curlen
      


      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      leigenvec=dvdsbsp_get_curlen(dvdsbsp)
      nnew=min(nroots,dvdsbsp_get_curlen(dvdsbsp))
      allocate(eigenvecs(leigenvec,nroots))
      
      call dvdsbsp_get_eigenvec(dvdsbsp, eigenvecs, rvals, eigi,
     &     nroots, leigenvec)
      if (ntest.ge.20) then
         write (lulog,*) 'eigenvec'
         call wrtmat2(eigenvecs,leigenvec,nroots,
     &        leigenvec,nroots)
      end if
      
      do iroot=1,nnew
         do ilist=1,nlists
            call switch_mel_record(me_res(ilist)%mel,iroot)
         end do
         call dvdsbsp_assemble_residual(dvdsbsp,
     &        eigenvecs(1:leigenvec,iroot),
     &        rvals(iroot), leigenvec, me_res, nlists, lxnrm2,
     &        xbuf1, xbuf2, nincore, lbuf)
         do ilist=1,nlists
            xnrm(iroot,ilist)=sqrt(lxnrm2(ilist))
         end do
      end do
      deallocate(eigenvecs)
      end subroutine
