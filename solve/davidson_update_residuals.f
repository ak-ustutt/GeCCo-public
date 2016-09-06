
      subroutine davidson_update_residuals(dvdsbsp,
     &     rvals, 
     $     me_res, nlists, nroots, xnrm,
     &     update,
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
      
      type(me_list_array),intent(in):: !lists change, array doesn't
     &     me_res(*)

      real(8),intent(out)::
     &     xnrm(nroots,nlists)
      
      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*),
     &     rvals(nroots)
      logical,intent(in)::
     &     update(nroots)

      real(8),allocatable::
     &     eigenvecs(:,:)
      
      real(8) ::
     &     eigi(nroots),
     &     lxnrm(nlists)

      integer::
     &     leigenvec,
     &     iroot,
     &     ilist
      integer,external::
     &     dvdsbsp_get_curlen
      


      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      leigenvec=dvdsbsp_get_curlen(dvdsbsp)
      allocate(eigenvecs(leigenvec,nroots))
      
      call dvdsbsp_get_eigenvec(dvdsbsp, eigenvecs, rvals, eigi,
     &     nroots, leigenvec)

      do iroot=1,nroots
         if (.not. update(iroot)) cycle
         call dvdsbsp_assemble_residual(dvdsbsp,
     &        eigenvecs(1:leigenvec,iroot),
     &        rvals(iroot), leigenvec, me_res, nlists, lxnrm,
     &        xbuf1, xbuf2, nincore, lbuf)
         do ilist=1,nlists
            xnrm(iroot,ilist)=lxnrm(ilist)
         end do
      end do
      deallocate(eigenvecs)
      end subroutine
