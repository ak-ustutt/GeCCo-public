*----------------------------------------------------------------------*
!>    ???????
!!     
*----------------------------------------------------------------------*
      subroutine dvdsbsp_append_vvec(dvdsbsp, me_lists, nlists,
     &     xbuf1, xbuf2, nincore, lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 1000
      character(len=*),parameter::
     &     i_am="dvdsbsp_append_vvec"
      real(8),parameter::
     &     thresh=1e-8
      integer,intent(in)::
     &     nlists,
     &     lbuf,nincore
     
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      type(me_list_array),intent(in)::
     &     me_lists(*)

      real(8)::
     &     xbuf1(*), xbuf2(*)
      real(8)::
     &     xnrm
      integer::
     &     nmaxsub, ncursub,icursub,
     &     ilist,
     &     idxdbg
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
         write(lulog,*) 'ncursub ',dvdsbsp%ncursub
      end if
      
      ncursub=dvdsbsp%ncursub
      icursub=dvdsbsp%icursub
      nmaxsub=dvdsbsp%nmaxsub
      
      if (nincore .lt.2) call quit(0,i_am, "at least 2 buffer required")
      call vecsp_orthvec(dvdsbsp%vspace, me_lists, nlists,
     &     xbuf1, xbuf2, lbuf)
      xnrm=vec_get_norm(me_lists,nlists,xbuf1,lbuf)
      if(ntest.ge.100)
     &     write(lulog,*)"old norm was",xnrm
      if(xnrm.lt.thresh)then
         call warn(
     &        i_am,"linear dependend trialvector rejected")
         return
      end if
      call vec_multiply(me_lists, nlists, 1/xnrm, xbuf1, lbuf)
      
      if (ncursub .ne. nmaxsub)then
         ncursub=ncursub+1
      end if
      icursub=mod(icursub,nmaxsub)+1
      if (ntest .gt. 30) then
         write (lulog,*) "subspace dimensions: pointer, current, max",
     &        icursub, ncursub, nmaxsub
      end if
      do ilist=1,nlists
         call vecsp_set_list_mel(dvdsbsp%vspace,me_lists(ilist)%mel,
     &        icursub, ilist, xbuf2,lbuf)
        if (ntest.ge.100)then 
           write(lulog, *) "trialvector",me_lists(ilist)%mel%label
           do idxdbg=1,me_lists(ilist)%mel%len_op
              write(lulog, *) idxdbg,xbuf2(idxdbg)
           end do
        end if
      end do
      dvdsbsp%ncursub=ncursub
      dvdsbsp%icursub=icursub

      contains
*----------------------------------------------------------------------*
      function vec_get_norm(me_lists, nlists, xbuf, lbuf)
*----------------------------------------------------------------------*
      implicit none
      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vec_get_norm"
      
      integer,intent(in)::
     &     nlists,
     &     lbuf
      type(me_list_array),intent(in)::
     &     me_lists(*)
      real(8),intent(inout)::
     &     xbuf(*)
      real(8)::vec_get_norm
      
      integer::
     &     ilist,
     &     lenlist,
     &     irec,
     &     ii
      real(8)::
     &     xnrm2
      type(filinf),pointer::
     &     ffme
      xnrm2=0
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf,lenlist)
         do ii=1,lenlist
            xnrm2=xnrm2+xbuf1(ii)**2
         end do
      end do
      vec_get_norm=sqrt(xnrm2)
      return
      end function
*----------------------------------------------------------------------*
      subroutine vec_multiply(me_lists, nlists,fac,  xbuf, lbuf)
*----------------------------------------------------------------------*
      implicit none
      integer, parameter::
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="vec_normalize"

      integer,intent(in)::
     &     nlists,
     &     lbuf
      type(me_list_array),intent(in)::
     &     me_lists(*)
      real(8),intent(in)::
     &     fac
      real(8),intent(inout)::
     &     xbuf(*)

      integer::
     &     ilist,
     &     lenlist,
     &     irec,
     &     ii
   
      type(filinf),pointer::
     &     ffme
      
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf1,lenlist)
         xbuf(1:lenlist)=xbuf(1:lenlist)*fac
         call vec_to_da(ffme,irec,xbuf1,lenlist)
      end do



      end subroutine
      
      end subroutine
