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
     &     ntest = 00
      character(len=*),parameter::
     &     i_am="dvdsbsp_orthvec"

      integer,intent(in)::
     &     nlists,
     &     lbuf,nincore
     
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      type(me_list_array),intent(in)::
     &     me_lists(*)

      real(8)::
     &     xbuf1(*), xbuf2(*)
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
      call vec_normalize(me_lists, nlists, xbuf1, lbuf)
      
      if (ncursub .ne. nmaxsub)then
         ncursub=ncursub+1
      end if
      icursub=mod(icursub+1,nmaxsub)
      
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
      subroutine vec_normalize(me_lists, nlists, xbuf, lbuf)
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
      real(8),intent(inout)::
     &     xbuf(*)

      integer::
     &     ilist,
     &     lenlist,
     &     irec,
     &     ii
      real(8)::
     &     xnrm
   
      type(filinf),pointer::
     &     ffme
      
      
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf,lenlist)
         do ii=1,lenlist
            xnrm=xnrm+xbuf1(ii)**2
         end do
      end do
      xnrm=sqrt(xnrm)
      
      do ilist=1,nlists
         ffme=> me_lists(ilist)%mel%fhand
         irec=me_lists(ilist)%mel%fhand%current_record
         lenlist= me_lists(ilist)%mel%len_op
         call vec_from_da(ffme,irec,xbuf1,lenlist)
         xbuf(1:lenlist)=xbuf(1:lenlist)/xnrm
         call vec_to_da(ffme,irec,xbuf1,lenlist)
      end do



      end subroutine
      
      end subroutine
