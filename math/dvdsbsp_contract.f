


      subroutine dvdsbsp_compress(dvdsbsp, nretain, nlists, me_scr,
     &     xbuf1, xbuf2,lbuf, nincore)
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 1000
      character(len=*),parameter::
     &     i_am="dvdsbsp_compress"




      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      real(8),intent(inout)::
     &     xbuf1(*), xbuf2(*)

      integer,intent(in)::
     &     nlists,nretain,
     &     lbuf,nincore
      type(me_list_array),intent(inout)::
     &     me_scr(*)
      
      real(8), allocatable::
     &     eigi(:),eigr(:),
     &     vecs(:,:)
      integer::
     &     ndim,idxdbg

      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if

      if(dvdsbsp%ncursub.ne.dvdsbsp%lcursub)
     &     call quit(1,i_am,"Subspace must be be up to date") !same number of Mv vectors as of Vvectors
      if(dvdsbsp%ncursub.lt.nretain)
     &     call quit(1,i_am,"I can't expand the subspace") !same number of Mv vectors as of Vvectors

      
      ndim=dvdsbsp%ncursub
      
      allocate(eigi(nretain),eigr(nretain),vecs(ndim,nretain))

      
      
      call dvdsbsp_get_eigenvec(dvdsbsp, vecs, eigr, eigi, nretain,
     &     dvdsbsp%ncursub)
      
      call vecsp_compress(dvdsbsp%vspace,vecs,nretain,nlists,ndim,
     &     me_scr, xbuf1, xbuf2, lbuf)
      call vecsp_compress(dvdsbsp%Mvspace,vecs,nretain,nlists,ndim,
     &     me_scr, xbuf1, xbuf2, lbuf)
      if (dvdsbsp%with_metric)then
         call vecsp_compress(dvdsbsp%Svspace,vecs,nretain,nlists,ndim,
     &        me_scr, xbuf1, xbuf2, lbuf)
      end if
      
      maxsub=dvdsbsp%nmaxsub
      dvdsbsp%vMv_mat(1:maxsub*maxsub) = 0d0
      
      call vecsp_calculate_overlapp(dvdsbsp%vspace,dvdsbsp%Mvspace,
     &     xbuf1,xbuf2, lbuf, dvdsbsp%vMv_mat, maxsub)
      
      deallocate(eigi, eigr,vecs)
      dvdsbsp%ncursub=nretain
      dvdsbsp%icursub=nretain
      dvdsbsp%lcursub=nretain
      return
      contains
      
      subroutine vecsp_compress(vecsp, vecs, nretain, nlists, ndim,
     &     me_scr,xbuf1, xbuf2,lbuf)
      implicit none
      integer,parameter ::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="vecsp_compress"
      
      type(vector_space_t), intent(inout)::
     &     vecsp
      integer,intent(in)::
     &     nlists,nretain,
     &     lbuf, ndim
      type(me_list_array),intent(inout)::
     &     me_scr(nlists)
      real(8),intent(inout)::
     &     vecs(ndim,nretain),
     &     xbuf1(*),xbuf2(*)
      integer::
     &     ilist, lenlist, iretain, ivec,idx
      
      if (ntest.ge.100) then
         call write_title(lulog,wst_dbg_subr,i_am)
      end if
      if (ntest.ge.20) then
         write (lulog,*) 'retained vectors'
         call wrtmat2(vecs,ndim,nretain,
     &        ndim,nretain)
      end if
      
      do ilist=1,nlists
         lenlist=vecsp%me_lists(ilist)%mel%len_op
         do iretain=1,nretain
            xbuf2(1:lenlist)=0d0
            do ivec=1,ndim
               call vecsp_get_list_buf(vecsp,ivec,ilist,
     &              lenlist, xbuf1, lbuf)
               
               xbuf2(1:lenlist)=xbuf2(1:lenlist)
     &              +xbuf1(1:lenlist)*vecs(ivec,iretain)
               
               if (ntest.gt.100)then
                  write (lulog,*) "added vec:",ivec," factor",
     &                 vecs(ivec,iretain)
                  write(lulog,*) "accumulated   | added (wo_factor)"
                  do idx=1,lenlist
                     write (lulog,*)  xbuf2(idx),xbuf1(idx)
                  end do
               end if
            end do
            call switch_mel_record(me_scr(ilist)%mel,iretain) 
            call vec_to_da(me_scr(ilist)%mel%fhand,iretain,
     &           xbuf2,lenlist) !intermediate save
         end do 
      end do
      do ilist=1,nlists
         lenlist=vecsp%me_lists(ilist)%mel%len_op
         do iretain=1,nretain
            call switch_mel_record(me_scr(ilist)%mel,iretain) 
            call vecsp_set_list_mel(vecsp,me_scr(ilist)%mel, iretain,
     &            ilist, xbuf1, lbuf)
         end do
      end do
      vecsp%nvec=nretain
      end subroutine

      subroutine vecsp_calculate_overlapp(  vecsp1, vecsp2,
     &     xbuf1,xbuf2, lbuf ,
     &     mat, matdim)
      implicit none
      integer,parameter ::
     &     ntest=00
      character(len=*),parameter::
     &     i_am="vecsp_calculate_overlapp"


      type(vector_space_t), intent(inout)::
     &     vecsp1, vecsp2
      
      integer, intent(in)::
     &     lbuf, matdim

      real(8),intent(inout)::
     &     xbuf1(*),xbuf2(*),
     &     mat(matdim*matdim)

      integer::
     &     nvec1, nvec2,
     &     ivec1, ivec2,
     &     nlist, ilist,
     &     lenlist
      real(8),external::
     &     ddot
      
      if (ntest.ge.100)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      nvec1=vecsp1%nvec
      nvec2=vecsp2%nvec
      nlist=vecsp1%nlists
      if (vecsp2%nlists.ne.vecsp1%nlists)then
         call quit(1,i_am, "vecspaces have to be on "//
     &        "vectors of the same number of lists")
      end if
      do ilist=1,nlist
         do ivec2=1,nvec2
            
            lenlist=vecsp2%me_lists(ilist)%mel%len_op
         
            call vecsp_get_list_buf(vecsp2,ivec2,ilist,
     &           lenlist, xbuf2, lbuf)
            do ivec1=1,nvec1
               call vecsp_get_list_buf(vecsp1,ivec1,ilist,
     &              lenlist, xbuf1, lbuf)
               mat(matdim*(ivec2-1)+ivec1)=
     &              mat(matdim*(ivec2-1)+ivec1)
     &              + ddot(lenlist,xbuf1,1,xbuf2,1)
            end do
         end do
      end do
      end subroutine
      end subroutine 
