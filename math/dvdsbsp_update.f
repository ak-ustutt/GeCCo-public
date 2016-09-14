


*----------------------------------------------------------------------*
!> adds a new_vector to the  davidson subspace dvdsbsp and updates the v*Mv matrix
!!
!!   @param dvdsbsp a davidson subspace
!!   @param vvec a list of me-lists representing the v-vector (currently active record)
!!   @param mvvec same for the Mv-vector (currently active record)
!!   @param nlist number of lists
!!   @param buf1,buf2,buf3 !scratch space buf3 is currently not used no guarantee whats on it on exit
!!   @param nincore how many buffers are usable
!!   @param lbuf length of those buffers
*----------------------------------------------------------------------*
      subroutine dvdsbsp_update(dvdsbsp, vvec, mvvec, nlist, buf1, buf2, 
     &     buf3, nincore,  lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 1000
      character(len=*),parameter::
     &     i_am="dvdsbsp_update"

      integer, intent(in)::
     &     nlist                !number of lists in new vectors

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      type(me_list_array),dimension(nlist)::
     &     vvec, mvvec          ! arrays with all me_lists signifiing the new vector and Mv vector

      integer,intent(in)::
     &     nincore,                 !how many buffer may be used?
     &     lbuf                  ! length of buffer
      real(8),dimension(*),intent(inout)::
     &     buf1, buf2, buf3

      
!     alialising the davidson subspace fields 
      real(8),dimension(:),pointer::
     &     vMv_mat
      integer::
     &     lcursub,ncursub,nmaxsub

      type(filinf)::
     &     vfile, mvfile        !handles for the files for the input me-lists

      integer ::
     &     nold,                !number of vectors the product is to be calculated to
     &     ivvec,imvvec,         !index of the vector and M*v vector
     &     ii,jj,ilist,          !running inices
     &     lenlist,
     &     idxdbg
      real(8), external ::
     &     ddot, da_ddot

! alialising for shorter names
      lcursub=dvdsbsp%lcursub
      ncursub=dvdsbsp%ncursub
      nmaxsub=dvdsbsp%nmaxsub
      vMv_mat=>dvdsbsp%vMv_mat
     

      if (ntest.ge.10)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on input:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)
      end if



      lcursub=mod(lcursub,dvdsbsp%nmaxsub)+1

      if (nincore.ge.2)then
         do ilist=1,nlist
            lenlist=mvvec(ilist)%mel%len_op

            !copy the new lists into the vector spaces
            call vecsp_set_list_mel(dvdsbsp%mvspace,mvvec(ilist)%mel,
     &           lcursub, ilist, buf1,lbuf)
            
            if (ntest.ge.100)then
               write(lulog, *) "Mv Product",Mvvec(ilist)%mel%label
               do idxdbg=1,Mvvec(ilist)%mel%len_op
                  write(lulog, *) idxdbg,buf1(idxdbg)
               end do
            end if
            ! update the vMv matrix

!unneccessary as currently buf1 now holds the mvvec(ilist)%mel elements
!            call vecsp_get_list_buf(dvdsbsp%mvspace, icursub, ilist, buf1, 
!     &              lenlist)
            do jj=1,ncursub
               call vecsp_get_list_buf(dvdsbsp%vspace, jj, ilist, 
     &              lenlist, buf2, lbuf)
               if (ilist.eq.1) vMv_mat(nmaxsub*(lcursub-1)+jj)=0 !if there was anything on that vector erase it 
               vMv_mat(nmaxsub*(lcursub-1)+jj) =
     &              vMv_mat(nmaxsub*(lcursub-1)+jj)
     &              + ddot(lenlist,buf1,1,buf2,1)
            end do

            call vecsp_get_list_buf(dvdsbsp%vspace, lcursub,ilist,
     &      lenlist, buf2, lbuf) 

            do jj=1,dvdsbsp%ncursub
               if (jj.eq.lcursub) cycle !don't calculated v_icursub * Mv_icursub twice
               if (jj.gt.lcursub .and. jj.le. dvdsbsp%icursub) cycle
               ! if in multiple vectors are appended, only one is in this call appended. The other (between lcursub  icursub) will be set in a later call. more importantly they are not yet initiated if the mvspace is not fully filled.
               call vecsp_get_list_buf(dvdsbsp%mvspace, jj, ilist, 
     &              lenlist, buf1, lbuf)


               if (ilist.eq.1)  vMv_mat(nmaxsub*(jj-1)+lcursub)=0
               vMv_mat(nmaxsub*(jj-1)+lcursub) =
     &              vMv_mat(nmaxsub*(jj-1)+lcursub)
     &              + ddot(lenlist,buf1,1,buf2,1)
            end do
         end do
      else
         call quit(1,i_am, "not prepared for nincore <2.")
      end if
      dvdsbsp%lcursub=lcursub
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on output:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)
      end if
      return



      end subroutine
