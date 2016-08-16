


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
     &     ntest = 00
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
      real(8),dimension(:,:),pointer::
     &     vMv_mat
      integer::
     &     ncursub,icursub,nmaxsub

      type(filinf)::
     &     vfile, mvfile        !handles for the files for the input me-lists

      integer ::
     &     nold,                !number of vectors the product is to be calculated to
     &     ivvec,imvvec,         !index of the vector and M*v vector
     &     ii,jj,ilist,          !running inices
     &     lenlist
      real(8), external ::
     &     ddot, da_ddot

      ! alialising for shorter names
      ncursub=dvdsbsp%ncursub
      icursub=dvdsbsp%icursub
      nmaxsub=dvdsbsp%nmaxsub
      vMv_mat=> dvdsbsp%vMv_mat(:nmaxsub,:nmaxsub)
     

      if (ntest.ge.10)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on input:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)
      end if





      nold=ncursub
      if (ncursub .eq. nmaxsub) nold=nold-1  !one existing row/column is over written
      
      icursub=mod(icursub+1,nmaxsub) !icursub now indexes the to be overwritten vector


      if (nincore.ge.2)then
         do ilist=1,nlist
            lenlist=mvvec(ilist)%mel%len_op

            !copy the new lists into the vector spaces
            call vecsp_set_list_mel(dvdsbsp%mvspace,mvvec(ilist)%mel,
     &           icursub, ilist, buf1,lbuf) 

!unneccessary as currently buf1 now holds the mvvec(ilist)%mel elements
!            call vecsp_get_list_buf(dvdsbsp%mvspace, icursub, ilist, buf1, 
!     &              lenlist)
            
            call vecsp_set_list_mel(dvdsbsp%mvspace,mvvec(ilist)%mel,
     &           icursub, ilist, buf2,lbuf) 

            ! update the vMv matrix
            do jj=1,ncursub
               call vecsp_get_list_buf(dvdsbsp%vspace, jj, ilist, 
     &              lenlist, buf2, lbuf)
               if (ilist.eq.1) vMv_mat(jj,icursub)=0
               vMv_mat(jj,icursub) = vMv_mat(jj,icursub)
     &              + ddot(lenlist,buf1,1,buf2,1)
            end do

            call vecsp_get_list_buf(dvdsbsp%vspace, icursub,ilist,
     &      lenlist, buf2, lbuf) 

            do jj=1,ncursub              
               call vecsp_get_list_buf(dvdsbsp%vspace, jj, ilist, 
     &              lenlist, buf1, lbuf)

               if (jj.eq.icursub) cycle !don't calculated v_icursub * Mv_icursub twice

               if (ilist.eq.1) vMv_mat(icursub,jj)=0
               vMv_mat(icursub,jj) = vMv_mat(icursub,jj)
     &              + ddot(lenlist,buf1,1,buf2,1)
            end do
         end do
      else
         call quit(1,i_am, "not prepared for nincore <2.")
      end if
      return



      end subroutine
