


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
      subroutine dvdsbsp_update(dvdsbsp, svvec, mvvec, nlist, buf1,buf2, 
     &     buf3, nincore,  lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      integer, parameter::
     &     ntest = 00,          !20 some scalar values; 60 input lists; 100 all used vectors (large)
     &     maxelem = 100        ! maximal maxelem elements per list are printed
      character(len=*),parameter::
     &     i_am="dvdsbsp_update"

      integer, intent(in)::
     &     nlist                !number of lists in new vectors

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      type(me_list_array),dimension(nlist)::
     &     svvec, mvvec          ! arrays with all me_lists signifiing the new vector and Mv vector

      integer,intent(in)::
     &     nincore,                 !how many buffer may be used?
     &     lbuf                  ! length of buffer
      real(8),dimension(*),intent(inout)::
     &     buf1, buf2, buf3

      
!     alialising the davidson subspace fields 
      real(8),dimension(:),pointer::
     &     vMv_mat,vSv_mat
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
      ncursub=dvdsbsp%ncursub
      nmaxsub=dvdsbsp%nmaxsub
      vMv_mat=>dvdsbsp%vMv_mat
      vSv_mat=>dvdsbsp%vSv_mat
     

      if (ntest.ge.10)
     &     call write_title(lulog,wst_dbg_subr,i_am)
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on input:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)
         if (dvdsbsp%with_metric)then
            write (lulog,*) 'overlapp matrix on input:'
            call wrtmat2(vSv_mat,ncursub,ncursub,
     &           nmaxsub,nmaxsub)
         end if
      end if



      dvdsbsp%lcursub=
     &     mod(dvdsbsp%lcursub,dvdsbsp%nmaxsub)+1

      if (nincore.ge.2)then
         do ilist=1,nlist
            lenlist=mvvec(ilist)%mel%len_op

            if (lenlist .eq. 0) then ! optimization, following code works also with empty lists (10.11.2016 Arne)
               if (ntest.ge.20)
     &              write (lulog,*) "Skipping empty list no.",ilist
               continue
            end if
            
            call dvdsbsp_set_lists(dvdsbsp,
     &           Mvvec(ilist)%mel, Svvec(ilist)%mel, ilist,
     &           buf1, lbuf )


            call dvdsbsp_update_matrices(dvdsbsp, ilist, lenlist,
     &           buf1, buf2, lbuf)
!calculate subspace and metric matrix
            
         end do 
      else !nincore
         call quit(1,i_am, "not prepared for nincore <2.")
      end if
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on output:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)
         
         if (dvdsbsp%with_metric)then
            write (lulog,*) 'overlapp matrix on input:'
            call wrtmat2(vSv_mat,ncursub,ncursub,
     &           nmaxsub,nmaxsub)
         end if
      end if
      return
      contains

!----------------------------------------------------------------------
!>   prints a buffer that contains elements of ME-list
!!
!!    @param lu output unit (open, formatted)
!!    @param mel me_list
!!    @param xbuf buffer with the elements
!!    @param maxelem number of elements to be printed
!!    @param msg additional message postfixed with the label of the me-list
!----------------------------------------------------------------------
      subroutine print_list_from_buffer(lu, mel,xbuf,maxelem, msg)
!----------------------------------------------------------------------
      implicit none

      integer, intent(in)::
     &     maxelem, lu
      type(me_list),intent(in)::
     &     mel
      real(8),dimension(*),intent(in)::
     &     xbuf
      character(len=*)::
     &     msg
      
      integer::
     &     idxdbg
      
      write(lu, *) msg,mel%label
      
      do idxdbg=1,maxelem
         write(lu, *) idxdbg,buf1(idxdbg)
      end do
      return
      end subroutine
!!----------------------------------------------------------------------
!>    add Mv and metric vector product vectors to a davidson subspace
!!
!!    @param dvdsbsp a davidson subspace handle
!!    @param me_mvp list of the matrix vector product
!!    @param me_svp pointer to the metric vectro product list (may be not associated )
!!    @param ilist index of the lists in the vectors
!!    @param xbuf scratch buffer
!!    @param lbuf length of scratch buffer
!-----------------------------------------------------------------------
      subroutine dvdsbsp_set_lists(dvdsbsp,me_mvp, me_svp, ilist, xbuf,
     &     lbuf )
!!----------------------------------------------------------------------

      implicit none
!! using includes of parent routine (lulog,
!     ! using ntest maxelem of parent routine
      character(len=*),parameter::
     &     i_am="dvdsbsp_update:set_list"
      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      type(me_list), intent(inout):: 
     &     me_mvp
      type(me_list),pointer, intent(in):: 
     &     me_svp
      integer, intent(in)::
     &     ilist,lbuf
      real(8),Dimension(*),intent(in)::
     &     xbuf
      
      integer::
     &     lcursub,
     &     lenlist
      integer::
     &     idxdbg
      
      lcursub=dvdsbsp%lcursub
      lenlist=me_mvp%len_op
      
      call vecsp_set_list_mel(dvdsbsp%mvspace,me_mvp,
     &     lcursub, ilist, xbuf, lbuf)
            
      if (ntest.ge.60)
     &     call print_list_from_buffer(lulog, me_mvp,
     &     xbuf, min(lenlist, maxelem), "Mv Product: ")

      if (dvdsbsp%with_metric)then
         if(associated(me_svp))then ! if there is a mel on svvec
            call vecsp_set_list_mel(dvdsbsp%Svspace,
     &           me_svp, lcursub, ilist, xbuf,lbuf)
            
            if (ntest.ge.60)
     &           call print_list_from_buffer(lulog, me_svp,
     &           xbuf, min(lenlist, maxelem), "Sv Product: ")
            
         else                   !otherwise use vvec
            call vecsp_get_list_buf(dvdsbsp%vspace, lcursub,ilist, 
     &           lenlist, xbuf, lbuf)
            
            if (ntest.ge.60)then
               write(lulog, *) "using v as Sv Product "
               do idxdbg=1,min(lenlist, maxelem)
                  write(lulog, *) idxdbg,xbuf(idxdbg)
               end do
            end if
            
            call vecsp_set_list_buf(dvdsbsp%Svspace,
     &           xbuf, lcursub, ilist, lbuf)
         end if   
      end if
      return
      end subroutine
!!----------------------------------------------------------------------

!!----------------------------------------------------------------------
      subroutine dvdsbsp_update_matrices(dvdsbsp, ilist, lenlist,
     &     buf1, buf2, lbuf)
!!----------------------------------------------------------------------
      implicit none
!! using includes of parent routine (lulog,
!     ! using ntest,maxelem of parent routine
      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer, intent(in)::
     &     ilist,
     &     lenlist,
     &     lbuf
      real(8),dimension(*),intent(inout)::
     &     buf1, buf2

      integer::
     &     ncursub,
     &     lcursub,
     &     nmaxsub

      integer::
     &     jj,
     &     idxdbg
      
      real(8),dimension(:),pointer::
     &     vMv_mat,vSv_mat
      
      nmaxsub=dvdsbsp%nmaxsub
      ncursub=dvdsbsp%ncursub
      lcursub=dvdsbsp%lcursub
      vMv_mat=>dvdsbsp%vMv_mat
      vSv_mat=>dvdsbsp%vSv_mat
      
      do jj=1,ncursub
         call vecsp_get_list_buf(dvdsbsp%vspace, jj, ilist,
     &        lenlist, buf2, lbuf)
         
         if (ntest.ge.100)then 
            write(lulog, *) "v vector", jj
            do idxdbg=1,min(lenlist, maxelem)
               write(lulog, *) idxdbg,buf2(idxdbg)
            end do 
         end if
         
         call vecsp_get_list_buf(dvdsbsp%mvspace, lcursub, ilist, 
     &        lenlist, buf1, lbuf)
         
         if (ilist.eq.1) vMv_mat(nmaxsub*(lcursub-1)+jj)=0 !if there was anything on that vector erase it
         
         vMv_mat(nmaxsub*(lcursub-1)+jj) =
     &        vMv_mat(nmaxsub*(lcursub-1)+jj)
     &        + ddot(lenlist,buf1,1,buf2,1)
         
         if (dvdsbsp%with_metric)then 
            if (ilist.eq.1) vSv_mat(nmaxsub*(lcursub-1)+jj)=0 !if there was anything on that vector erase it
            call vecsp_get_list_buf(dvdsbsp%Svspace, lcursub,
     &           ilist, lenlist, buf1, lbuf)
            
            vSv_mat(nmaxsub*(lcursub-1)+jj) =
     &           vSv_mat(nmaxsub*(lcursub-1)+jj)
     &           + ddot(lenlist,buf1,1,buf2,1)
         end if
      end do                    ! jj
      
      call vecsp_get_list_buf(dvdsbsp%vspace, lcursub,ilist,
     &     lenlist, buf2, lbuf)
      
      do jj=1,dvdsbsp%ncursub
         if (jj.eq.lcursub) cycle !don't calculated v_icursub * Mv_icursub twice
         if (jj.gt.lcursub .and. jj.le. dvdsbsp%icursub) cycle
         call vecsp_get_list_buf(dvdsbsp%mvspace, jj, ilist, 
     &        lenlist, buf1, lbuf)
         if (ilist.eq.1)  vMv_mat(nmaxsub*(jj-1)+lcursub)=0
         vMv_mat(nmaxsub*(jj-1)+lcursub) =
     &        vMv_mat(nmaxsub*(jj-1)+lcursub)
     &        + ddot(lenlist,buf1,1,buf2,1)
      end do
      if (dvdsbsp%with_metric)then
         do jj=1,dvdsbsp%ncursub
            if (jj.eq.lcursub) cycle !don't calculated v_lcursub * Sv_lcursub twice
            if (jj.gt.lcursub .and. jj.le. dvdsbsp%icursub)
     &           cycle
            
            call vecsp_get_list_buf(dvdsbsp%Svspace, jj, ilist, 
     &           lenlist, buf1, lbuf)
            if (ilist.eq.1)  vSv_mat(nmaxsub*(jj-1)+lcursub)=0
            vSv_mat(nmaxsub*(jj-1)+lcursub) =
     &           vSv_mat(nmaxsub*(jj-1)+lcursub)
     &           + ddot(lenlist,buf1,1,buf2,1)
         end do
      end if
      end subroutine
      
      end subroutine
