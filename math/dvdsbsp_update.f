


*----------------------------------------------------------------------*
!> appends covariant vectors (Mv) and possibly Sv to the dvdsbsp
!!
!!   also updates the reduced subspace and overlap matrices
!!   @param dvdsbsp a davidson subspace
!!   @param vvec a list of me-lists representing the v-vector (currently active record)
!!   @param mvvec same for the Mv-vector (currently active record)
!!   @param nlist number of lists
!!   @param buf1,buf2,buf3 !scratch space buf3 is currently not used !no guarantee whats on buffers on exit
!!   @param nincore how many buffers are usable
!!   @param lbuf length of those buffers
*----------------------------------------------------------------------*
      subroutine dvdsbsp_append_vecs(dvdsbsp, svvec, mvvec, vvec, nrec,
     &     nlists, buf1, buf2, buf3, nincore,  lbuf)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_file_array.h'
      include 'def_davidson_subspace.h'

      real(8),parameter::
     &     thresh = 1.0d-16
      integer, parameter::
     &     ntest = 000,          !20 some scalar values; 60 input lists; 100 all used vectors (large)
     &     maxelem = 100        ! maximal maxelem elements per list are printed
      character(len=*),parameter::
     &     i_am="dvdsbsp_update"

      integer, intent(in)::
     &     nlists                !number of lists in new vectors
      integer,intent(inout)::
     &     nrec
      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp

      type(me_list_array),dimension(nlists),intent(inout)::
     &     svvec, mvvec, vvec
! arrays with all me_lists signifiing the new vector and Mv vector

      integer,intent(in)::
     &     nincore,                 !how many buffer may be used?
     &     lbuf                  ! length of buffer
      real(8),dimension(*),intent(inout)::
     &     buf1, buf2, buf3 !buf3 currently not used


!     alialising the davidson subspace fields
      real(8),dimension(:),pointer::
     &     vMv_mat,vSv_mat
      integer::
     &     lcursub,ncursub,nmaxsub

      type(filinf)::
     &     vfile, mvfile        !handles for the files for the input me-lists

      integer ::
     &     nnew,                !number of vectors the product is to be calculated to
     &     irec,
     &     ivvec,imvvec,         !index of the vector and M*v vector
     &     ii,jj,ilist,          !running inices
     &     lenlist,
     &     idxdbg
      real(8), external ::
     &     ddot, da_ddot
      logical::
     &     accepted

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
 



      nnew=0
      if (nincore.ge.2)then
         do irec=1,nrec
            do ilist=1,nlists
               call switch_mel_record(vvec(ilist)%mel,irec)
               call switch_mel_record(mvvec(ilist)%mel,irec)
               if(dvdsbsp%use_metric(ilist))then
                  call switch_mel_record(svvec(ilist)%mel,irec)
               end if
            end do
            call dvdsbsp_append_lists(dvdsbsp, 
     &           Mvvec,
     &           Svvec,
     &           vvec,
     &           nlists,
     &           accepted,
     &           buf1,buf2,lbuf)
            if (accepted)then
               nnew=nnew+1
            else
               call warn(i_am, "linear dependend trialvector rejected")
               cycle
            end if 
            do ilist=1,nlists
               lenlist=mvvec(ilist)%mel%len_op

               if (lenlist .eq. 0) then ! optimization, following code works also with empty lists (10.11.2016 Arne)
                  if (ntest.ge.20)
     &                 write (lulog,*) "Skipping empty list no.",ilist
                  cycle
               end if
               call dvdsbsp_update_matrices(dvdsbsp, ilist, lenlist,
     &              buf1, buf2, lbuf)
!calculate subspace and metric matrix
            end do
         end do
      else                      !nincore
         call quit(1,i_am, "not prepared for nincore <2.")
      end if
      ncursub=dvdsbsp%ncursub
      if (ntest.ge.20) then
         write (lulog,*) 'subspace matrix on output:'
         call wrtmat2(vMv_mat,ncursub,ncursub,
     &        nmaxsub,nmaxsub)

         if (dvdsbsp%with_metric)then
            write (lulog,*) 'overlapp matrix on output:'
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
!----------------------------------------------------------------------
      subroutine dvdsbsp_append_lists(dvdsbsp, 
     &     me_MVlists,
     &     me_svlists,
     &     me_vlists,
     &     nlists,
     &     accepted,
     &     buf1,buf2,lbuf)
!----------------------------------------------------------------------
      implicit none
      
      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer,intent(in)::
     &     lbuf,
     &     nlists
      logical, intent(out)::
     &     accepted
      real(8),intent(inout)::
     &     buf1(lbuf), buf2(lbuf)
      
      type(me_list_array),intent(inout)::
     &     me_vlists(nlists), me_svlists(nlists),me_MVlists(nlists)
      
      type(filinf),pointer::
     &     ffme,ffcme
      real(8)::
     &     norm2
      integer::
     &     ilist,ivec,idx
      
      call dvdsbsp_orthvecs(dvdsbsp, 
     &     me_MVlists,
     &     me_svlists,
     &     me_vlists,
     &     nlists,
     &     buf1,buf2,lbuf)

      
      norm2=0
      if (ntest.ge.20) write(lulog,*) 'norm of input list:'
      do ilist=1,nlists
         ffme=> me_vlists(ilist)%mel%fhand
         irec= me_vlists(ilist)%mel%fhand%current_record
         lenlist= me_vlists(ilist)%mel%len_op
         if (dvdsbsp%use_metric(ilist)) then
            ffcme => me_svlists(ilist)%mel%fhand
         else
            ffcme => me_vlists(ilist)%mel%fhand
         end if
         norm2 =norm2+ da_ddot(ffme,irec, ffcme,irec,
     &        lenlist,buf1,buf2,lbuf)
         if (ntest.ge.20) write(lulog,*) ilist,norm2
      end do
      if (norm2 .gt. thresh) then
         accepted =.true.
         ivec = mod(dvdsbsp%ncursub,
     &        dvdsbsp%nmaxsub)+1
         dvdsbsp%ncursub=ivec
         dvdsbsp%lcursub=
     &        mod(dvdsbsp%lcursub,dvdsbsp%nmaxsub)+1
         do ilist=1,nlists
            irec= me_vlists(ilist)%mel%fhand%current_record
            lenlist= me_vlists(ilist)%mel%len_op
            call vec_from_da(me_vlists(ilist)%mel%fhand,
     &           irec,buf1,lenlist)
            do idx=1,lenlist
               buf1(idx) =1/sqrt(norm2)*buf1(idx)
            end do
            call  vecsp_set_list_buf(dvdsbsp%vspace,  buf1, ivec, ilist,
     &           lbuf)
            if(dvdsbsp%with_metric)then
               if (dvdsbsp%use_metric(ilist))then
                  ffme=> me_Svlists(ilist)%mel%fhand
               else
                  ffme=> me_vlists(ilist)%mel%fhand
               end if
               call vec_from_da(ffme,
     &              irec,buf1,lenlist)
               do idx=1,lenlist
                  buf1(idx) =1/sqrt(norm2)*buf1(idx)
               end do
               call  vecsp_set_list_buf(dvdsbsp%Svspace,  buf1, ivec, 
     &              ilist,lbuf)
            end if
            call vec_from_da(me_mvlists(ilist)%mel%fhand,
     &           irec,buf1,lenlist)
            do idx=1,lenlist
               buf1(idx) =1/sqrt(norm2)*buf1(idx)
            end do
            call  vecsp_set_list_buf(dvdsbsp%mvspace,  buf1, ivec, 
     &           ilist,lbuf)
            end do
         else
            accepted=.false.
            return
         end if
      end subroutine
!----------------------------------------------------------------------
      subroutine dvdsbsp_orthvecs(dvdsbsp,
     &     me_MVlists,
     &     me_svlists,
     &     me_vlists,
     &     nlists,
     &     buf1,buf2,lbuf)
!----------------------------------------------------------------------

      type(davidson_subspace_t),intent(inout)::
     &     dvdsbsp
      integer,intent(in)::
     &     lbuf,
     &     nlists
      real(8),intent(inout)::
     &     buf1(lbuf), buf2(lbuf)
      real(8)::
     &     overlapp
      type(me_list_array),intent(inout)::
     &     me_vlists(nlists), me_Svlists(nlists),me_MVlists(nlists)
      integer::
     &     ivec,irec, lenlist,
     &     idx
      
      type(filinf),pointer::
     &     ffme,ffcme

      if (dvdsbsp%with_metric)then
         if (ntest.ge.20) write(lulog,*) 'orthogonalizing using metric'
         do ivec= 1,dvdsbsp%ncursub
            if (ntest.ge.20) write(lulog,*) ' ivec = ',ivec
            overlapp=0
            do ilist=1,nlists
               ffme=> me_vlists(ilist)%mel%fhand
               irec= me_vlists(ilist)%mel%fhand%current_record
               lenlist= me_vlists(ilist)%mel%len_op
               call vec_from_da(ffme,irec,buf1,lenlist)
               if (dvdsbsp%use_metric(ilist))then
                  call  vecsp_get_list_buf(dvdsbsp%Svspace, ivec, ilist,
     &                 lenlist,
     &                 buf2, lbuf)
               else
                  call  vecsp_get_list_buf(dvdsbsp%vspace, ivec, ilist,
     &                 lenlist,
     &                 buf2, lbuf)
               end if
               overlapp = overlapp + ddot(lenlist,buf1,1,buf2,1)
               if (ntest.ge.20) then
                 write(lulog,"(1x,i4,3(1x,f16.12))") 
     &                ilist,ddot(lenlist,buf1,1,buf1,1),
     &                      ddot(lenlist,buf1,1,buf2,1),overlapp
               end if
            end do
            
            do ilist=1,nlists
               ffme=> me_vlists(ilist)%mel%fhand
               irec= me_vlists(ilist)%mel%fhand%current_record
               lenlist= me_vlists(ilist)%mel%len_op
               call vec_from_da(ffme,irec,buf1,lenlist)
               call  vecsp_get_list_buf(dvdsbsp%vspace, ivec, ilist,
     &              lenlist,
     &              buf2, lbuf)
               do idx=1,lenlist
                  buf1(idx) = buf1(idx) - overlapp*buf2(idx)
               end do
               call vec_to_da(ffme, irec,buf1,lenlist)

               ffme=> me_mvlists(ilist)%mel%fhand
               call vec_from_da(ffme,irec,buf1,lenlist)
               call  vecsp_get_list_buf(dvdsbsp%Mvspace, ivec, ilist,
     &              lenlist,
     &              buf2, lbuf)
               do idx=1,lenlist
                  buf1(idx) = buf1(idx) - overlapp*buf2(idx)
               end do
               call vec_to_da(ffme, irec,buf1,lenlist)
               if (dvdsbsp%use_metric(ilist)) then
                  ffme=> me_Svlists(ilist)%mel%fhand
                  call vec_from_da(ffme,irec,buf1,lenlist)
                  call  vecsp_get_list_buf(dvdsbsp%Svspace, ivec, ilist,
     &                 lenlist,
     &                 buf2, lbuf)
                  do idx=1,lenlist
                     buf1(idx) = buf1(idx) - overlapp*buf2(idx)
                  end do
                  call vec_to_da(ffme, irec,buf1,lenlist)
               end if
            end do
         end do
      else
         continue !was already done in dvdsbsp_orth_vvec
      end if
      end subroutine

      
!!----------------------------------------------------------------------
      subroutine dvdsbsp_update_matrices(dvdsbsp, ilist, lenlist,
     &     buf1, buf2, lbuf)
!!----------------------------------------------------------------------
      implicit none
!! using includes of parent routine (lulog,
!     ! using ntest,maxelem of parent routine
      integer,parameter::
     &     ntest=00,
     &     maxelem=100
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
      lcursub=dvdsbsp%ncursub
      vMv_mat=>dvdsbsp%vMv_mat
      vSv_mat=>dvdsbsp%vSv_mat


      !add a row
      do jj=1,ncursub
         call vecsp_get_list_buf(dvdsbsp%vspace, jj, ilist,
     &        lenlist, buf2, lbuf)

         if (ntest.ge.60)then
            write(lulog, *) "v vector", jj
            do idxdbg=1,min(lenlist, maxelem)
               write(lulog, *) idxdbg,buf2(idxdbg)
            end do
         end if

         call vecsp_get_list_buf(dvdsbsp%mvspace, lcursub, ilist,
     &        lenlist, buf1, lbuf)

         if (ilist.eq.1) vMv_mat(nmaxsub*(lcursub-1)+jj)=0 !if there was anything on that vector erase it

         if(ntest.ge.20 )then
            write(lulog, *) "vMv, lcursub, jj",lcursub, jj
            write(lulog, *) "old:",vMv_mat(nmaxsub*(lcursub-1)+jj)
            write(lulog, *) "added:",ddot(lenlist,buf1,1,buf2,1)
         end if
         vMv_mat(nmaxsub*(lcursub-1)+jj) =
     &        vMv_mat(nmaxsub*(lcursub-1)+jj)
     &        + ddot(lenlist,buf1,1,buf2,1)

         if (dvdsbsp%with_metric)then
            if (ilist.eq.1) vSv_mat(nmaxsub*(lcursub-1)+jj)=0 !if there was anything on that vector erase it
            call vecsp_get_list_buf(dvdsbsp%Svspace, lcursub,
     &           ilist, lenlist, buf1, lbuf)

            if(ntest.ge.20 )then
               write(lulog, *) "vSv, lcursub, jj",lcursub, jj
               write(lulog, *) "old:",vSv_mat(nmaxsub*(lcursub-1)+jj)
               write(lulog, *) "added:",ddot(lenlist,buf1,1,buf2,1)
            end if
            vSv_mat(nmaxsub*(lcursub-1)+jj) =
     &           vSv_mat(nmaxsub*(lcursub-1)+jj)
     &           + ddot(lenlist,buf1,1,buf2,1)
         end if
      end do                    ! jj

      call vecsp_get_list_buf(dvdsbsp%vspace, lcursub,ilist,
     &     lenlist, buf2, lbuf)

      ! add a column
      do jj=1,dvdsbsp%ncursub
         if (jj.eq.lcursub) cycle !don't calculated v_icursub * Mv_icursub twice
         if (jj.gt.lcursub .and. jj.le. dvdsbsp%icursub) cycle
         call vecsp_get_list_buf(dvdsbsp%mvspace, jj, ilist,
     &        lenlist, buf1, lbuf)
         if (ilist.eq.1)  vMv_mat(nmaxsub*(jj-1)+lcursub)=0
         if(ntest.ge.20 )then
            write(lulog, *) "vMv, lcursub, jj",lcursub, jj
            write(lulog, *) "old:",vMv_mat(nmaxsub*(lcursub-1)+jj)
            write(lulog, *) "added:",ddot(lenlist,buf1,1,buf2,1)
         end if
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
            if(ntest.ge.20 )then
               write(lulog, *) "vSv, lcursub, jj",lcursub, jj
               write(lulog, *) "old:",vSv_mat(nmaxsub*(lcursub-1)+jj)
               write(lulog, *) "added:",ddot(lenlist,buf1,1,buf2,1)
            end if
            vSv_mat(nmaxsub*(jj-1)+lcursub) =
     &           vSv_mat(nmaxsub*(jj-1)+lcursub)
     &           + ddot(lenlist,buf1,1,buf2,1)
         end do
      end if
      end subroutine

      end subroutine
