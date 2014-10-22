      subroutine analyze_list_core(me_cov,me_contrv,nvec,idxrec,mode,
     &     orb_info,str_info)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'mdef_operator_info.h'

       integer, parameter ::
     &     ntest = 0
      character(len=17), parameter ::
     &     i_am='analyze_list_core'

      integer, intent(in) ::
     &     nvec, idxrec
      type(me_list_array), intent(inout) ::
     &     me_cov(nvec), me_contrv(nvec) 
      character(len=*), intent(in) ::
     &     mode

      type(orbinf), intent(inout) ::
     &     orb_info
      type(strinf), intent(inout) ::
     &     str_info     

      integer ::
     &     maxblk, idx
      
      real(8), pointer ::
     &     xnorm(:,:)


      maxblk=1
      do idx = 1, nvec
        maxblk = max(maxblk,me_cov(idx)%mel%op%n_occ_cls)
      end do
      allocate(xnorm(maxblk,nvec))

      ! get all info
      if (mode(1:4).eq.'NORM'.or.mode(1:4).eq.'norm') then
        do idx = 1, nvec
          call blk_norm_for_list(xnorm(1,idx),maxblk,
     &         me_cov(idx)%mel,me_contrv(idx)%mel)
        end do
      else
        call quit(1,i_am,'not yet programmed!')
      end if

      ! print it:
      write(luout,'(1x,60("-"))')
      write(luout,'(2x,a,i4)') 'record # ',idxrec
      write(luout,'(1x,60("-"))')
      do idx = 1, nvec
        call print_mel_contribs(luout,me_cov(idx)%mel,'e10.4',
     &       xnorm(1,idx),maxblk)
        write(luout,'(1x,60("-"))')
      end do

      deallocate(xnorm)

      end
