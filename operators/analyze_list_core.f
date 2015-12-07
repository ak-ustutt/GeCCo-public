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
      integer, parameter ::
     &     maxlist = 50

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
     &     maxblk, idx, ioff, iii
      real(8) ::
     &     xnrm, weight_last
      
      real(8), pointer ::
     &     xnorm(:,:), coeff_list(:), weight_list(:)
      integer, pointer ::
     &     idx_list(:)
      character(len=80), pointer ::
     &     label_list(:)


      maxblk=1
      do idx = 1, nvec
        maxblk = max(maxblk,me_cov(idx)%mel%op%n_occ_cls)
      end do

      ! get all info
      if (mode(1:4).eq.'NORM'.or.mode(1:4).eq.'norm') then

        allocate(xnorm(maxblk,nvec))
      
        do idx = 1, nvec
          call blk_norm_for_list(xnorm(1,idx),maxblk,
     &         me_cov(idx)%mel,me_contrv(idx)%mel)
        end do

        ! print it:
        write(luout,'(1x,60("-"))')
        write(luout,'(2x,a,i4)') 'record # ',idxrec
        write(luout,'(1x,60("-"))')
        do idx = 1, nvec
          call print_mel_contribs(luout,me_cov(idx)%mel,'e10.4',
     &         xnorm(1,idx),maxblk)
          write(luout,'(1x,60("-"))')
        end do

        deallocate(xnorm)

      else if (mode(1:7).eq.'CONTRIB'.or.mode(1:7).eq.'contrib') then

        allocate(label_list(maxlist),coeff_list(maxlist),
     &       weight_list(maxlist),idx_list(maxlist))

        ioff = 0
        do idx = 1, nvec
          call maxcontrb_for_list(
     &         coeff_list,weight_list,idx_list,
     &         maxlist,ioff,idx.eq.1,
     &         me_cov(idx)%mel,me_contrv(idx)%mel)
          ioff = ioff + me_cov(idx)%mel%len_op
        end do

        ! reduce list to 99% of norm (but cut behind degen. dets.)
        xnrm = 0d0
        iii = 1
        weight_last = 0d0
        do while(iii.lt.maxlist.and.
     &       (xnrm.lt.0.99d0.or.
     &       abs(weight_list(iii)-weight_last).lt.1d-10))
c          write(luout,'(i7,2x,f12.8,2x,f12.8,2x,f12.8)')
c     &         idx_list(iii),coeff_list(iii),weight_list(iii),
c     &         xnrm+weight_list(iii)
          weight_last = weight_list(iii)
          xnrm=xnrm+weight_list(iii)
          iii = iii+1
        end do
c        lenlist = iii-1

        ! idx_list -> labels
        ! (label_list,me_cov(idx)%mel,str_info,orb_info)
        
        ! must assemble indices for the contributing list and inquire labels
c        do iii = 1, lenlist
c          ioff = 0
c          do idx = 1, nvec
c            if (idx_list(iii).le.ioff+me_cov(idx)%mel%len_op) then
c              idx_list2(1,iii) = idx
c              idx_list2(2,iii) = idx_list(iii)-ioff
c            end if
c          end do
c        end do
c        
c        do idx = 1, nvec
c          jjj = 0
c          do iii = 1, nvec
c            if (idx_list2(1,iii).eq.idx) then
c              jjj = jjj+1
c              list3(jjj) = idx_list2(2,iii)
c            end if
c          end do
c        end do

        deallocate(label_list,coeff_list,
     &             weight_list,idx_list)
        
      else
        call quit(1,i_am,'not yet programmed!')
      end if


      end
