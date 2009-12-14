*----------------------------------------------------------------------*
      subroutine special_restr(irestr_op1op2,
     &     iocc_op1op2,njoined_op1op2,
     &     merge_op1op2,
     &     iocc_op1,iocc_ex1,irestr_op1,njoined_op1,
     &     iocc_op2,iocc_ex2,irestr_op2,njoined_op2,
     &     hpvxgas,nspin,ngas)
*----------------------------------------------------------------------*
*     more advanced (but simple) version to generate the restrictions
*     on intermediates. so far, we are allowed to just add up 
*     restrictions (works for simple frozen core stuff, w/o masking)
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nspin, ngas, njoined_op1, njoined_op2, njoined_op1op2
      integer, intent(inout) ::
     &     irestr_op1op2(2,ngas,2,2,nspin,njoined_op1op2)
      integer, intent(in) ::
     &     merge_op1op2(*), iocc_op1op2(ngastp,2,njoined_op1op2),
     &     iocc_op1(ngastp,2,njoined_op1),
     &     iocc_ex1(ngastp,2,njoined_op1),
     &     irestr_op1(2,ngas,2,2,nspin,njoined_op1),
     &     iocc_op2(ngastp,2,njoined_op2),
     &     iocc_ex2(ngastp,2,njoined_op2),
     &     irestr_op2(2,ngas,2,2,nspin,njoined_op2),
     &     hpvxgas(ngas,nspin)

      integer ::
     &     ivtx12, ivtx1, ivtx2, idx, ica, ispin, hpvx, igas,
     &     n1, n2, ii
      integer ::
     &     irestr_ex1(2,ngas,2,2,nspin,njoined_op1),
     &     irestr_ex2(2,ngas,2,2,nspin,njoined_op2)

      logical, external ::
     &     list_cmp
     
      call fit_restr2(irestr_ex1,
     &     iocc_ex1,njoined_op1,irestr_op1,njoined_op1,hpvxgas,ngas)
      call fit_restr2(irestr_ex2,
     &     iocc_ex2,njoined_op2,irestr_op2,njoined_op2,hpvxgas,ngas)

      irestr_op1op2 = 0
 
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'special_restr')
        write(luout,*) 'on input: 12'
        do ivtx12 = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op1op2(:,:,ivtx12),
     &         irestr_op1op2(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
        write(luout,*) ' 1'
        do ivtx12 = 1, njoined_op1
          call wrt_occ_rstr(luout,ivtx12,iocc_op1(:,:,ivtx12),
     &         irestr_op1(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
        write(luout,*) ' 1ex'
        do ivtx12 = 1, njoined_op1
          call wrt_occ_rstr(luout,ivtx12,iocc_ex1(:,:,ivtx12),
     &         irestr_ex1(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
        write(luout,*) ' 2'
        do ivtx12 = 1, njoined_op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op2(:,:,ivtx12),
     &         irestr_op2(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
        write(luout,*) ' 2ex'
        do ivtx12 = 1, njoined_op2
          call wrt_occ_rstr(luout,ivtx12,iocc_ex2(:,:,ivtx12),
     &         irestr_ex2(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
        write(luout,*) 'map:'
        call print_mapinfo(luout,merge_op1op2,njoined_op1op2)
      end if


      idx = 1
      do ivtx12 = 1, njoined_op1op2
        n1 = merge_op1op2(idx)
        n2 = merge_op1op2(idx+n1+1)

        do ii = 1, n1
          ivtx1 = merge_op1op2(idx+ii)
          do ica = 1, 2
            do hpvx = 1, ngastp
              if (iocc_ex1(hpvx,ica,ivtx1).eq.0) cycle
              do ispin = 1, nspin
                do igas = 1, ngas
                  if (hpvxgas(igas,ispin).eq.hpvx) then
                      irestr_op1op2(1:2,igas,ica,1:2,ispin,ivtx12) =
     &                irestr_op1op2(1:2,igas,ica,1:2,ispin,ivtx12) +
     &                   irestr_ex1(1:2,igas,ica,1:2,ispin,ivtx1)
                  end if
                end do
              end do
            end do
          end do
        end do

        do ii = 1, n2
          ivtx2 = merge_op1op2(idx+n1+1+ii)
          do ica = 1, 2
            do hpvx = 1, ngastp
              if (iocc_ex2(hpvx,ica,ivtx2).eq.0) cycle
              do ispin = 1, nspin
                do igas = 1, ngas
                  if (hpvxgas(igas,ispin).eq.hpvx) then
                      irestr_op1op2(1:2,igas,ica,1:2,ispin,ivtx12) =
     &                irestr_op1op2(1:2,igas,ica,1:2,ispin,ivtx12) +
     &                   irestr_ex2(1:2,igas,ica,1:2,ispin,ivtx2)
                  end if
                end do
              end do
            end do
          end do
        end do

        idx = idx + n1 + n2 + 2

      end do

      if (ntest.ge.100) then
        write(luout,*) 'on output:'
        do ivtx12 = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op1op2(:,:,ivtx12),
     &         irestr_op1op2(:,:,:,:,:,ivtx12),
     &         ngas,nspin)
        end do
      end if

      return
      end
