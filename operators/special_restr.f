*----------------------------------------------------------------------*
      subroutine special_restr(irestr_op1op2,
     &     iocc_op1op2,njoined_op1op2,
     &     merge_op1op2,
     &     iocc_op1,iocc_ex1,irestr_op1,njoined_op1,
     &     iocc_op2,iocc_ex2,irestr_op2,njoined_op2,
     &     nspin,ngas)
*----------------------------------------------------------------------*
*     quick and dirty fix to pass restrictions in selected cases
*     (that in fact only occur for fixed-geminal intermediates with
*      frozen core orbitals)
*     if the C or A part of OP1OP2 is exactly the C or A part of 
*     OP1 and OP2 --> inherit their restriction
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
     &     irestr_op2(2,ngas,2,2,nspin,njoined_op2)

      integer ::
     &     ivtx12, ivtx1, ivtx2, idx, ica
      logical, external ::
     &     list_cmp
     
      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'special_restr')
        print *,'on input: 12'
        do ivtx12 = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op1op2,irestr_op1op2,
     &         ngas,nspin)
        end do
        print *,' 1'
        do ivtx12 = 1, njoined_op1
          call wrt_occ_rstr(luout,ivtx12,iocc_op1,irestr_op1,
     &         ngas,nspin)
        end do
        print *,' 1ex'
        do ivtx12 = 1, njoined_op1
          call wrt_occ_rstr(luout,ivtx12,iocc_ex1,irestr_op1,
     &         ngas,nspin)
        end do
        print *,' 2'
        do ivtx12 = 1, njoined_op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op2,irestr_op2,
     &         ngas,nspin)
        end do
        print *,' 2ex'
        do ivtx12 = 1, njoined_op2
          call wrt_occ_rstr(luout,ivtx12,iocc_ex2,irestr_op2,
     &         ngas,nspin)
        end do
      end if

      idx = 1
      do ivtx12 = 1, njoined_op1op2
        ivtx1 = merge_op1op2(idx)
        ivtx2 = merge_op1op2(idx+ivtx1+1)
        idx = idx + ivtx1 + ivtx2 + 2
        if (ivtx1.gt.0) then
          do ica = 1, 2
            if (list_cmp(iocc_op1(1,ica,ivtx1),
     &                  iocc_ex1(1,ica,ivtx1),ngastp) .and.
     &          list_cmp(iocc_op1op2(1,ica,ivtx12),
     &                  iocc_ex1(1,ica,ivtx1),ngastp))
     &      then
              irestr_op1op2(1:2,1:ngas,ica,1:2,1:nspin,ivtx12) =
     &             irestr_op1(1:2,1:ngas,ica,1:2,1:nspin,ivtx1)
            end if
          end do
        end if
        if (ivtx2.gt.0) then
          do ica = 1, 2
            if (list_cmp(iocc_op2(1,ica,ivtx2),
     &                   iocc_ex2(1,ica,ivtx2),ngastp) .and.
     &          list_cmp(iocc_op1op2(1,ica,ivtx12),
     &                   iocc_ex2(1,ica,ivtx2),ngastp))
     &      then
              irestr_op1op2(1:2,1:ngas,ica,1:2,1:nspin,ivtx12) =
     &             irestr_op2(1:2,1:ngas,ica,1:2,1:nspin,ivtx2)
            end if
          end do
        end if

      end do

      if (ntest.ge.100) then
        print *,'on output:'
        do ivtx12 = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ivtx12,iocc_op1op2,irestr_op1op2,
     &         ngas,nspin)
        end do
      end if

      return
      end
