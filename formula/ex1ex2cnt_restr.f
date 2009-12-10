*----------------------------------------------------------------------*
      subroutine ex1ex2cnt_restr(
     &     rst_ex1, rst_ex2, rst_cnt,
     &     rst_cnt_in,
     &     self,
     &     occ_op1, occ_op2, occ_op1op2, occ_op1op2tmp,
     &     occ_ex1, occ_ex2, occ_cnt,
     &     rst_op1, rst_op2, rst_op1op2, rst_op1op2tmp,
     &     merge_map1, merge_map2, merge_map12, merge_map21,
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*
*     obtain restrictions for (EX1,CNT), (EX2,CNT)
*     based on restrictions on OP1, OP2, OP1OP2 and possibly
*     a prototype restriction on rst_cnt_in
*
*     exracted from condense_bc_info
*
*     dec 2009, andreas
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'hpvxseq.h'
      include 'def_orbinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_operators.h'
      include 'def_contraction_info.h'

      type(strinf) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      logical, intent(in) ::
     &     self
      integer, intent(in) ::
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     occ_op1(ngastp,2,njoined_op1),
     &     occ_op2(ngastp,2,njoined_op2),
     &     occ_op1op2(ngastp,2,njoined_op1op2),
     &     occ_op1op2tmp(ngastp,2,njoined_op1op2),
     &     occ_ex1(ngastp,2,njoined_op1),
     &     occ_ex2(ngastp,2,njoined_op2),
     &     occ_cnt(ngastp,2,njoined_cnt),
     &     rst_op1(2,orb_info%ngas,2,2,njoined_op1),
     &     rst_op2(2,orb_info%ngas,2,2,njoined_op2),
     &     rst_op1op2(2,orb_info%ngas,2,2,njoined_op1op2),
     &     rst_op1op2tmp(2,orb_info%ngas,2,2,njoined_op1op2),
     &     rst_cnt_in(2,orb_info%ngas,2,2,njoined_cnt),
     &     merge_map1(*), merge_map2(*), merge_map12(*), merge_map21(*)

      integer, intent(out) ::
     &     rst_ex1(2,orb_info%ngas,2,2,njoined_op1),
     &     rst_ex2(2,orb_info%ngas,2,2,njoined_op2),
     &     rst_cnt(2,orb_info%ngas,2,2,njoined_cnt)

      logical ::
     &     ok, equal
      integer ::
     &     ijoin, ngas, nspin, ii
      integer ::
     &     rst_common(2,orb_info%ngas,2,2),
     &     rst_cnt_dagger(2,orb_info%ngas,2,2,njoined_cnt),
     &     occ_cnt_dagger(ngastp,2,njoined_cnt)
      integer, pointer ::
     &     ihpvgas(:,:)
      logical, external ::
     &     irestr_equal, common_restr

      ngas = orb_info%ngas
      nspin = orb_info%nspin
      ihpvgas => orb_info%ihpvgas

      do ijoin = 1, njoined_cnt
        occ_cnt_dagger(1:ngastp,1:2,ijoin) =
     &         iocc_dagger(occ_cnt(1:ngastp,1:2,ijoin))
      end do

      call fit_restr3(rst_ex1,rst_cnt,
     &     occ_op1,occ_ex1,occ_cnt,
     &     rst_op1,merge_map1,
     &     njoined_op1,njoined_cnt,
     &     ihpvgas,ngas,nspin)
      if (.not.self) then
        call fit_restr3(rst_ex2,rst_cnt_dagger,
     &       occ_op2,occ_ex2,occ_cnt_dagger,
     &       rst_op2,merge_map2,
     &       njoined_op2,njoined_cnt,
     &       ihpvgas,ngas,nspin)

        call fit_restr3(rst_ex2,rst_ex1,
     &       occ_op1op2tmp,occ_ex2,occ_ex1,
     &       rst_op1op2tmp,merge_map12,
     &       njoined_op1op2,njoined_op2,
     &       ihpvgas,ngas,nspin)
      else
        do ijoin = 1, njoined_cnt/2
          rst_cnt_dagger(1:2,1:ngas,1:2,1:2,ijoin) =
     &         rst_cnt(1:2,1:ngas,1:2,1:2,ijoin+njoined_cnt/2)
        end do
        do ijoin = 1, njoined_cnt/2
          rst_cnt_dagger(1:2,1:ngas,1:2,1:2,ijoin+njoined_cnt/2) =
     &         rst_cnt(1:2,1:ngas,1:2,1:2,ijoin)
        end do
      end if
      
      ok = .true.
      do ijoin = 1, njoined_cnt
        equal = irestr_equal(rst_cnt(1,1,1,1,ijoin),.false.,
     &         rst_cnt_dagger(1,1,1,1,ijoin),.true.,ngas)
        ok = ok.and.equal

        ! use the maximum overlap of restrictions
        if (.not.equal) then
          if (.not.common_restr(rst_common,
     &           rst_cnt(1,1,1,1,ijoin),.false.,
     &           rst_cnt_dagger(1,1,1,1,ijoin),.true.,
     &           ihpvgas,ngas,nspin)) then
            do ii = 1, njoined_cnt
              call wrt_occ_rstr(luout,ii,occ_cnt(:,:,ii),
     &                                     rst_cnt(:,:,:,:,ijoin),
     &                                     ngas,nspin)
              call wrt_occ_rstr(luout,ii,occ_cnt_dagger(:,:,ii),
     &                                   rst_cnt_dagger(:,:,:,:,ijoin),
     &                                     ngas,nspin)
            end do
            call quit(1,'ex1ex2cnt_restr',
     &             'cannot merge restrictions')
          end if
          ! set rst_cnt to maximum common restriction
          rst_cnt(:,:,:,:,ijoin) = rst_common(:,:,:,:)
          ! rst_cnt_dagger is not used further!
        end if
      end do

      if (rst_cnt_in(1,1,1,1,1).gt.-1) then

        ! restrict to input restriction
        do ijoin = 1, njoined_cnt
          equal = irestr_equal(rst_cnt(1,1,1,1,ijoin),.false.,
     &                        rst_cnt_in(1,1,1,1,ijoin),.false.,ngas)
          ok = ok.and.equal

          ! use the maximum overlap of restrictions
          if (.not.equal) then
            if (.not.common_restr(rst_common,
     &           rst_cnt(1,1,1,1,ijoin),.false.,
     &           rst_cnt_in(1,1,1,1,ijoin),.false.,
     &             ihpvgas,ngas,nspin)) then
              do ii = 1, njoined_cnt
                call wrt_occ_rstr(luout,ii,occ_cnt(:,:,ii),
     &                                     rst_cnt(:,:,:,:,ijoin),
     &                                     ngas,nspin)
                call wrt_occ_rstr(luout,ii,occ_cnt(:,:,ii),
     &                                   rst_cnt_in(:,:,:,:,ijoin),
     &                                     ngas,nspin)
              end do
              call quit(1,'ex1ex2cnt_restr',
     &               'cannot merge restrictions (2)')
            end if
            ! set rst_cnt to maximum common restriction
            rst_cnt(:,:,:,:,ijoin) = rst_common(:,:,:,:)
          end if
        end do

      end if

      if (ntest.ge.100) then
        write(luout,*) 'generated restrictions: CNT'
        do ijoin = 1, njoined_cnt
          if (njoined_cnt.gt.1) write(luout,*) 'pair # ',ijoin
          call wrt_rstr(luout,rst_cnt(1,1,1,1,ijoin),ngas)
          call wrt_rstr(luout,rst_cnt_dagger(1,1,1,1,ijoin),ngas)
        end do
        write(luout,*) 'OP1 is'
        do ijoin = 1, njoined_op1
          call wrt_occ_rstr(luout,ijoin,occ_op1(:,:,ijoin),
     &         rst_op1(:,:,:,:,ijoin),
     &         ngas,nspin)
        end do
        write(luout,*) 'EX1 is'
        do ijoin = 1, njoined_op1
          call wrt_occ_rstr(luout,ijoin,occ_ex1(:,:,ijoin),
     &         rst_ex1(:,:,:,:,ijoin),
     &         ngas,nspin)
        end do
        if (.not.self) then
          write(luout,*) 'OP2 is'
          do ijoin = 1, njoined_op2
            call wrt_occ_rstr(luout,ijoin,occ_op2(:,:,ijoin),
     &                                    rst_op2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do
          write(luout,*) 'EX2 is'
          do ijoin = 1, njoined_op2
            call wrt_occ_rstr(luout,ijoin,occ_ex2(:,:,ijoin),
     &                                    rst_ex2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do
        end if
        write(luout,*) 'OP1OP2TMP is'
        do ijoin = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ijoin,occ_op1op2tmp(:,:,ijoin),
     &                                    rst_op1op2tmp(:,:,:,:,ijoin),
     &                                    ngas,nspin)
        end do
        write(luout,*) 'OP1OP2 is'
        do ijoin = 1, njoined_op1op2
          call wrt_occ_rstr(luout,ijoin,occ_op1op2(:,:,ijoin),
     &                                    rst_op1op2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
        end do
      end if

      return
      end
