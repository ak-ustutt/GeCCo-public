*----------------------------------------------------------------------*
      subroutine condense_bc_info(
     &     cnt_info,self,
     &     iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
     &     irst_ex1_io, irst_ex2_io, irst_cnt_io, mode_rst,
     &     merge_map1, merge_map2, merge_map12, merge_map21,
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     obtain contraction info in condensed form
*     
*      OP1: cinfo_op1c/a(1..nc/a,1) -> condensed occupations
*           cinfo_op1c/a(1..nc/a,2) -> condensed graph indices
*           cinfo_op1c/a(1..nc/a,2) -> condensed hpvx info
*      dto. for OP2, OP1OP2, EX1, EX2, CNT
*
*      mapping info (cf. set_mapping_info):
*      OP1 -> EX1/CNT :  map_info_1c/a 
*      OP2 -> EX2/CNT :  map_info_2c/a
*      EX1/EX2 -> OP1OP2 : map_info_12c/a
*
*     mode_rst:  0 -- ignore rst_cnt_io
*                    1 -- return rst_cnt on rst_cnt_io
*                    2 -- set rst_cnt to rst_cnt_io
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

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
      type(contraction_info), intent(inout) ::
     &     cnt_info
      logical, intent(in) ::
     &     self
      integer, intent(in) ::
     &     mode_rst,
     &     njoined_op1, njoined_op2, njoined_op1op2, njoined_cnt,
     &     iocc_op1(ngastp,2,njoined_op1),
     &     iocc_op2(ngastp,2,njoined_op2),
     &     iocc_op1op2(ngastp,2,njoined_op1op2),
     &     iocc_op1op2tmp(ngastp,2,njoined_op1op2),
     &     iocc_ex1(ngastp,2,njoined_op1),
     &     iocc_ex2(ngastp,2,njoined_op2),
     &     iocc_cnt(ngastp,2,njoined_cnt),
     &     irst_op1(2,orb_info%ngas,2,2,njoined_op1),
     &     irst_op2(2,orb_info%ngas,2,2,njoined_op2),
     &     irst_op1op2(2,orb_info%ngas,2,2,njoined_op1op2),
     &     irst_op1op2tmp(2,orb_info%ngas,2,2,njoined_op1op2),
     &     merge_map1(*), merge_map2(*), merge_map12(*), merge_map21(*)

      integer, intent(inout) ::
     &     irst_ex1_io(2,orb_info%ngas,2,2,njoined_op1),
     &     irst_ex2_io(2,orb_info%ngas,2,2,njoined_op2),
     &     irst_cnt_io(2,orb_info%ngas,2,2,njoined_cnt)

      logical ::
     &     ok, equal
      integer ::
     &     ijoin, ngas, nspin, ii
      integer ::
     &     irst_ex1(2,orb_info%ngas,2,2,njoined_op1),
     &     irst_ex2(2,orb_info%ngas,2,2,njoined_op2),
     &     irst_cnt(2,orb_info%ngas,2,2,njoined_cnt),
     &     irst_common(2,orb_info%ngas,2,2),
     &     irst_cnt_dagger(2,orb_info%ngas,2,2,njoined_cnt),
     &     iocc_cnt_dagger(ngastp,2,njoined_cnt)
      integer, pointer ::
     &     igrph(:,:,:), ihpvgas(:,:)
      logical, external ::
     &     irestr_equal, common_restr

      ngas = orb_info%ngas
      nspin = orb_info%nspin
      ihpvgas => orb_info%ihpvgas

      if (mode_rst.ne.2) then

        do ijoin = 1, njoined_cnt
          iocc_cnt_dagger(1:ngastp,1:2,ijoin) =
     &         iocc_dagger(iocc_cnt(1:ngastp,1:2,ijoin))
        end do

        call fit_restr3(irst_ex1,irst_cnt,
     &     iocc_op1,iocc_ex1,iocc_cnt,
     &     irst_op1,merge_map1,
     &     njoined_op1,njoined_cnt,
     &     ihpvgas,ngas,nspin)
        if (.not.self) then
          call fit_restr3(irst_ex2,irst_cnt_dagger,
     &       iocc_op2,iocc_ex2,iocc_cnt_dagger,
     &       irst_op2,merge_map2,
     &       njoined_op2,njoined_cnt,
     &       ihpvgas,ngas,nspin)

c test
          write(luout,*) 'EX1 before'
          do ijoin = 1, njoined_op1
            call wrt_occ_rstr(luout,ijoin,iocc_ex1(:,:,ijoin),
     &                                    irst_ex1(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do
          write(luout,*) 'EX2 before'
          do ijoin = 1, njoined_op2
            call wrt_occ_rstr(luout,ijoin,iocc_ex2(:,:,ijoin),
     &                                    irst_ex2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do

          call fit_restr3(irst_ex2,irst_ex1,
     &       iocc_op1op2tmp,iocc_ex2,iocc_ex1,
     &       irst_op1op2tmp,merge_map12,
     &       njoined_op1op2,njoined_op2,
     &       ihpvgas,ngas,nspin)

          write(luout,*) 'EX1 after'
          do ijoin = 1, njoined_op1
            call wrt_occ_rstr(luout,ijoin,iocc_ex1(:,:,ijoin),
     &                                    irst_ex1(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do
          write(luout,*) 'EX2 after'
          do ijoin = 1, njoined_op2
            call wrt_occ_rstr(luout,ijoin,iocc_ex2(:,:,ijoin),
     &                                    irst_ex2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
          end do
c test end

      
          ok = .true.
          do ijoin = 1, njoined_cnt
            equal = irestr_equal(irst_cnt(1,1,1,1,ijoin),.false.,
     &         irst_cnt_dagger(1,1,1,1,ijoin),.true.,ngas)
            ok = ok.and.equal
c patch
            ! use the maximum overlap of restrictions
            if (.not.equal) then
              if (.not.common_restr(irst_common,
     &           irst_cnt(1,1,1,1,ijoin),.false.,
     &           irst_cnt_dagger(1,1,1,1,ijoin),.true.,
     &           ihpvgas,ngas,nspin)) then
                do ii = 1, njoined_cnt
                  call wrt_occ_rstr(luout,ii,iocc_cnt(:,:,ii),
     &                                     irst_cnt(:,:,:,:,ijoin),
     &                                     ngas,nspin)
                  call wrt_occ_rstr(luout,ii,iocc_cnt_dagger(:,:,ii),
     &                                   irst_cnt_dagger(:,:,:,:,ijoin),
     &                                     ngas,nspin)
                end do
                call quit(1,'condense_bc_info',
     &             'cannot merge restrictions')
              end if
              ! set irst_cnt to maximum common restriction
              irst_cnt(:,:,:,:,ijoin) = irst_common(:,:,:,:)
              ! irst_cnt_dagger is not used further!
            end if
c end patch
          end do
          if (ntest.ge.100) then
            write(luout,*) 'generated restrictions: CNT'
            do ijoin = 1, njoined_cnt
              if (njoined_cnt.gt.1) write(luout,*) 'pair # ',ijoin
              call wrt_rstr(luout,irst_cnt(1,1,1,1,ijoin),ngas)
              call wrt_rstr(luout,irst_cnt_dagger(1,1,1,1,ijoin),ngas)
            end do
            write(luout,*) 'OP1 is'
            do ijoin = 1, njoined_op1
              call wrt_occ_rstr(luout,ijoin,iocc_op1(:,:,ijoin),
     &                                    irst_op1(:,:,:,:,ijoin),
     &                                    ngas,nspin)
            end do
            write(luout,*) 'EX1 is'
            do ijoin = 1, njoined_op1
              call wrt_occ_rstr(luout,ijoin,iocc_ex1(:,:,ijoin),
     &                                    irst_ex1(:,:,:,:,ijoin),
     &                                    ngas,nspin)
            end do
            if (.not.self) then
              write(luout,*) 'OP2 is'
              do ijoin = 1, njoined_op2
                call wrt_occ_rstr(luout,ijoin,iocc_op2(:,:,ijoin),
     &                                    irst_op2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
              end do
              write(luout,*) 'EX2 is'
              do ijoin = 1, njoined_op2
                call wrt_occ_rstr(luout,ijoin,iocc_ex2(:,:,ijoin),
     &                                    irst_ex2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
              end do
            end if
            write(luout,*) 'OP1OP2TMP is'
            do ijoin = 1, njoined_op1op2
              call wrt_occ_rstr(luout,ijoin,iocc_op1op2tmp(:,:,ijoin),
     &                                    irst_op1op2tmp(:,:,:,:,ijoin),
     &                                    ngas,nspin)
            end do
            write(luout,*) 'OP1OP2 is'
            do ijoin = 1, njoined_op1op2
              call wrt_occ_rstr(luout,ijoin,iocc_op1op2(:,:,ijoin),
     &                                    irst_op1op2(:,:,:,:,ijoin),
     &                                    ngas,nspin)
            end do
c          call quit(1,'condense_bc_info','problem with restrictions !')
          end if
        end if                  ! .not.self
      end if                    ! mode_cnt_rstr.ne.2

      if (mode_rst.eq.1) then ! return what has been set
        irst_ex1_io = irst_ex1
        irst_ex2_io = irst_ex2
        irst_cnt_io = irst_cnt
      end if
      if (mode_rst.eq.2) then ! set to input restrictions
        irst_ex1 = irst_ex1_io
        irst_ex2 = irst_ex2_io
        irst_cnt = irst_cnt_io
      end if

      allocate(igrph(ngastp,2,
     &     max(njoined_op1,njoined_op2,njoined_op1op2,njoined_cnt)))
     
      ! OP1
      call get_grph4occ(igrph,iocc_op1,irst_op1,njoined_op1,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_op1c, cnt_info%cinfo_op1a,
     &                  cnt_info%cinfo_op1c(1,3),
     &                  cnt_info%cinfo_op1a(1,3),
     &                  iocc_op1,njoined_op1,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_op1c(1,2),
     &                  cnt_info%cinfo_op1a(1,2),
     &                  cnt_info%cinfo_op1c(1,3),
     &                  cnt_info%cinfo_op1a(1,3),
     &                  igrph,njoined_op1,hpvxblkseq)

      ! EX1
      call get_grph4occ(igrph,iocc_ex1,irst_ex1,njoined_op1,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_ex1c, cnt_info%cinfo_ex1a,
     &                  cnt_info%cinfo_ex1c(1,3),
     &                  cnt_info%cinfo_ex1a(1,3),
     &                  iocc_ex1,njoined_op1,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_ex1c(1,2),
     &                  cnt_info%cinfo_ex1a(1,2),
     &                  cnt_info%cinfo_ex1c(1,3),
     &                  cnt_info%cinfo_ex1a(1,3),
     &                  igrph,njoined_op1,hpvxblkseq)

      ! OP2
      call get_grph4occ(igrph,iocc_op2,irst_op2,njoined_op2,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_op2c, cnt_info%cinfo_op2a,
     &                  cnt_info%cinfo_op2c(1,3),
     &                  cnt_info%cinfo_op2a(1,3),
     &                  iocc_op2,njoined_op2,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_op2c(1,2),
     &                  cnt_info%cinfo_op2a(1,2),
     &                  cnt_info%cinfo_op2c(1,3),
     &                  cnt_info%cinfo_op2a(1,3),
     &                  igrph,njoined_op2,hpvxblkseq)

      ! EX2
      call get_grph4occ(igrph,iocc_ex2,irst_ex2,njoined_op2,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_ex2c, cnt_info%cinfo_ex2a,
     &                  cnt_info%cinfo_ex2c(1,3),
     &                  cnt_info%cinfo_ex2a(1,3),
     &                  iocc_ex2,njoined_op2,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_ex2c(1,2),
     &                  cnt_info%cinfo_ex2a(1,2),
     &                  cnt_info%cinfo_ex2c(1,3),
     &                  cnt_info%cinfo_ex2a(1,3),
     &                  igrph,njoined_op2,hpvxblkseq)

      ! CNT
      call get_grph4occ(igrph,iocc_cnt,irst_cnt,njoined_cnt,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_cntc,
     &                  cnt_info%cinfo_cnta,
     &                  cnt_info%cinfo_cntc(1,3),
     &                  cnt_info%cinfo_cnta(1,3),
     &                  iocc_cnt,njoined_cnt,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_cntc(1,2),
     &                  cnt_info%cinfo_cnta(1,2),
     &                  cnt_info%cinfo_cntc(1,3),
     &                  cnt_info%cinfo_cnta(1,3),
     &                  igrph,njoined_cnt,hpvxblkseq)

      ! OP1OP2
      call get_grph4occ(igrph,iocc_op1op2,irst_op1op2,njoined_op1op2,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_op1op2c, cnt_info%cinfo_op1op2a,
     &                  cnt_info%cinfo_op1op2c(1,3),
     &                  cnt_info%cinfo_op1op2a(1,3),
     &                  iocc_op1op2,njoined_op1op2,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_op1op2c(1,2),
     &                  cnt_info%cinfo_op1op2a(1,2),
     &                  cnt_info%cinfo_op1op2c(1,3),
     &                  cnt_info%cinfo_op1op2a(1,3),
     &                  igrph,njoined_op1op2,hpvxblkseq)

      call get_grph4occ(igrph,iocc_op1op2tmp,irst_op1op2tmp,
     &                                       njoined_op1op2,
     &                  str_info,orb_info,.true.)
      call condense_occ(cnt_info%cinfo_op1op2tmpc,
     &                  cnt_info%cinfo_op1op2tmpa,
     &                  cnt_info%cinfo_op1op2tmpc(1,3),
     &                  cnt_info%cinfo_op1op2tmpa(1,3),
     &                  iocc_op1op2tmp,njoined_op1op2,hpvxblkseq)
      call condense_occ(cnt_info%cinfo_op1op2tmpc(1,2),
     &                  cnt_info%cinfo_op1op2tmpa(1,2),
     &                  cnt_info%cinfo_op1op2tmpc(1,3),
     &                  cnt_info%cinfo_op1op2tmpa(1,3),
     &                  igrph,njoined_op1op2,hpvxblkseq)

      deallocate(igrph)

      ! OP1 -> CNT/EX1 map
      call set_mapping_info(cnt_info%map_info_1c,cnt_info%map_info_1a,
     &                  0,
     &                  iocc_cnt,njoined_cnt,.false.,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_op1,merge_map1,njoined_op1,hpvxblkseq)
      ! OP2 -> CNT^+/EX2 map
      call set_mapping_info(cnt_info%map_info_2c,cnt_info%map_info_2a,
     &                  0,
     &                  iocc_cnt,njoined_cnt,.true.,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_op2,merge_map2,njoined_op2,hpvxblkseq)

      ! OP1OP2 -> EX1/EX2 map
      ! note: map goes to OP1OP2tmp, actually
      !  if OP1OP2 differs, this is taken care of by the additional
      !  reordering step
      ! EX1/EX2 for C
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  1,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_op1op2tmp,merge_map12,
     &                                  njoined_op1op2,hpvxblkseq)
      ! EX2/EX1 for A
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  2,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_op1op2tmp,merge_map21,
     &                                  njoined_op1op2,hpvxblkseq)

      return
      end
