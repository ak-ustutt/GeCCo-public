*----------------------------------------------------------------------*
      subroutine condense_bc_info(
     &     cnt_info,self,
     &     iocc_op1, iocc_op2, iocc_op1op2, iocc_op1op2tmp,
     &     iocc_ex1,iocc_ex2,iocc_cnt,
     &     irst_op1, irst_op2, irst_op1op2, irst_op1op2tmp,
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
*----------------------------------------------------------------------*
      implicit none

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

      logical ::
     &     ok
      integer ::
     &     ijoin, ngas, nspin
      integer ::
     &     irst_ex1(2,orb_info%ngas,2,2,njoined_op1),
     &     irst_ex2(2,orb_info%ngas,2,2,njoined_op2),
     &     irst_cnt(2,orb_info%ngas,2,2,njoined_cnt),
     &     irst_cnt_dagger(2,orb_info%ngas,2,2,njoined_cnt),
     &     iocc_cnt_dagger(ngastp,2,njoined_cnt)
      integer, pointer ::
     &     igrph(:,:,:), ihpvgas(:,:)
      logical, external ::
     &     irestr_equal

c dbgmh
      print *,'m1'
c dbgend
      ngas = orb_info%ngas
      nspin = orb_info%nspin
      ihpvgas => orb_info%ihpvgas
      ! preliminary treatment of restrictions for EX, CNT:
c      call fit_restr(irst_ex1,iocc_ex1,
c     &     irst_op1,ihpvgas,ngas)
c      call fit_restr(irst_ex2,iocc_ex2,
c     &     irst_op2,ihpvgas,ngas)
c      call fit_restr(irst_cnt,iocc_cnt,
c     &     irst_op1,ihpvgas,ngas)
      ! QUICK FIX FOR FROZEN-CORE-R12
      do ijoin = 1, njoined_cnt
        iocc_cnt_dagger(1:ngastp,1:2,ijoin) =
     &       iocc_dagger(iocc_cnt(1:ngastp,1:2,ijoin))
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
      
        ok = .true.
        do ijoin = 1, njoined_cnt
          ok = ok.and.irestr_equal(irst_cnt(1,1,1,1,ijoin),.false.,
     &         irst_cnt_dagger(1,1,1,1,ijoin),.true.,ngas)
        end do
        if (.not.ok) then
          write(luout,*) 'generated restrictions: CNT != CNT^+'
          do ijoin = 1, njoined_cnt
            if (njoined_cnt.gt.1) write(luout,*) 'pair # ',ijoin
            call wrt_rstr(luout,irst_cnt(1,1,1,1,ijoin),ngas)
            call wrt_rstr(luout,irst_cnt_dagger(1,1,1,1,ijoin),ngas)
          end do
          call quit(1,'condense_bc_info','problem with restrictions !')
        end if
      end if

c      if (njoined_op1.eq.1) then
c        call fit_restr(irst_ex1,iocc_ex1,
c     &       irst_op1,ihpvgas,ngas)
c      else
c        call dummy_restr(irst_ex1,
c     &       iocc_ex1,njoined_op1,orb_info)
c      end if
c      if (njoined_op2.eq.1) then
c        call fit_restr(irst_ex2,iocc_ex2,
c     &       irst_op2,ihpvgas,ngas)
c      else
c        call dummy_restr(irst_ex2,
c     &       iocc_ex2,njoined_op2,orb_info)
c      end if
c      if (njoined_cnt.eq.1.and.njoined_op1.eq.1) then
c        call fit_restr(irst_cnt,iocc_cnt,
c     &     irst_op1,ihpvgas,ngas)
c      else if (njoined_cnt.eq.1.and.njoined_op2.eq.1) then
c        call fit_restr(irst_cnt,iocc_cnt,
c     &     irst_op2,ihpvgas,ngas)
c      else
c        call dummy_restr(irst_cnt,
c     &       iocc_cnt,njoined_cnt,orb_info)
c      end if
      ! end of preliminary code

      allocate(igrph(ngastp,2,
     &     max(njoined_op1,njoined_op2,njoined_op1op2,njoined_cnt)))
     
c dbgmh
      print *,'m2. cnt_info%cinfo_op1c: ',cnt_info%cinfo_op1c
      print *,'m2. (1,3): ',cnt_info%cinfo_op1c(1,3)
      print *,'m2. cnt_info%cinfo_op1a: ',cnt_info%cinfo_op1a
      print *,'m2. (1,3): ',cnt_info%cinfo_op1a(1,3)
c dbgend 
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

c dbgmh
      print *,'m3. cnt_info%cinfo_op1c: ',cnt_info%cinfo_op1c
      print *,'m3. (1,3): ',cnt_info%cinfo_op1c(1,3)
      print *,'m3. cnt_info%cinfo_op1a: ',cnt_info%cinfo_op1a
      print *,'m3. (1,3): ',cnt_info%cinfo_op1a(1,3)
c dbgend
      ! EX1
      call get_grph4occ(igrph,iocc_ex1,irst_ex1,njoined_op1,
     &                  str_info,orb_info,.true.)
c dbgmh
      print *,'m3a. cnt_info%cinfo_ex1c: ',cnt_info%cinfo_ex1c
c      print *,'m3a. (1,3): ',cnt_info%cinfo_ex1c(1,3)
      print *,'m3a. cnt_info%cinfo_ex1a: ',cnt_info%cinfo_ex1a
c      print *,'m3a. (1,3): ',cnt_info%cinfo_ex1a(1,3)
c dbgend
      call condense_occ(cnt_info%cinfo_ex1c, cnt_info%cinfo_ex1a,
     &                  cnt_info%cinfo_ex1c(1,3),
     &                  cnt_info%cinfo_ex1a(1,3),
     &                  iocc_ex1,njoined_op1,hpvxblkseq)
c dbgmh
      print *,'m3b'
c dbgend
      call condense_occ(cnt_info%cinfo_ex1c(1,2),
     &                  cnt_info%cinfo_ex1a(1,2),
     &                  cnt_info%cinfo_ex1c(1,3),
     &                  cnt_info%cinfo_ex1a(1,3),
     &                  igrph,njoined_op1,hpvxblkseq)

c dbgmh
      print *,'m4'
c dbgend
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

c dbgmh
      print *,'m5'
c dbgend
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

c dbgmh
      print *,'m6'
c dbgend
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

c dbgmh
      print *,'m7'
c dbgend
      ! OP1OP2
      call get_grph4occ(igrph,iocc_op1op2,irst_op1op2,njoined_op1op2,
     &                  str_info,orb_info,.true.)
c dbg
c      print *,'relevant call to condense_occ (OCC,OP1OP2)'
c dbg
      call condense_occ(cnt_info%cinfo_op1op2c, cnt_info%cinfo_op1op2a,
     &                  cnt_info%cinfo_op1op2c(1,3),
     &                  cnt_info%cinfo_op1op2a(1,3),
     &                  iocc_op1op2,njoined_op1op2,hpvxblkseq)
c dbg
c      print *,'relevant call to condense_occ (GRPH,OP1OP2)'
c dbg
      call condense_occ(cnt_info%cinfo_op1op2c(1,2),
     &                  cnt_info%cinfo_op1op2a(1,2),
     &                  cnt_info%cinfo_op1op2c(1,3),
     &                  cnt_info%cinfo_op1op2a(1,3),
     &                  igrph,njoined_op1op2,hpvxblkseq)

c dbgmh
      print *,'m8'
c dbgend
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
c dbg
c      print *,'relevant call to set_mapping_info C'
c dbg
      ! EX1/EX2 for C
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  1,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_op1op2tmp,merge_map12,
     &                                  njoined_op1op2,hpvxblkseq)
c dbg
c      print *,'relevant call to set_mapping_info A'
c dbg
      ! EX2/EX1 for A
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  2,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_op1op2tmp,merge_map21,
     &                                  njoined_op1op2,hpvxblkseq)

c dbgmh
      print *,'mend'
c dbgend
      return
      end
