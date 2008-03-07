*----------------------------------------------------------------------*
      subroutine condense_bc_info(
     &     cnt_info,
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

      integer ::
     &     ijoin, ngas
      integer ::
     &     irst_ex1(2,orb_info%ngas,2,2,njoined_op1),
     &     irst_ex2(2,orb_info%ngas,2,2,njoined_op2),
     &     irst_cnt(2,orb_info%ngas,2,2,njoined_cnt),
     &     iocc_cnt_dagger(ngastp,2,njoined_cnt)
      integer, pointer ::
     &     igrph(:,:,:), ihpvgas(:,:)


      ngas = orb_info%ngas
      ihpvgas => orb_info%ihpvgas
      ! preliminary treatment of restrictions for EX, CNT:
c      call fit_restr(irst_ex1,iocc_ex1,
c     &     irst_op1,ihpvgas,ngas)
c      call fit_restr(irst_ex2,iocc_ex2,
c     &     irst_op2,ihpvgas,ngas)
c      call fit_restr(irst_cnt,iocc_cnt,
c     &     irst_op1,ihpvgas,ngas)
      call dummy_restr(irst_ex1,
     &     iocc_ex1,njoined_op1,orb_info)
      call dummy_restr(irst_ex2,
     &     iocc_ex2,njoined_op2,orb_info)
      call dummy_restr(irst_cnt,
     &     iocc_cnt,njoined_cnt,orb_info)
      ! end of preliminary code

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
c dbg
c      print *,'OP1:'
c      call wrt_occ_n(6,iocc_op1,njoined_op1)
c      print *,' C: ',cnt_info%cinfo_op1c(1:nca_blk(1,1),1)
c      print *,'    ',cnt_info%cinfo_op1c(1:nca_blk(1,1),3)
c      print *,' A: ',cnt_info%cinfo_op1a(1:nca_blk(2,1),1)
c      print *,'    ',cnt_info%cinfo_op1a(1:nca_blk(2,1),3)
c dbg

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
c dbg
c      print *,'OP2:'
c      call wrt_occ_n(6,iocc_op2,njoined_op2)
c      print *,' C: ',cnt_info%cinfo_op2c(1:nca_blk(1,2),1)
c      print *,'    ',cnt_info%cinfo_op2c(1:nca_blk(1,2),3)
c      print *,' A: ',cnt_info%cinfo_op2a(1:nca_blk(2,2),1)
c      print *,'    ',cnt_info%cinfo_op2a(1:nca_blk(2,2),3)
c dbg

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
c dbg
c      print *,'call for OP1OP2'
c dbg
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
c dbg
c      print *,'call for ex1ex2'
c dbg
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  1,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_op1op2tmp,merge_map12,
     &                                  njoined_op1op2,hpvxblkseq)
      ! EX2/EX1 for A
c dbg
c      print *,'call for ex2ex1'
c dbg
      call set_mapping_info(cnt_info%map_info_12c,cnt_info%map_info_12a,
     &                  2,
     &                  iocc_ex2,njoined_op2,.false.,
     &                  iocc_ex1,njoined_op1,.false.,
     &                  iocc_op1op2tmp,merge_map21,
     &                                  njoined_op1op2,hpvxblkseq)

      return
      end
