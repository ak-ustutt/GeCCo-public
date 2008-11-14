*----------------------------------------------------------------------*
      logical function cmp_contr(contr1,contr2,ignore_fac)
*----------------------------------------------------------------------*
*     compare contractions
*     we assume canonical ordering, so the contractions are equivalent
*     if all entries are equivalent
*     for some cases, we might ignore the factor
*     NEW: no canonical ordering necessary!
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr1, contr2
      logical, intent(in) ::
     &     ignore_fac

      integer ::
     &     iarc, ivtx, nvtx, nj,
     &     sum_blk1, sum_blk2, sum_op1, sum_op2,
     &     occ(ngastp,2)
      type(cntr_vtx), pointer ::
     &     vtx1(:), vtx2(:)
      type(cntr_arc), pointer ::
     &     arc1(:), arc2(:)

      integer, pointer ::
     &     scr(:)
      integer(8), pointer ::
     &     ivtx1(:),topo1(:,:),xlines1(:,:),
     &     ivtx2(:),topo2(:,:),xlines2(:,:)

      integer, external ::
     &     list_cmp, i8list_cmp, njres_contr
      logical, external ::
     &     iocc_zero

c dbg
c      if (contr1%idx_res.eq.15) print *,'comparing'
c dbg
      cmp_contr = contr1%nvtx.eq.contr2%nvtx .and.
     &            contr1%narc.eq.contr2%narc .and.
     &            contr1%nxarc.eq.contr2%nxarc .and.
     &            contr1%nsupvtx.eq.contr2%nsupvtx
      cmp_contr = cmp_contr.and.
     &            contr1%idx_res.eq.contr2%idx_res .and.
     &            contr1%iblk_res.eq.contr2%iblk_res

      if (.not.cmp_contr) return

      if (.not.ignore_fac)
     &     cmp_contr = abs(contr1%fac-contr2%fac).lt.1d-12

      if (.not.cmp_contr) return
c dbg
c      if (contr1%idx_res.eq.15) print *,'survived basic tests'
c dbg

      ! as the exact comparison takes some time, pre-screen with
      ! checksum
      vtx1 => contr1%vertex
      vtx2 => contr2%vertex
      
      sum_op1  = 0
      sum_blk1 = 0
      do ivtx = 1, contr1%nvtx
        sum_op1  = sum_op1  + vtx1(ivtx)%idx_op
        sum_blk1 = sum_blk1 + vtx1(ivtx)%iblk_op
      end do

      sum_op2  = 0
      sum_blk2 = 0
      do ivtx = 1, contr2%nvtx
        sum_op2  = sum_op2  + vtx2(ivtx)%idx_op
        sum_blk2 = sum_blk2 + vtx2(ivtx)%iblk_op
      end do

      cmp_contr = sum_op1.eq.sum_op2.and.sum_blk1.eq.sum_blk2
c dbg
c      if (contr1%idx_res.eq.15.and..not.cmp_contr) then
c        print *,'vtx1: ',vtx1(1:contr1%nvtx)%idx_op
c        print *,'vtx2: ',vtx2(1:contr2%nvtx)%idx_op
c        print *,'sum_op1, sum_op2: ',sum_op1,sum_op2
c        print *,'sum_blk1, sum_blk2: ',sum_blk1,sum_blk2
c      end if
c dbg

      if (.not.cmp_contr) return

      arc1 => contr1%arc
      arc2 => contr2%arc

      occ = 0
      do iarc = 1, contr1%narc
        occ = occ + arc1(iarc)%occ_cnt
      end do
      do iarc = 1, contr2%narc
        occ = occ - arc2(iarc)%occ_cnt
      end do

      cmp_contr = iocc_zero(occ)
c dbg
c      if (contr1%idx_res.eq.15.and..not.cmp_contr) then
c        call wrt_occ(6,occ)
c      end if
c dbg

      if (.not.cmp_contr) return
c dbg
c      print *,'survived initial tests'
c dbg

      ! NEW:
      nj = njres_contr(contr1)  ! idxres was already compared
      nvtx = contr1%nvtx        ! nvtx dto.
      allocate(ivtx1(nvtx),topo1(nvtx,nvtx),xlines1(nvtx,nj),
     &         ivtx2(nvtx),topo2(nvtx,nvtx),xlines2(nvtx,nj),
     &         scr(nvtx))
      call pack_contr(scr,ivtx1,topo1,xlines1,contr1,nj)
      call pack_contr(scr,ivtx2,topo2,xlines2,contr2,nj)

      call topo_make_unique(scr,ivtx1,topo1,xlines1,nvtx,nj)
      call topo_make_unique(scr,ivtx2,topo2,xlines2,nvtx,nj)

      cmp_contr = i8list_cmp(ivtx1,ivtx2,nvtx).eq.0
      if (ntest.ge.100) write(luout,*) 'cmp_contr > (1): ',cmp_contr
      cmp_contr = cmp_contr.and.
     &            i8list_cmp(xlines1,xlines2,nvtx*nj).eq.0
      if (ntest.ge.100) write(luout,*) 'cmp_contr > (2): ',cmp_contr
      cmp_contr = cmp_contr.and.
     &            i8list_cmp(topo1,topo2,nvtx*nvtx).eq.0
      if (ntest.ge.100) write(luout,*) 'cmp_contr > (3): ',cmp_contr
c dbg
c      if (.not.cmp_contr.and.contr1%idx_res.eq.15) then
c        print *,'topo1'
c        call prt_contr_p(6,ivtx1,ivtx1,topo1,xlines1,nvtx,nj)
c        print *,'topo2'
c        call prt_contr_p(6,ivtx2,ivtx2,topo2,xlines2,nvtx,nj)
c      end if
c dbg

      deallocate(ivtx1,topo1,xlines1,
     &           ivtx2,topo2,xlines2,scr)

      return
      ! OLD:
      do ivtx = 1, contr1%nsupvtx
        cmp_contr = cmp_contr.and.
     &       contr1%svertex(ivtx).eq.contr2%svertex(ivtx)
      end do

      if (.not.cmp_contr) return
      
      vtx1 => contr1%vertex
      vtx2 => contr2%vertex
      
      do ivtx = 1, contr1%nvtx
        cmp_contr = cmp_contr.and.
     &       vtx1(ivtx)%idx_op .eq.vtx2(ivtx)%idx_op .and.
     &       vtx1(ivtx)%iblk_op.eq.vtx2(ivtx)%iblk_op
      end do

      if (.not.cmp_contr) return

      arc1 => contr1%arc
      arc2 => contr2%arc

      do iarc = 1, contr1%narc
        cmp_contr = cmp_contr.and.
     &       arc1(iarc)%link(1).eq.arc2(iarc)%link(1) .and.
     &       arc1(iarc)%link(2).eq.arc2(iarc)%link(2) .and.
     &       list_cmp(arc1(iarc)%occ_cnt,arc2(iarc)%occ_cnt,ngastp*2)
      end do

      return
      end
