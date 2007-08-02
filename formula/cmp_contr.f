*----------------------------------------------------------------------*
      logical function cmp_contr(contr1,contr2,ignore_fac)
*----------------------------------------------------------------------*
*     compare contractions
*     we assume canonical ordering, so the contractions are equivalent
*     if all entries are equivalent
*     for some cases, we might ignore the factor
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr1, contr2
      logical, intent(in) ::
     &     ignore_fac

      integer ::
     &     iarc, ivtx
      type(cntr_vtx), pointer ::
     &     vtx1(:), vtx2(:)
      type(cntr_arc), pointer ::
     &     arc1(:), arc2(:)

      integer, external ::
     &     list_cmp

      cmp_contr = contr1%nvtx.eq.contr2%nvtx .and.
     &            contr1%narc.eq.contr2%narc .and.
     &            contr1%nsupvtx.eq.contr2%nsupvtx
      cmp_contr = cmp_contr.and.
     &            contr1%idx_res.eq.contr2%idx_res .and.
     &            contr1%iblk_res.eq.contr2%iblk_res

      if (.not.cmp_contr) return

      if (.not.ignore_fac)
     &     cmp_contr = abs(contr1%fac-contr2%fac).lt.1d-12

      if (.not.cmp_contr) return
      
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
