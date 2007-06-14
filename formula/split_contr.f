*----------------------------------------------------------------------*
      subroutine split_contr(contr_rem,contr_spl,contr)
*----------------------------------------------------------------------*
*     given a contraction contr and a contraction contr_spl (which
*     must be contained in contr, check with contr_in_contr()!),
*     obtain the remainder on contr_rem
*
*     a empty node (idx 0, block 0) is used to indicate the position
*     of the part that has been "cut out"
*
*  ?? maybe not: ??
*     in contrast to usual contractions, we also store external 
*     contraction arcs for contr_rem, i.e. those which connected
*     contr_rem with contr_spl
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) ::
     &     contr, contr_spl
      type(contraction), intent(out) ::
     &     contr_rem

      integer ::
     &     nvtx_t, nvtx_spl, nvtx_rem, idx_spl, idx, idx_rem,
     &     narc_t, narc_spl, narc_rem, idx1, idx2, idx_t,
     &     idxop, iblkop

      integer, pointer ::
     &     ivtx4ivtx_t(:,:)

      ! prefactor
      if (abs(contr_spl%fac).lt.1d-100)
     &     call quit(1,'split_contr','division by zero encountered')
      contr_rem%fac = contr%fac / contr_spl%fac
      
      ! ensure sufficient memory for vertices, arcs
      nvtx_t = contr%nvtx
      narc_t = contr%narc

      ! vertices
      nvtx_t = contr%nvtx
      nvtx_spl = contr_spl%nvtx
      nvtx_rem = nvtx_t - nvtx_spl + 1 ! one additional for
                                   ! position of removed part
      
      ! arcs
      narc_t = contr%narc
      narc_spl = contr_spl%narc
      narc_rem = narc_t - narc_spl

      call resize_contr(contr_rem,nvtx_rem,narc_rem,0)
      contr_rem%nvtx = nvtx_rem
      contr_rem%narc = narc_rem
c
c      allocate(contr_rem%vertex(nvtx_rem))

      allocate(ivtx4ivtx_t(2,nvtx_t))

      idx_spl = 1
      idx_rem = 0
      do idx_t = 1, nvtx_t
        idxop  = contr%vertex(idx_t)%idx_op
        iblkop = contr%vertex(idx_t)%iblk_op
        if (idxop.eq.contr_spl%vertex(idx_spl)%idx_op .and.
     &      iblkop.eq.contr_spl%vertex(idx_spl)%iblk_op) then
          if (idx_spl.eq.1) then
            idx_rem = idx_rem+1
            contr_rem%vertex(idx_rem)%idx_op = 0
            contr_rem%vertex(idx_rem)%iblk_op = 0
          end if
          ivtx4ivtx_t(1,idx_t) = 1
          ivtx4ivtx_t(2,idx_t) = idx_spl
          idx_spl = idx_spl+1
        else          
          idx_rem = idx_rem+1
          contr_rem%vertex(idx_rem)%idx_op = idxop
          contr_rem%vertex(idx_rem)%iblk_op = iblkop
          ivtx4ivtx_t(1,idx_t) = 2
          ivtx4ivtx_t(2,idx_t) = idx_rem
        end if
      end do

      ! arcs
      idx_rem = 0
      arc_loop: do idx_t = 1, narc_t
        idx1 = contr%arc(idx_t)%link(1)
        idx2 = contr%arc(idx_t)%link(2)
        ! both links on contr_spl? ignore
        if (ivtx4ivtx_t(1,idx1).eq.1.and.ivtx4ivtx_t(1,idx2).eq.1) cycle
        ! else: store information
        if (ivtx4ivtx_t(1,idx1).eq.2.and.ivtx4ivtx_t(1,idx2).eq.2) then
          idx_rem = idx_rem+1
          contr_rem%arc(idx_rem)%link(1) = ivtx4ivtx_t(2,idx1)
          contr_rem%arc(idx_rem)%link(2) = ivtx4ivtx_t(2,idx2)
          contr_rem%arc(idx_rem)%occ_cnt = contr%arc(idx_t)%occ_cnt
        else 
          if (ivtx4ivtx_t(1,idx1).eq.2) then
            ! look whether already one external line exists:
            do idx = 1, idx_rem
              if (contr_rem%arc(idx)%link(1).eq.ivtx4ivtx_t(2,idx1))then
                contr_rem%arc(idx)%occ_cnt =
     &               contr_rem%arc(idx)%occ_cnt +
     &               contr%arc(idx_t)%occ_cnt
                cycle arc_loop
              end if
            end do
            idx_rem = idx_rem+1
            contr_rem%arc(idx_rem)%link(1) = ivtx4ivtx_t(2,idx1)
            contr_rem%arc(idx_rem)%link(2) = 0 ! external link
            contr_rem%arc(idx_rem)%occ_cnt = contr%arc(idx_t)%occ_cnt
          else ! if (ivtx4ivtx_t(1,idx2).eq.2) then
            do idx = 1, idx_rem
              if (contr_rem%arc(idx)%link(2).eq.ivtx4ivtx_t(2,idx2))then
                contr_rem%arc(idx)%occ_cnt =
     &               contr_rem%arc(idx)%occ_cnt +
     &               contr%arc(idx_t)%occ_cnt
                cycle arc_loop
              end if
            end do
            idx_rem = idx_rem+1
            contr_rem%arc(idx_rem)%link(1) = 0 ! external link
            contr_rem%arc(idx_rem)%link(2) = ivtx4ivtx_t(2,idx2)
            contr_rem%arc(idx_rem)%occ_cnt = contr%arc(idx_t)%occ_cnt
          end if
          
        end if
      end do arc_loop

      contr_rem%narc = idx_rem  ! set to actually counted number of arcs

      deallocate(ivtx4ivtx_t)
      
      return
      end
