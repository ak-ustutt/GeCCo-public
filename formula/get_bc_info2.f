*----------------------------------------------------------------------*
      subroutine get_bc_info2(idx_op,iblk_op,
     &     iocc_ext,iocc_cnt,
     &     iocc_op,iocc_int,
     &     irestr_op,irestr_int,
     &     mst_op,mst_int,
     &     igamt_op,igamt_int,
     &     contr,occ_vtx,irestr_vtx,info_vtx,iarc,
     &     irestr_res,ihpvgas,ngas)
*----------------------------------------------------------------------*
*     set info for binary contraction
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'ifc_operators.h'
      include 'def_contraction.h'
      include 'multd2h.h'

      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     ngas, iarc
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx+1),
     &     irestr_vtx(2,ngas,2,2,contr%nvtx+1),
     &     info_vtx(2,contr%nvtx+1),
     &     irestr_res(2,ngas,2,2),
     &     ihpvgas(ngas)
      integer, intent(out) ::
     &     idx_op(2), iblk_op(2),
     &     iocc_ext(ngastp,2,2), iocc_cnt(ngastp,2),
     &     iocc_op(ngastp,2,2),  iocc_int(ngastp,2),
     &     irestr_op(2,ngas,2,2,2), irestr_int(2,ngas,2,2),
     &     mst_op(2), mst_int, igamt_op(2), igamt_int

      integer ::
     &     ivtx1, ivtx2

      ! set up operator 1 and 2
      ivtx1 = contr%arc(iarc)%link(1)
      ivtx2 = contr%arc(iarc)%link(2)

      if (ivtx1.le.contr%nvtx) then
        idx_op(1) = contr%vertex(ivtx1)%idx_op
        iblk_op(1) = contr%vertex(ivtx1)%iblk_op
      end if

      if (ivtx2.le.contr%nvtx) then
        idx_op(2) = contr%vertex(ivtx2)%idx_op
        iblk_op(2) = contr%vertex(ivtx2)%iblk_op
      end if

      iocc_op(1:ngastp,1:2,1) = occ_vtx(1:ngastp,1:2,ivtx1+1)
      iocc_op(1:ngastp,1:2,2) = occ_vtx(1:ngastp,1:2,ivtx2+1)
      irestr_op(1:2,1:ngas,1:2,1:2,1) =
     &     irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx1+1)
      irestr_op(1:2,1:ngas,1:2,1:2,2) =
     &     irestr_vtx(1:2,1:ngas,1:2,1:2,ivtx2+1)
      mst_op(1) = info_vtx(1,ivtx1+1)
      mst_op(2) = info_vtx(1,ivtx2+1)
      igamt_op(1) = info_vtx(2,ivtx1+1)
      igamt_op(2) = info_vtx(2,ivtx2+1)

      ! set external and contraction indices
      iocc_cnt = contr%arc(iarc)%occ_cnt
      
      iocc_ext(1:ngastp,1:2,1) = iocc_op(1:ngastp,1:2,1)
     &                         - iocc_cnt
      iocc_ext(1:ngastp,1:2,2) = iocc_op(1:ngastp,1:2,2)
     &                         - iocc_dagger(iocc_cnt)

      ! set intermediate/result
      if (contr%narc.eq.1) then
        iocc_int(1:ngastp,1:2) = occ_vtx(1:ngastp,1:2,1)
        irestr_int(1:2,1:ngas,1:2,1:2) =
     &       irestr_vtx(1:2,1:ngas,1:2,1:2,1)
        mst_int = info_vtx(1,1)
        igamt_int = info_vtx(2,1)
      else
        iocc_int = iocc_ext(1:ngastp,1:2,1)
     &           + iocc_ext(1:ngastp,1:2,2)
        ! yes, the preliminary code
        call fit_restr(irestr_int,iocc_int,
     &     irestr_res,ihpvgas,ngas)
        mst_int = mst_op(1) + mst_op(2)
        igamt_int = multd2h(igamt_op(1),igamt_op(2)) 
      end if

      return
      end
