*----------------------------------------------------------------------*
      subroutine check_xarcs(contr,op_info)
*----------------------------------------------------------------------*
*     check consistency of xarc info
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     idxres, nvtx, narc, nxarc, njoined_res, iarc, ivtx1, ivtx2

      type(cntr_arc), pointer ::
     &     arc(:), xarc(:)
      type(operator_array), pointer ::
     &     op_arr(:)
      integer, pointer ::
     &     occ_vtx(:,:,:), occ_vtx2(:,:,:)

c      logical, external ::
c     &     iocc_equal_n

      op_arr => op_info%op_arr

      nvtx   = contr%nvtx
      narc   = contr%narc
      nxarc   = contr%nxarc
      idxres = contr%idx_res

      arc => contr%arc
      xarc => contr%xarc

      njoined_res = op_arr(idxres)%op%njoined

      allocate(occ_vtx(ngastp,2,nvtx+njoined_res),
     &        occ_vtx2(ngastp,2,nvtx+njoined_res))

      call occvtx4contr(0,occ_vtx,contr,op_info)

      write(luout,*) 'checking xarcs for:'
      call prt_contr3(luout,contr,occ_vtx(1,1,njoined_res+1))

      ! check consistency of xarcs with result occupation
      occ_vtx2 = 0
      do iarc = 1, nxarc
        ivtx1 = xarc(iarc)%link(1)+njoined_res
        ivtx2 = xarc(iarc)%link(2)
        occ_vtx2(1:ngastp,1:2,ivtx2) =
     &       occ_vtx2(1:ngastp,1:2,ivtx2) +
     &       xarc(iarc)%occ_cnt
        occ_vtx2(1:ngastp,1:2,ivtx1) =
     &       occ_vtx2(1:ngastp,1:2,ivtx1) +
     &       xarc(iarc)%occ_cnt
      end do

      write(luout,*) 'result: 1) from op_info, 2) from xarcs'
      call wrt_occ_n(luout,occ_vtx ,njoined_res)
      call wrt_occ_n(luout,occ_vtx2,njoined_res)

      if (.not.iocc_equal_n(occ_vtx,.false.,
     &                     occ_vtx2,.false.,njoined_res))
     &     call quit(1,'check_xarcs','result is not consistent!')

      ! check consistency of xarcs with vertex occupations
      do iarc = 1, narc
        ivtx1 = arc(iarc)%link(1)+njoined_res
        ivtx2 = arc(iarc)%link(2)+njoined_res
        occ_vtx(1:ngastp,1:2,ivtx1) =
     &       occ_vtx(1:ngastp,1:2,ivtx1) -
     &       arc(iarc)%occ_cnt
        occ_vtx(1:ngastp,1:2,ivtx2) =
     &       occ_vtx(1:ngastp,1:2,ivtx2) -
     &       iocc_dagger(arc(iarc)%occ_cnt)
      end do

      write(luout,*)
     &     'vertex occupations after stripping contracted indices'
      write(luout,*) '1) from contraction,  2) from xarcs'
      call wrt_occ_n(luout,occ_vtx (1,1,njoined_res+1),nvtx)
      call wrt_occ_n(luout,occ_vtx2(1,1,njoined_res+1),nvtx)

      if (.not.iocc_equal_n(
     &    occ_vtx (1:ngastp,1:2,njoined_res+1:njoined_res+nvtx),.false.,
     &    occ_vtx2(1:ngastp,1:2,njoined_res+1:njoined_res+nvtx),.false.,
     &                      nvtx))
     &     call quit(1,'check_xarcs','vertices are not consistent!')
        
      deallocate(occ_vtx, occ_vtx2)

      return
      end
