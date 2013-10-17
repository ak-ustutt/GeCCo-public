*----------------------------------------------------------------------*
      subroutine get_unconnected4vertex2(iocc,ivtx,
&                                        contr,occ_vtx,op_info)
*----------------------------------------------------------------------*
*     return (in occupation form) the number of unconnected indices 
*     of vertex #ivtx in contraction contr
*     version with predefined occ_vtx (old version could not handle
*     vertices that come from result operators)
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     iocc(ngastp,2)
      integer, intent(in) ::
     &     ivtx
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,contr%nvtx)

      logical ::
     &     dag
      integer ::
     &     idx_op, iblk_op, iarc
      type(operator_array), pointer ::
     &     op_arr(:)
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)

      if (ntest.ge.100) then
        call prt_contr2(luout,contr,op_info)
      end if

      ! get vertex occupation
      iocc(:,:) = occ_vtx(:,:,ivtx)
      
      if (ntest.ge.100) then
        write(luout,*) 'initial'
        call wrt_occ(luout,iocc)
      end if

      ! loop over all arcs and remove contractions
      arc => contr%arc
      do iarc = 1, contr%narc
        ! for proto-contractions: ignore certain arcs
        if (arc(iarc)%occ_cnt(1,1).lt.0) cycle
        if (arc(iarc)%link(1).eq.ivtx) then
          iocc = iocc - arc(iarc)%occ_cnt
        else if (arc(iarc)%link(2).eq.ivtx) then
          iocc = iocc - iocc_dagger(arc(iarc)%occ_cnt)
        end if
        if (ntest.ge.100) then
          write(luout,*) 'after arc ',iarc
          call wrt_occ(luout,iocc)
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'final'
        call wrt_occ(luout,iocc)
      end if

      return
      end
