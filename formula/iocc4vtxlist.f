*----------------------------------------------------------------------*
      subroutine iocc4vtxlist(iocc,icnt,contr,vtxlist,nlist,op_info)
*----------------------------------------------------------------------*
*     carry out all contraction described in contr, that involve
*     the vertices on list vtxlist (in the given sequence)
*     the resulting intermediate is represented by iocc
*     contraction arcs to vertices not include in list are gathered
*     on icnt
*
*     this version allows only contraction of consecutive operators
*     (i.e. no density-type "inner" open lines must remain)
*
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'
      include 'ifc_operators.h'
      include 'ifc_baserout.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     iocc(ngastp,2), icnt(ngastp,2)
      integer, intent(in) ::
     &     nlist, vtxlist(nlist)
      type(contraction) ::
     &     contr
      type(operator_info) ::
     &     op_info

      integer ::
     &     jocc(ngastp,2)
      integer ::
     &     nvtx, narc, idx, idx_op, iblk_op
      logical ::
     &     dag
      type(cntr_vtx), pointer ::
     &     vertex(:)
      type(cntr_arc), pointer ::
     &     arc(:)
      type(operator_array), pointer ::
     &     op_arr(:)


      op_arr => op_info%op_arr
      if (ntest.ge.100) then
        write(lulog,*) '----------------------'
        write(lulog,*) ' iocc4vtxlist at work'
        write(lulog,*) '----------------------'
        call prt_contr2(lulog,contr,op_arr)
        write(lulog,*) 'vertex list: ',vtxlist(1:nlist)
      end if

      ! initialize
      iocc(1:ngastp,1:2) = 0
      icnt(1:ngastp,1:2) = 0

      nvtx = contr%nvtx
      vertex => contr%vertex
      ! add all occupations of involved operators
      do idx = 1, nlist
        idx_op = vertex(vtxlist(idx))%idx_op
        iblk_op = vertex(vtxlist(idx))%iblk_op
        ! get occupation of operator block and add
        jocc = op_arr(idx_op)%op%ihpvca_occ(1:ngastp,1:2,iblk_op)
        dag  = op_arr(idx_op)%op%dagger
        iocc = iocc_add(1,iocc,.false.,1,jocc,dag)
      end do

      narc = contr%narc
      arc => contr%arc
      ! remove all contraction parts
      do idx = 1, narc
        ! both vertices included in list?
        if (idxlist(arc(idx)%link(1),vtxlist,nlist,1).gt.0.and.
     &      idxlist(arc(idx)%link(2),vtxlist,nlist,1).gt.0) then
          ! yes: modify result of iocc by
          ! removing [C]+[C^+]
          jocc = arc(idx)%occ_cnt
          jocc = iocc_add(1,jocc,.false.,1,jocc,.true.)
          iocc = iocc_add(1,iocc,.false.,1,jocc,.false.)
        else
          ! no: collect on icnt
          jocc = arc(idx)%occ_cnt
          icnt = iocc_add(1,icnt,.false.,1,jocc,.false.)
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) 'resulting occupation'
        call wrt_occ(lulog,iocc)
        write(lulog,*) 'there-of already contracted'
        call wrt_occ(lulog,icnt)
      end if

      return
      end
