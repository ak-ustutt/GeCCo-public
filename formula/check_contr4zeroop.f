*----------------------------------------------------------------------*
      logical function check_contr4zeroop(contr,op_info)
*----------------------------------------------------------------------*
*     check, whether the formula contains operator blocks that are
*     associated with zero length matrix elements
*
*     andreas, nov 2012
*----------------------------------------------------------------------*
      implicit none
      
      include 'opdim.h'
      include 'def_contraction.h'
      include 'mdef_operator_info.h'

      integer, parameter ::
     &     mxfound = 10

      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      integer ::
     &     nsupvtx, nj, ivtx, isupvtx, idxmel, idxop, iblkop
      character(len=mxlen_melabel) :: 
     &     me_name
      type(cntr_vtx), pointer ::
     &     vertex(:), pvertex(:)

      integer, external ::
     &     idx_mel_list


      check_contr4zeroop = .true.
      
      ! loop over super operator vertices
      nsupvtx = contr%nsupvtx
      do isupvtx = 1, nsupvtx
        ! we visit only the first vertex of the super operator
        nj   = contr%joined(0,isupvtx)
        ivtx = contr%joined(1,isupvtx)

        ! get operator number and block
        idxop  = contr%vertex(ivtx)%idx_op
        iblkop = (contr%vertex(ivtx)%iblk_op-1)/nj+1

        ! get name of associated ME list
        me_name = op_info%op_arr(idxop)%op%assoc_list

        ! look up info on that list
        idxmel = idx_mel_list(me_name,op_info)

        ! is the length == 0 ??
        if (op_info%mel_arr(idxmel)%mel%len_op_occ(iblkop).eq.0)
     &     return  ! then return

      end do

      check_contr4zeroop = .false.

      return
      end

