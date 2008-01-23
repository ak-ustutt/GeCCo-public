*----------------------------------------------------------------------*
      subroutine prt_contr2(luout,contr,op_info)
*----------------------------------------------------------------------*
*     write info on contraction onto unit luout
*     new version with op_info
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     luout
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator_array), pointer ::
     &     ops(:)
      
      integer ::
     &     idx, idxph

      ops => op_info%op_arr
      write(luout,*) '+++ contraction info +++'
      if (contr%idx_res.gt.0) then
        write(luout,*) ' name (index) and block of result: ',
     &     trim(ops(contr%idx_res)%op%name),
     &       '(',contr%idx_res,')', contr%iblk_res
      else
        write(luout,*) ' index and block of result: ',
     &     contr%idx_res, contr%iblk_res
      end if
      write(luout,*) ' factor: ',contr%fac
      write(luout,'(x,a,3i5)')
     &     ' number of prim.vertices/sup.vertices/arcs: ',
     &     contr%nvtx,contr%nsupvtx,contr%narc
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.0) then
          write(luout,'(2x,"v",i2.2,x,"0")') contr%svertex(idx)
          cycle
        end if
        idxph = 1
        if (ops(contr%vertex(idx)%idx_op)%op%dagger) idxph=2
        write(luout,'(x,"v",i2.2,x,a,i4,2x,4i3)')
     &       contr%svertex(idx),
     &       ops(contr%vertex(idx)%idx_op)%op%name(1:8),
     &       contr%vertex(idx)%iblk_op,
     &       ops(contr%vertex(idx)%idx_op)%op%
     &       ihpvca_occ(1:ngastp,idxph,contr%vertex(idx)%iblk_op)
        idxph = 2
        if (ops(contr%vertex(idx)%idx_op)%op%dagger) idxph=1
        write(luout,'(x,a,14x,4i3)')
     &       '    ',ops(contr%vertex(idx)%idx_op)%op%
     &       ihpvca_occ(1:ngastp,idxph,contr%vertex(idx)%iblk_op)
      end do
      do idx = 1, contr%narc
        if (contr%arc(idx)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
          write(luout,'(x,a,x,2i4)')
     &       ' cp ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2)
        else
          write(luout,'(x,a,x,2i4,5x,4i3)')
     &       ' c  ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2),
     &       contr%arc(idx)%occ_cnt(1:ngastp,1)
          write(luout,'(x,a,14x,4i3)')
     &       '    ',contr%arc(idx)%occ_cnt(1:ngastp,2)
        end if
      end do
      do idx = 1, contr%nxarc
        if (contr%xarc(idx)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
          write(luout,'(x,a,x,2i4)')
     &       ' xp ',contr%xarc(idx)%link(1),
     &       contr%xarc(idx)%link(2)
        else
          write(luout,'(x,a,x,2i4,5x,4i3)')
     &       ' x  ',contr%xarc(idx)%link(1),
     &       contr%xarc(idx)%link(2),
     &       contr%xarc(idx)%occ_cnt(1:ngastp,1)
          write(luout,'(x,a,14x,4i3)')
     &       '    ',contr%xarc(idx)%occ_cnt(1:ngastp,2)
        end if
      end do
      do idx = 1, contr%nfac
        write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &       contr%inffac(1:4,idx)
      end do

      return
      end
