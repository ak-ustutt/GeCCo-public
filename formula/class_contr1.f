*----------------------------------------------------------------------*
      subroutine class_contr1(lulog,contr,op_info,termnr)
*----------------------------------------------------------------------*
*     classify terms:
*     1 : only one-electron operators except hamiltonian, no R12 ops
*     2 : only zero order two-electron operators, no R12 ops
*     3 : higher order two-electron ops, no R12 ops
*     4 : with R12 ops, only zero order two-electron ops
*     5 : with R12 ops, higher order two-electron ops
*     
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog, termnr
      type(contraction), intent(in) ::
     &     contr
      type(operator_info), intent(in) ::
     &     op_info

      type(operator_array), pointer ::
     &     ops(:)

      character(2) ::
     &     cdag
      character(len_opname) ::
     &     opstr
      character(1) ::
     &     charord
      integer, parameter ::
     &     r12ops = 13, singler12 = 3
      integer ::
     &     idx, idxph(2),idxr12(r12ops), ii, class, opidx(2)
      logical ::
     &     r12, twoel, twoelcur, highord, top(2), lop(2), r12nor(2),
     &     r12dag(2)
      integer, external ::
     &     idx_oplist2

      idxr12(1) = idx_oplist2('V',op_info)
      idxr12(2) = idx_oplist2('C-INT',op_info)
      idxr12(3) = idx_oplist2('R12-INT',op_info)
      idxr12(4) = idx_oplist2('X',op_info)
      idxr12(5) = idx_oplist2('Xh',op_info)
      idxr12(6) = idx_oplist2('B',op_info)
      idxr12(7) = idx_oplist2('Bh',op_info)
      idxr12(8) = idx_oplist2('BVX',op_info)
      idxr12(9) = idx_oplist2('BVY',op_info)
      idxr12(10) = idx_oplist2('BVZ',op_info)
      idxr12(11) = idx_oplist2('CVX',op_info)
      idxr12(12) = idx_oplist2('CVY',op_info)
      idxr12(13) = idx_oplist2('CVZ',op_info)

      ops => op_info%op_arr

      r12 = .false.
      ! contains R12 operators?
      do idx = 1, contr%nvtx
        do ii = 1,r12ops
          r12 = r12.or.contr%vertex(idx)%idx_op.eq.idxr12(ii)
        end do
      end do

      twoel = .false.
      highord = .false.
      ! contains conventional 2-El. amplitudes or multipliers of higher order?
      do idx = 1, contr%nvtx
        if (ops(contr%vertex(idx)%idx_op)%op%species.eq.1.or.
     &      ops(contr%vertex(idx)%idx_op)%op%species.eq.2) then
          twoelcur = (sum(ops(contr%vertex(idx)%idx_op)%op%
     &       ihpvca_occ(1:ngastp,1,contr%vertex(idx)%iblk_op)).gt.1)
          twoel = twoel.or.twoelcur
          if (twoelcur) highord = highord.or.(ops(contr%vertex(idx)%
     &         idx_op)%op%order.gt.0)
        end if
      end do

      if (r12.and..not.highord) then
        ! any composite r12-amplitude of higher order?
        do idx = 1, contr%narc
          opidx(1:2) = contr%vertex(contr%arc(idx)%link(1:2))%idx_op
          idxph(1:2) = 1
          if (contr%vertex(contr%arc(idx)%link(1))%dagger) idxph(1) = 2
          if (contr%vertex(contr%arc(idx)%link(2))%dagger) idxph(2) = 2
          r12nor(1:2) = .false.
          r12dag(1:2) = .false.
          do ii = 1,singler12
            r12nor(1) = r12nor(1).or.(opidx(1).eq.idxr12(ii).and.
     &                                idxph(1).eq.1)
            r12nor(2) = r12nor(2).or.(opidx(2).eq.idxr12(ii).and.
     &                                idxph(2).eq.1)
            r12dag(1) = r12dag(1).or.(opidx(1).eq.idxr12(ii).and.
     &                                idxph(1).eq.2)
            r12dag(2) = r12dag(2).or.(opidx(2).eq.idxr12(ii).and.
     &                                idxph(2).eq.2)
          end do
          do ii = singler12+1,r12ops
            r12nor(1) = r12nor(1).or.(opidx(1).eq.idxr12(ii))
            r12nor(2) = r12nor(2).or.(opidx(2).eq.idxr12(ii))
            r12dag(1) = r12dag(1).or.(opidx(1).eq.idxr12(ii))
            r12dag(2) = r12dag(2).or.(opidx(2).eq.idxr12(ii))
          end do

          top(1) = (ops(opidx(1))%op%species.eq.1.and.
     &              ops(opidx(1))%op%order.gt.0)
          top(2) = (ops(opidx(2))%op%species.eq.1.and.
     &              ops(opidx(2))%op%order.gt.0)
          lop(1) = (ops(opidx(1))%op%species.eq.2.and.
     &              ops(opidx(1))%op%order.gt.0)
          lop(2) = (ops(opidx(2))%op%species.eq.2.and.
     &              ops(opidx(2))%op%order.gt.0)
          highord = highord.or.
     &              (lop(1).and.r12dag(2).and.
     &               contr%arc(idx)%occ_cnt(2,2).gt.0).or.
     &              (lop(2).and.r12dag(1).and.
     &               contr%arc(idx)%occ_cnt(2,1).gt.0).or.
     &              (top(2).and.r12nor(1).and.
     &               contr%arc(idx)%occ_cnt(2,2).gt.0).or.
     &              (top(1).and.r12nor(2).and.
     &               contr%arc(idx)%occ_cnt(2,1).gt.0)
        end do
      end if

      ! define term class
      if (r12.and.highord) then
        class = 5
      else if (r12) then
        class = 4
      else if (twoel.and.highord) then
        class = 3
      else if (twoel) then
        class = 2
      else
        class = 1
      end if

      ! print class
      write(lulog,'(x,i8,x,i8)') termnr,class

      return
      end
