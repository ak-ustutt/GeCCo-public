

*----------------------------------------------------------------------*
      subroutine prt_contr(lulog,contr,ops)
*----------------------------------------------------------------------*
*     write info on contraction onto unit lulog
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(contraction), intent(in) ::
     &     contr
      type(operator), intent(in) ::
     &     ops(*)
      
      integer ::
     &     idx, idxph

      call quit(1,'prt_contr','call to obsolete routine')
      write(lulog,*) '+++ contraction info +++'
      if (contr%idx_res.gt.0) then
        write(lulog,*) ' name (index) and block of result: ',
     &     trim(ops(contr%idx_res)%name),
     &       '(',contr%idx_res,')', contr%iblk_res
      else
        write(lulog,*) ' index and block of result: ',
     &     contr%idx_res, contr%iblk_res
      end if
      write(lulog,*) ' factor: ',contr%fac
      write(lulog,*) ' number of vertices/arcs: ',
     &     contr%nvtx,contr%narc
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.0) then
          write(lulog,'(x,a)') ' v   0'
          cycle
        end if
        idxph = 1
        if (ops(contr%vertex(idx)%idx_op)%dagger) idxph=2
        write(lulog,'(x,a,x,a,i4,2x,4i3)')
     &       ' v  ',ops(contr%vertex(idx)%idx_op)%name(1:4),
     &       contr%vertex(idx)%iblk_op,
     &       ops(contr%vertex(idx)%idx_op)%
     &       ihpvca_occ(1:ngastp,idxph,contr%vertex(idx)%iblk_op)
        idxph = 2
        if (ops(contr%vertex(idx)%idx_op)%dagger) idxph=1
        write(lulog,'(x,a,11x,4i3)')
     &       '    ',ops(contr%vertex(idx)%idx_op)%
     &       ihpvca_occ(1:ngastp,idxph,contr%vertex(idx)%iblk_op)
      end do
      do idx = 1, contr%narc
        write(lulog,'(x,a,x,2i4,2x,4i3)')
     &       ' c  ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2),
     &       contr%arc(idx)%occ_cnt(1:ngastp,1)
        write(lulog,'(x,a,11x,4i3)')
     &       '    ',contr%arc(idx)%occ_cnt(1:ngastp,2)
      end do
      do idx = 1, contr%nfac
        write(lulog,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &       contr%inffac(1:4,idx)
      end do

      return
      end
