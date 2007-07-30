*----------------------------------------------------------------------*
      subroutine prt_contr3(luout,contr,occ_vtx)
*----------------------------------------------------------------------*
*     write info on contraction onto unit luout
*     alternative version if no op_info available
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     luout
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,*)
      
      integer ::
     &     idx, idxph

      write(luout,*) '+++ contraction info +++'
      if (contr%idx_res.gt.0) then
        write(luout,*) ' index and block of result: ',
     &     contr%idx_res, contr%iblk_res
      else
        write(luout,*) ' index and block of result: ',
     &     contr%idx_res, contr%iblk_res
      end if
      write(luout,*) ' factor: ',contr%fac
      write(luout,*) ' number of vertices/arcs: ',
     &     contr%nvtx,contr%narc
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.0) then
          write(luout,'(x,a)') ' v   0'
          cycle
        end if
        write(luout,'(x,a,x,i5,x,i4,2x,4i3)')
     &       ' v  ',contr%vertex(idx)%idx_op,
     &       contr%vertex(idx)%iblk_op,
     &       occ_vtx(1:ngastp,1,idx+1)
        write(luout,'(x,a,13x,4i3)')
     &       '    ',
     &       occ_vtx(1:ngastp,2,idx+1)
      end do
      do idx = 1, contr%narc
        if (contr%arc(idx)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
          write(luout,'(x,a,x,2i4)')
     &       ' cp ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2)
        else
          write(luout,'(x,a,x,2i4,2x,4i3)')
     &       ' c  ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2),
     &       contr%arc(idx)%occ_cnt(1:ngastp,1)
          write(luout,'(x,a,11x,4i3)')
     &       '    ',contr%arc(idx)%occ_cnt(1:ngastp,2)
        end if
      end do
      do idx = 1, contr%nfac
        write(luout,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &       contr%inffac(1:4,idx)
      end do

      return
      end
