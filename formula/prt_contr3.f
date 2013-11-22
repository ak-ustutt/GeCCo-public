*----------------------------------------------------------------------*
      subroutine prt_contr3(lulog,contr,occ_vtx_in)
*----------------------------------------------------------------------*
*     write info on contraction onto unit lulog
*     alternative version if no op_info available
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lulog
      type(contraction), intent(in) ::
     &     contr
      integer, intent(in), target ::
     &     occ_vtx_in(ngastp,2,contr%nvtx)
      
      character(1) ::
     &     cdag
      integer ::
     &     idx, idxph

      integer, pointer ::
     &     occ_vtx(:,:,:)

      if (occ_vtx_in(1,1,1).lt.0) then
        allocate(occ_vtx(ngastp,2,contr%nvtx))
        call occvtx_from_arcs(1,occ_vtx,contr,0)
      else
        occ_vtx => occ_vtx_in
      end if

      write(lulog,*) '+++ contraction info +++'
      cdag = '  '
      if (contr%dagger) cdag = '^+'
      if (contr%idx_res.gt.0) then
        write(lulog,*) ' index and block of result: ',
     &     contr%idx_res, cdag, contr%iblk_res
      else
        write(lulog,*) ' index and block of result: ',
     &     contr%idx_res, cdag, contr%iblk_res
      end if
      write(lulog,*) ' factor: ',contr%fac
      write(lulog,'(x,a,3i5)')
     &     ' number of prim.vertices/sup.vertices/arcs: ',
     &     contr%nvtx,contr%nsupvtx,contr%narc
      do idx = 1, contr%nvtx
        if (contr%vertex(idx)%idx_op.eq.0) then
          write(lulog,'(2x,"v",i2.2,x,"0")') contr%svertex(idx)
          cycle
        end if
        cdag = ' '
        if (contr%vertex(idx)%dagger) cdag = '+'
        write(lulog,'(x,"v",i2.2,x,i5,a,x,i4,2x,4i3)')
     &       contr%svertex(idx),contr%vertex(idx)%idx_op,cdag,
     &       contr%vertex(idx)%iblk_op,
     &       occ_vtx(1:ngastp,1,idx)
        write(lulog,'(x,a,13x,4i3)')
     &       '    ',
     &       occ_vtx(1:ngastp,2,idx)
      end do
      do idx = 1, contr%narc
        if (contr%arc(idx)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
          write(lulog,'(x,a,x,2i4)')
     &       ' cp ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2)
        else
          write(lulog,'(x,a,x,2i4,2x,4i3)')
     &       ' c  ',contr%arc(idx)%link(1),
     &       contr%arc(idx)%link(2),
     &       contr%arc(idx)%occ_cnt(1:ngastp,1)
          write(lulog,'(x,a,11x,4i3)')
     &       '    ',contr%arc(idx)%occ_cnt(1:ngastp,2)
        end if
      end do
      do idx = 1, contr%nxarc
        if (contr%xarc(idx)%occ_cnt(1,1).lt.0) then
          ! prototype connection:
          write(lulog,'(x,a,x,2i4)')
     &       ' xp ',contr%xarc(idx)%link(1),
     &       contr%xarc(idx)%link(2)
        else
          write(lulog,'(x,a,x,2i4,5x,4i3)')
     &       ' x  ',contr%xarc(idx)%link(1),
     &       contr%xarc(idx)%link(2),
     &       contr%xarc(idx)%occ_cnt(1:ngastp,1)
          write(lulog,'(x,a,14x,4i3)')
     &       '    ',contr%xarc(idx)%occ_cnt(1:ngastp,2)
        end if
      end do
      do idx = 1, contr%nfac
        write(lulog,'(x,i5,"*",i5,"->",i5,"(",i5,")")')
     &       contr%inffac(1:4,idx)
      end do

      if (occ_vtx_in(1,1,1).lt.0) then
        deallocate(occ_vtx)
      end if

      return
      end
