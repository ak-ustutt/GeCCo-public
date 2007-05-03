*----------------------------------------------------------------------*
      subroutine reduce_graph(ignew,narc_new,iarceq,
     &     icnt,ivtxnew,ig,narc,narc_total)
*----------------------------------------------------------------------*
*     generate graph after contraction icnt giving new vertex ivtxnew
*     important: keep sequence of operators !
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     narc, narc_total, icnt, ivtxnew, ig(3,narc)

      integer, intent(out) ::
     &     ignew(3,*), narc_new, iarceq(2,narc_total)

      logical ::
     &     new
      integer ::
     &     idx, jdx, kdx, ivtx1, ivtx2, idxls, iinew

      if (ntest.ge.100) then
        write(luout,*) '======================'
        write(luout,*) ' reduce_graph at work'
        write(luout,*) '======================'
        write(luout,*) ' icnt, ivtxnew: ', icnt, ivtxnew
        write(luout,*) ' graph described by: '
        do idx = 1, narc
          write(luout,*) '  ',ig(1,idx),':',ig(2,idx),'-',ig(3,idx)
        end do
      end if

      ! get indices of contracted vertices
      do idx = 1, narc
        if (ig(1,idx).eq.icnt) then
          ivtx1 = ig(2,idx)
          ivtx2 = ig(3,idx)
          exit
        end if
      end do

      ! loop over arcs of old graph
      narc_new = 0
      do idx = 1, narc
        ! is this the contracted arc?
        if (ig(1,idx).eq.icnt) cycle ! skip it
        ! else we copy
        jdx = narc_new+1
        ignew(1:3,jdx) = ig(1:3,idx)
        iinew = -1
        if (ignew(2,jdx).eq.ivtx1.or.ignew(2,jdx).eq.ivtx2) then
c          ignew(2,jdx) = ignew(3,jdx)
c          ignew(3,jdx) = ivtxnew  ! set new index always on 2nd position
          ignew(3,jdx) = ignew(3,jdx)
          ignew(2,jdx) = ivtxnew
          iinew = 3
        else if (ignew(3,jdx).eq.ivtx1.or.ignew(3,jdx).eq.ivtx2) then
          ignew(3,jdx) = ivtxnew
          iinew = 2
        end if
        ! goes arc to one of the contracted vertices?
        new = .true.
        ! fix: treat ignew(2,jdx) case
        if (iinew.gt.0) then
          ! compare with previous arcs on new graph
          do kdx = 1, narc_new
            if ((ignew(2,kdx).eq.ignew(iinew,jdx).and.
     &           ignew(3,kdx).eq.ivtxnew) .or.
     &          (ignew(3,kdx).eq.ignew(iinew,jdx).and.
     &           ignew(2,kdx).eq.ivtxnew)) then
              ! not new -> fused connection
              ! ignew(1,jdx) is now equivalent to ignew(1,kdx)
              ! replace in equivalence array all references to that node:
              do idxls = 1, narc_total
                if (iarceq(2,idxls).eq.ignew(1,jdx))
     &               iarceq(2,idxls) = ignew(1,kdx)
              end do
              new = .false.
            end if
          end do
        end if
        ! indeed new: increment counter
        if (new) narc_new = narc_new+1

      end do

      if (ntest.ge.100) then
        write(luout,*) ' new graph described by: '
        do idx = 1, narc_new
          write(luout,*) '  ',ignew(1,idx),':',
     &         ignew(2,idx),'-',ignew(3,idx)
        end do
        write(luout,*) ' equivalent arcs:'
        write(luout,'(x,a,20i4)') ' node   > ',iarceq(1,1:narc_total)
        write(luout,'(x,a,20i4)') ' equiv. > ',iarceq(2,1:narc_total)
      end if

      return
      end
