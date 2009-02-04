*----------------------------------------------------------------------*
      subroutine weight_ssg(iyssg,iwssg,nelmax,nspc,mnmxspc)
*----------------------------------------------------------------------*
*     generate weight arrays for subspace distribution graphs
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(out) ::
     &     iyssg(1:nelmax,nspc), iwssg(0:nelmax,nspc)
!                ^                     ^  !!!
      integer, intent(in) ::
     &     nelmax, nspc, mnmxspc(2,nspc)

      integer ::
     &     ispc, iel

      if (ntest.gt.0) then
        write(luout,*) '------------------'
        write(luout,*) 'Entered weight_ssg'
        write(luout,*) '------------------'
      end if
      if (ntest.ge.100) then
        write(luout,*) 'restriction array:'
        do ispc = 1, nspc
          write(luout,'(x,i2,5x,2i4)') 
     &          ispc,mnmxspc(1,ispc),mnmxspc(2,ispc)
        end do
      end if

      iwssg(0:nelmax,1:nspc) = 0
      iyssg(1:nelmax,1:nspc) = 0

      ! Form vertex weights first.
      do ispc = 1, nspc
        do iel = 0, nelmax

          if (ispc.eq.1.and.iel.le.mnmxspc(2,1)) then
            iwssg(iel,ispc) = 1
c dbg fix by mh
c dbg original line          else if (ispc.gt.1.and.iel.ge.mnmxspc(1,ispc-1)
          else if (ispc.gt.1) then
          if (iel.ge.mnmxspc(1,ispc-1)
c dbg original
     &                      .and.iel.le.mnmxspc(2,ispc)) then
            iwssg(iel,ispc) = iwssg(iel,ispc-1)
            if (iel.gt.0) iwssg(iel,ispc)
     &           = iwssg(iel,ispc)+iwssg(iel-1,ispc)
c dbg resume fix
          end if
c dbg end fix
          end if

        end do
      end do

      ! make arc weight array ...
      do ispc = 2, nspc
        do iel = 1, nelmax
          if (iwssg(iel,ispc).ne.0) then
            iyssg(iel,ispc) = iwssg(iel,ispc-1)
          end if
        end do
      end do

      if (ntest.ge.100) then
        write(luout,*) 'Generated vertex weights:'
        do ispc = 1, nspc
          write(luout,'(2x,i2,5x,10i6)') ispc,iwssg(0:nelmax,ispc)
        end do
        write(luout,*) 'Generated arc weights:'
        do ispc = 1, nspc
          write(luout,'(2x,i2,8x,10i6)') ispc,iyssg(1:nelmax,ispc)
        end do
      end if

      return
      end
