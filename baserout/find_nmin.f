*----------------------------------------------------------------------*
      subroutine find_nmin(xlist,idxlist,nlist,vec,ndim,ioff,init)
*----------------------------------------------------------------------*
*     scan vec(1:ndim) and update the list of nlist lowest values
*     on xlist; idxlist contains the index of the respective values
*     (including the offset given by ioff); init initialises xlist
*
*     by subsequent calls we can scan over blocks vector vec
*----------------------------------------------------------------------*

      implicit none

      logical, intent(in) ::
     &     init
      integer, intent(in) ::
     &     ndim, nlist, ioff
      real(8), intent(in) ::
     &     vec(ndim)
      integer, intent(inout) ::
     &     idxlist(nlist)
      real(8), intent(inout) ::
     &     xlist(nlist)

      integer ::
     &     idx, jdx, kdx, mndx

      mndx = 0
      if (init) then
        mndx = min(nlist,ndim)
        do idx = 1, mndx
          xlist(idx) = vec(idx)
          idxlist(idx) = idx
        end do
        do idx = mndx+1, nlist
          xlist(idx) = huge(xlist(idx))
          idxlist(idx) = idx
        end do
        call idxsortx(xlist,idxlist,nlist,+1)
      end if

      ! scan vec
      do idx = mndx+1, ndim
        ! scan list
        do jdx = 1, nlist
          ! vec-element < list-element?
          if (vec(idx).lt.xlist(jdx)) then
            ! insert into list
            do kdx = nlist, jdx+1, -1
              xlist(kdx) = xlist(kdx-1)
              idxlist(kdx) = idxlist(kdx-1)
            end do
            xlist(jdx) = vec(idx)
            idxlist(jdx) = idx+ioff
            exit
          end if
        end do
      end do

      return

      end
