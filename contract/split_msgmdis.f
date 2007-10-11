*----------------------------------------------------------------------*
      subroutine split_msgmdis(msdis_1,gamdis_1,
     &                         msdis_2,gamdis_2,
     &                         msdis_r,gamdis_r,
     &                         nel_r,
     &                         map_inf,inv12)
*----------------------------------------------------------------------*
*     input: Ms and IRREP distributions msdis_2, gamdis_2 and
*            msdis_r, gmdis_r (all in condensed form)
*     get the msdis_1, gmdis_1 which (using the map map_inf) merge
*          with msdis_2/gmdis_2 to msdis_r/gmdis_r
*     inv12: if true 1 and 2 are reversed
*----------------------------------------------------------------------*
      implicit none

      include 'multd2h.h'

      logical ::
     &     inv12
      integer, intent(in) ::
     &     nel_r,
     &     msdis_r(nel_r), gamdis_r(nel_r),
     &     msdis_2(*), gamdis_2(*),
     &     map_inf(*)
      integer, intent(out) ::
     &     msdis_1(*), gamdis_1(*)

      logical ::
     &     err
      integer ::
     &     idx, jdx, nel, iel
      integer ::
     &     msdis_scr(nel_r), gamdis_scr(nel_r)

      err = .false.
      if (.not.inv12) then
        msdis_scr = msdis_r
        gamdis_scr = gamdis_r
c dbg
c        print *,'msdis_scr:',msdis_scr
c dbg

        idx = 0
        do jdx = 1, nel_r
          idx = idx+1
          nel = map_inf(idx)
          idx = idx+nel+1
          nel = map_inf(idx)
          do iel = 1, nel
            idx = idx+1
            msdis_scr(jdx) = msdis_scr(jdx) - msdis_2(map_inf(idx))
            gamdis_scr(jdx) = multd2h(gamdis_scr(jdx),
     &           gamdis_2(map_inf(idx)))
          end do
        end do
c dbg
c        print *,'msdis_scr:',msdis_scr
c dbg
        idx = 0
        do jdx = 1, nel_r
          idx = idx+1
          nel = map_inf(idx)
          err = err.or.nel.gt.1
          do iel = 1, nel
            idx = idx+1
            msdis_1(map_inf(idx)) = msdis_scr(jdx)
            gamdis_1(map_inf(idx)) = gamdis_scr(jdx) 
          end do
          idx = idx+1
          nel = map_inf(idx)
          idx = idx+nel
        end do
        if (err)
     &       call quit(1,'split_msgmdis','non-invertible map')
c dbg
c        print *,'msdis_1:',msdis_1(1:2)
c dbg
      else
        msdis_scr = msdis_r
        gamdis_scr = gamdis_r

        idx = 0
        do jdx = 1, nel_r
          idx = idx+1
          nel = map_inf(idx)
          do iel = 1, nel
            idx = idx+1
            msdis_scr(jdx) = msdis_scr(jdx) - msdis_2(map_inf(idx))
            gamdis_scr(jdx) = multd2h(gamdis_scr(jdx),
     &           gamdis_2(map_inf(idx)))
          end do
          idx = idx+1
          nel = map_inf(idx)
          err = err.or.nel.gt.1
          do iel = 1, nel
            idx = idx+1
            msdis_1(map_inf(idx)) = msdis_scr(jdx)
            gamdis_1(map_inf(idx)) = gamdis_scr(jdx) 
          end do
        end do
        if (err)
     &       call quit(1,'split_msgmdis','non-invertible map')
      end if
      
      return
      end
