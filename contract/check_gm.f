*----------------------------------------------------------------------*
      pure logical function check_gm(gm,occ,nel)
*----------------------------------------------------------------------*
*     check whether IRREP values on gm conform with max occ. on occ
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel, gm(nel), occ(nel)

      integer ::
     &     iel

      check_gm = .true.
      do iel = 1, nel
        check_gm = check_gm.and.(gm(iel).eq.1.or.occ(iel).gt.0)
      end do

      return
      end
c*----------------------------------------------------------------------*
c      logical function check_gm(gmt,gm,nel)
c*----------------------------------------------------------------------*
c*     check whether IRREPS on gm conform with total IRREP gmt
c*----------------------------------------------------------------------*
c      implicit none
c
c      include 'multd2h.h'
c
c      integer, intent(in) ::
c     &     nel, gm(nel), gmt
c
c      integer ::
c     &     iel, gm_check
c
c      gm_check = 1
c      do iel = 1, nel
c        gm_check = multd2h(gm_check,gm(iel))
c      end do
c      check_gm = gm_check.eq.gmt
c
c      return
c      end
