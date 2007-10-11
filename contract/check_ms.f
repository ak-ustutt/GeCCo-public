*----------------------------------------------------------------------*
      pure logical function check_ms(ms,occ,nel)
*----------------------------------------------------------------------*
*     check whether ms values on ms conform with max occ. on occ
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nel, ms(nel), occ(nel)

      integer ::
     &     iel

      check_ms = .true.
      do iel = 1, nel
        check_ms = check_ms.and.abs(ms(iel)).le.occ(iel)
      end do

      return
      end
