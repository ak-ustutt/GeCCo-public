      logical function occ_bound_p(cbound,occ1,occ2)

      implicit none

      include 'opdim.h'

       character, intent(in) ::
     &     cbound*(*)
      integer(8), intent(in) ::
     &     occ1, occ2

      integer(8) ::
     &     scr1, scr2
      integer ::
     &     iel


      occ_bound_p = .false.
      scr1 = occ1
      scr2 = occ2
      select case(trim(cbound))
      case('<','LT','lt')
        do iel = 1, ngastp*2
          if(mod(scr1,pack_base).ge.mod(scr2,pack_base)) exit
          scr1 = scr1 / pack_base
          scr2 = scr2 / pack_base
          occ_bound_p = iel.eq.ngastp*2
        end do
      case('<=','LE','le')
        do iel = 1, ngastp*2
          if(mod(scr1,pack_base).gt.mod(scr2,pack_base)) exit
          scr1 = scr1 / pack_base
          scr2 = scr2 / pack_base
          occ_bound_p = iel.eq.ngastp*2
        end do
      case('>','GT','gt')
        do iel = 1, ngastp*2
          if(mod(scr1,pack_base).le.mod(scr2,pack_base)) exit
          scr1 = scr1 / pack_base
          scr2 = scr2 / pack_base
          occ_bound_p = iel.eq.ngastp*2
        end do
      case('>=','GE','ge')
        do iel = 1, ngastp*2
          if(mod(scr1,pack_base).lt.mod(scr2,pack_base)) exit
          scr1 = scr1 / pack_base
          scr2 = scr2 / pack_base
          occ_bound_p = iel.eq.ngastp*2
        end do
      case default
        call quit(1,'occ_bound_p','Unknown relation: "'//cbound//'"')
      end select

      return
      end
