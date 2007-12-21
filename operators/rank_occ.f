      integer function rank_occ(modestr,occ,njoined)

      implicit none

      include 'opdim.h'

      character(*), intent(in) ::
     &     modestr
      integer, intent(in) ::
     &     njoined
      integer, intent(in) ::
     &     occ(ngastp,2,njoined)

      integer ::
     &     igastp, ijoin

      select case(trim(modestr))
      case ('C-A')
        rank_occ = 0
        do ijoin = 1, njoined
          do igastp = 1, ngastp
            rank_occ = occ(igastp,1,ijoin) - occ(igastp,2,ijoin)
          end do
        end do
      case default
        call quit(1,'rank_occ','unknown modestr: ',trim(modestr))
      end select

      return
      end
