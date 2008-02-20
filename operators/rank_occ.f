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
     &     igastp, ijoin, rank_occ_tmp

      select case(trim(modestr))
      case ('C-A')
        rank_occ = 0
        do ijoin = 1, njoined
          do igastp = 1, ngastp
            rank_occ = rank_occ
     &           +occ(igastp,1,ijoin) - occ(igastp,2,ijoin)
          end do
        end do
      case ('C')
        rank_occ = 0
        do ijoin = 1, njoined
          do igastp = 1, ngastp
            rank_occ = rank_occ+occ(igastp,1,ijoin)
          end do
        end do
      case ('A')
        rank_occ = 0
        do ijoin = 1, njoined
          do igastp = 1, ngastp
            rank_occ = rank_occ+occ(igastp,2,ijoin)
          end do
        end do
      case ('X')
        rank_occ = 0
        do ijoin = 1, njoined
          rank_occ = rank_occ+
     &        -occ(IHOLE,1,ijoin)
     &        +occ(IPART,1,ijoin)
     &        +occ(IEXTR,1,ijoin)
     &        +occ(IHOLE,2,ijoin)
     &        -occ(IPART,2,ijoin)
     &        -occ(IEXTR,1,ijoin)
        end do
        rank_occ = rank_occ/2
      case default
        call quit(1,'rank_occ','unknown modestr: ',trim(modestr))
      end select

      return
      end
