      pure function get_nidx4contr(contr) result(nidx)

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'

      type(contraction), intent(in) :: contr
      integer :: iarc
      integer :: nidx

      nidx = 0
      do iarc = 1, contr%narc
        nidx = nidx+2*sum(contr%arc(iarc)%occ_cnt(1:ngastp,1:2))
      end do

      do iarc = 1, contr%nxarc
        nidx = nidx+sum(contr%xarc(iarc)%occ_cnt(1:ngastp,1:2))
      end do

      return
      
      end function
