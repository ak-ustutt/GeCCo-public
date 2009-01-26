      subroutine get_nice_arc_seq(arc_seq,contr)

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, intent(out) ::
     &     arc_seq(*)
      type(contraction), intent(in) ::
     &     contr

      integer ::
     &     narc, iarc, jarc, iarc0
      type(cntr_arc), pointer ::
     &     arc(:)

      narc = contr%narc
      arc => contr%arc

      ! init with ascending sequence
      do iarc = 1, narc
        arc_seq(iarc) = iarc
      end do

      do iarc = 2, narc
        iarc0 = arc_seq(iarc)
        jarc = iarc-1        
        do while (jarc.gt.0.and.
     &            .not.correct_arc_seq(arc(arc_seq(max(jarc,1))),
     &            arc(iarc0)))
          arc_seq(jarc+1) = arc_seq(jarc)
          jarc = jarc-1
        end do
        arc_seq(jarc+1) = iarc0
      end do

      return

      contains

      logical function correct_arc_seq(arc1,arc2)

      type(cntr_arc) ::
     &     arc1, arc2

      integer ::
     &     n1, n2
      logical ::
     &     even1, even2

      n1 = sum(arc1%occ_cnt(1:ngastp,1:2))
      n2 = sum(arc2%occ_cnt(1:ngastp,1:2))

      even1 = mod(n1,2).eq.0
      even2 = mod(n2,2).eq.0

      correct_arc_seq = even1.or..not.even2

      if (.not.correct_arc_seq.or.even1.neqv.even2) return

      correct_arc_seq = n1.ge.n2

      return

      end function

      end
