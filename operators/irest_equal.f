*----------------------------------------------------------------------*
      logical function irestr_equal(irst_op1,irst_op2,ngas)
*----------------------------------------------------------------------*
*
*     compare the restrictions irst_op1 and irst_op2
*
*     return true if equal
*
*----------------------------------------------------------------------*

      implicit none
c      include 'opdim.h'

      integer, intent(in) ::
     &     ngas, irst_op1(2,ngas,2,2), irst_op2(2,ngas,2,2)

      logical ::
     &     same
      integer ::
     &     igas, ica, imask

      same = .true.
      cmp_loop: do imask = 1, 2
        do ica = 1, 2 
          do igas = 1, ngas
            same = same.and.
     &           (irst_op1(1,igas,ica,imask).eq.
     &            irst_op2(1,igas,ica,imask)).and.
     &           (irst_op1(2,igas,ica,imask).eq.
     &            irst_op2(2,igas,ica,imask))
            if (.not.same) exit cmp_loop
          end do
        end do
      end do cmp_loop

      irestr_equal = same

      return
      end
