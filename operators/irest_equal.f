*----------------------------------------------------------------------*
      logical function irestr_equal(irst_op1,dagger1,
     &                              irst_op2,dagger2,ngas)
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
      logical, intent(in) ::
     &     dagger1, dagger2

      logical ::
     &     same, reverse
      integer ::
     &     igas, ica1, ica2, imask

      reverse = dagger1.eq..not.dagger2

      same = .true.
      cmp_loop: do imask = 1, 2
        do ica1 = 1, 2
          ica2 = ica1
          if (reverse) ica2 = 3-ica1
          do igas = 1, ngas
            same = same.and.
     &           (irst_op1(1,igas,ica1,imask).eq.
     &            irst_op2(1,igas,ica2,imask)).and.
     &           (irst_op1(2,igas,ica1,imask).eq.
     &            irst_op2(2,igas,ica2,imask))
            if (.not.same) exit cmp_loop
          end do
        end do
      end do cmp_loop

      irestr_equal = same

      return
      end
