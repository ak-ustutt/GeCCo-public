*----------------------------------------------------------------------*
      subroutine set_spprjmap_cases(iocc,idspn_prm)
*----------------------------------------------------------------------*
*     loops over all possible permutations of a simple spin string
*     and assignes a number to each.
*     the reference permutation (with alternating + and -) is el. #0.
*     the permutation numbers are written to idspn_prm(0,0:2**iocc-1),
*     and the spin strings are written to idspn_prm(1:iocc,0:2**iocc-1)
*
*     matthias, feb 2014
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     iocc
      integer, intent(out) ::
     &     idspn_prm(0:iocc,0:2**iocc-1)

      integer
     &     icnt, icase, idigit

      icase = 0
      do icnt = 0, 2**iocc-1
        ! check if this is the reference element:
        ! +, -+, +-+, -+-+, ...
        ! the reference element has a + (0 bit) on the right
        ! and alternating + and - (0 and 1 bits)
        if (.not.btest(icnt,0).and.
     &      ieor(icnt,ishftc(icnt,-1,iocc)).eq.2**((iocc/2)*2)-1) then
          idspn_prm(0,0) = icnt
          do idigit = 1, iocc
            idspn_prm(idigit,0) = -2*ibits(icnt,iocc-idigit,1)+1
          end do
          cycle
        end if
        icase = icase+1
        do idigit = 1, iocc
          idspn_prm(idigit,icase) = -2*ibits(icnt,iocc-idigit,1)+1
        end do
        idspn_prm(0,icase) = icnt
      end do

      return
      end
