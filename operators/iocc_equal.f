
*----------------------------------------------------------------------*
      logical function iocc_equal(iocc,dagi,jocc,dagj)
*----------------------------------------------------------------------*
*     compare two occupations and return .true. if equal
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
*----------------------------------------------------------------------*

      logical, intent(in) ::
     &     dagi, dagj
      integer, intent(in) ::
     &     iocc(ngastp,2), jocc(ngastp,2)

      integer ::
     &     ica, ica_i, ica_j, ihpv

      iocc_equal = .true.

      do ica = 1,2
        ica_i = ica
        ica_j = ica
        if (dagi) ica_i = 3-ica
        if (dagj) ica_j = 3-ica
        
        do ihpv = 1, ngastp
          iocc_equal = iocc_equal.and.
     &         iocc(ihpv,ica_i).eq.jocc(ihpv,ica_j)
        end do
        if (iocc_equal.eq..false.) exit
      end do

      return
      end
