
*----------------------------------------------------------------------*
      function iocc_add(ifac,iocc,dagi,jfac,jocc,dagj)
*----------------------------------------------------------------------*
*     add two occupations
*     include interface file ifc_operators.h in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'
      
      integer :: iocc_add(ngastp,2)

      logical, intent(in) ::
     &     dagi, dagj

      integer, intent(in) ::
     &     ifac, jfac,
     &     iocc(ngastp,2), jocc(ngastp,2)
      
      integer ::
     &     ica, ica_i, ica_j, ihpv, iscr(ngastp,2)
      
      do ica = 1,2
        ica_i = ica
        ica_j = ica
        if (dagi) ica_i = 3-ica
        if (dagj) ica_j = 3-ica

        iscr(1:ngastp,ica) =
     &       ifac*iocc(1:ngastp,ica_i)+jfac*jocc(1:ngastp,ica_j)

      end do

      iocc_add = iscr

      return
      end
