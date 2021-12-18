
*----------------------------------------------------------------------*
      function iocc_overlap(iocc,dagi,jocc,dagj)
*----------------------------------------------------------------------*
*     return overlap (minimum common occupation per HPV/CA) 
*     include interface file ifc_ioccfunc.inc in calling routines!
*----------------------------------------------------------------------*
      
      implicit none
      include 'opdim.h'
      
      integer :: iocc_overlap(ngastp,2)

      logical, intent(in) ::
     &     dagi, dagj

      integer, intent(in) ::
     &     iocc(ngastp,2), jocc(ngastp,2)
      
      integer ::
     &     ica, ica_i, ica_j, ihpv
      integer ::
     &     iocc_scr(ngastp,2)

      do ica = 1,2
        ica_i = ica
        ica_j = ica
        if (dagi) ica_i = 3-ica
        if (dagj) ica_j = 3-ica
        do ihpv = 1,ngastp
          iocc_scr(ihpv,ica) =
     &         min(iocc(ihpv,ica_i),jocc(ihpv,ica_j))
        end do
      end do

      iocc_overlap = iocc_scr

      return
      end
