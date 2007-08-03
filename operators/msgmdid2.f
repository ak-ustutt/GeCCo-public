*----------------------------------------------------------------------*
      pure integer function msgmdid2(occ_c,idxms_c,gam_c,nc,
     &                               occ_a,idxms_a,gam_a,na,nsym)
*----------------------------------------------------------------------*
*     calculate integer-valued ID of (ms,Gamma) distribution
*     version for condensed representation
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'

      integer, intent(in) ::
     &     nc, na, nsym,
     &     occ_c(nc), idxms_c(nc), gam_c(nc),
     &     occ_a(na), idxms_a(na), gam_a(na) 

      integer ::
     &     ipatt, ibase, idx

      ipatt = 0
      ibase = 1
      do idx = 1, nc
        ipatt = ipatt + ((idxms_c(idx)-1)*nsym + gam_c(idx)-1)*ibase
        ibase = ibase * (occ_c(idx)+1)*nsym
      end do
      do idx = 1, na
        ipatt = ipatt + ((idxms_a(idx)-1)*nsym + gam_a(idx)-1)*ibase
        ibase = ibase * (occ_a(idx)+1)*nsym
      end do

      msgmdid2 = ipatt

      return
      end
