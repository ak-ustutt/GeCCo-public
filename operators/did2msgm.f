*----------------------------------------------------------------------*
      subroutine did2msgm(msd,gmd,did,iocc,nsym)
*----------------------------------------------------------------------*
*     calculate (ms,Gamma) distribution from integer-valued Distrib. ID 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     msd(ngastp,2), gmd(ngastp,2)
      integer, intent(in) ::
     &     did, iocc(ngastp,2), nsym

      integer ::
     &     ipatt, ica, ihpvdx, ihpv, idxms

      ipatt = did
      msd(1:ngastp,1:2) = 0
      gmd(1:ngastp,1:2) = 1
      do ica = 1,2
        do ihpvdx = 1, ngastp
          ihpv = hpvxseq(ihpvdx)
          if (iocc(ihpv,ica).eq.0) cycle
          gmd(ihpv,ica) = mod(ipatt,nsym)+1
          ipatt = ipatt/nsym
          idxms = mod(ipatt,iocc(ihpv,ica)+1)+1
          msd(ihpv,ica) = iocc(ihpv,ica)-(idxms-1)*2
          ipatt = ipatt/(iocc(ihpv,ica)+1)
        end do
      end do

      return
      end
*----------------------------------------------------------------------*
