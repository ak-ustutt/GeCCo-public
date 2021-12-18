*----------------------------------------------------------------------*
      subroutine did2msgm(msd,gmd,did,iocc,nsym,njoined)
*----------------------------------------------------------------------*
*     calculate (ms,Gamma) distribution from integer-valued Distrib. ID 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(in) ::
     &     did, iocc(ngastp,2,njoined), nsym, njoined
      integer, intent(out) ::
     &     msd(ngastp,2,njoined), gmd(ngastp,2,njoined)

      integer ::
     &     ipatt, ica, ihpvdx, ihpv, idxms, ijoin

      ipatt = did
      msd(1:ngastp,1:2,1:njoined) = 0
      gmd(1:ngastp,1:2,1:njoined) = 1
      do ica = 1,2
        do ihpvdx = 1, ngastp
          do ijoin = 1, njoined
            ihpv = hpvxseq(ihpvdx)
            if (iocc(ihpv,ica,ijoin).eq.0) cycle
            gmd(ihpv,ica,ijoin) = mod(ipatt,nsym)+1
            ipatt = ipatt/nsym
            idxms = mod(ipatt,iocc(ihpv,ica,ijoin)+1)+1
            msd(ihpv,ica,ijoin) = iocc(ihpv,ica,ijoin)-(idxms-1)*2
            ipatt = ipatt/(iocc(ihpv,ica,ijoin)+1)
          end do
        end do
      end do

      return
      end
*----------------------------------------------------------------------*
