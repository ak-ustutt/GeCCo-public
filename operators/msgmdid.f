*----------------------------------------------------------------------*
      integer function msgmdid(iocc,msd,gmd,nsym)
*----------------------------------------------------------------------*
*     calculate integer-valued ID of (ms,Gamma) distribution
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(in) ::
     &     iocc(ngastp,2), msd(ngastp,2), gmd(ngastp,2), nsym

      integer ::
     &     ipatt, ibase, ica, ihpvdx, ihpv

      ipatt = 0
      ibase = 1
      do ica = 1,2
        do ihpvdx = 1, ngastp
          ihpv = hpvxseq(ihpvdx)
          if (iocc(ihpv,ica).eq.0) cycle
          ipatt = ipatt +
     &        (((iocc(ihpv,ica)-msd(ihpv,ica))/2)*nsym+
     &                                  gmd(ihpv,ica)-1)*ibase
          ibase = ibase*(iocc(ihpv,ica)+1)*nsym
        end do
      end do

      msgmdid = ipatt

      return
      end
