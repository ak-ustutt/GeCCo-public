      integer function len_blk_est(occ,nj,orb_info)
      
      implicit none

      include 'opdim.h'
      include 'def_orbinf.h'

      integer, intent(in) ::
     &   nj
      integer, intent(in) ::
     &   occ(ngastp,2,nj)
      type(orbinf) ::
     &   orb_info

      integer ::
     &   nhole, nvale, npart, nextr
      integer, external ::
     &   len_est__

      nhole = orb_info%nactt_hpv(IHOLE)
      nvale = orb_info%nactt_hpv(IVALE)
      npart = orb_info%nactt_hpv(IPART)
      nextr = orb_info%nactt_hpv(IEXTR)
              
      len_blk_est = len_est__(nhole,sum(occ(IHOLE,1:2,1:nj)))
     &             *len_est__(nvale,sum(occ(IVALE,1:2,1:nj)))
     &             *len_est__(npart,sum(occ(IPART,1:2,1:nj)))
     &             *len_est__(nextr,sum(occ(IEXTR,1:2,1:nj)))

      end

      integer function len_est__(norb,nidx)

      integer, intent(in) ::
     &    norb,nidx
      integer ::
     &    ii

      len_est__ = 1
      do ii = norb, max(1,norb-nidx+1), -1
        len_est__ = len_est__*ii
      end do

      end 

