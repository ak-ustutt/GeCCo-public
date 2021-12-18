*----------------------------------------------------------------------*
      subroutine set_dis_tra_map(mapca,mapac,
     &     hpvxc,hpvxa,nc,na)
*----------------------------------------------------------------------*
*     set the map array for tranpositions of MS or GAMMA distributions
*----------------------------------------------------------------------*

      implicit none
      
      include 'opdim.h'
      include 'hpvxseq.h' ! <- determines the sequence of HPVX

      integer, intent(in) ::
     &     nc, na,
     &     hpvxc(nc), hpvxa(na)
      integer, intent(out) ::
     &     mapca(nc), mapac(na)

      integer ::
     &     idx, idxa, idxc, idx_hpvx, hpvx

      idxc = 0
      idxa = 0
      do idx_hpvx = 1, ngastp
        hpvx = hpvxseq(idx_hpvx)        
          
        do idx = nc, 1, -1
          if (hpvxc(idx).ne.hpvx) cycle
          idxc = idxc+1
          mapca(idx) = idxc
        end do
        do idx = na, 1, -1
          if (hpvxa(idx).ne.hpvx) cycle
          idxa = idxa+1
          mapac(idx) = idxa
        end do

      end do

      return
      end
      
