      integer function std_spsign_msdis(msdis,occ,nblk)
*
*     set up standard sign for standard spin-string arising from
*     a distribution of MS values over nblk blocks:
*
*     e.g.  msdis=(+2,-1,+1)
*           occ  =(2,3,1)
*
*     --> standard spin-string:  aa abb a
*
*     standard sign: determined by number of permutations to arrange to
*     a sequence with first all alphas and then all betas, in our example
*
*       aaaabb  -> +1
*     
*
      implicit none

      integer, intent(in) ::
     &     nblk, msdis(nblk), occ(nblk)

      integer ::
     &     nidx, iblk, idx, isp

      integer, pointer ::
     &     scr(:)

      integer, external ::
     &     std_spsign

      std_spsign_msdis = 1
      if (nblk.le.1) return

      nidx = 0
      do iblk = 1, nblk
        nidx = nidx+occ(iblk)
      end do
      allocate(scr(nidx))

      idx = 0
      do iblk = 1, nblk
        do isp = 1, (occ(iblk)+msdis(iblk))/2
          idx = idx+1
          scr(idx) = +1
        end do
        do isp = 1, (occ(iblk)-msdis(iblk))/2
          idx = idx+1
          scr(idx) = -1
        end do
      end do

      std_spsign_msdis = std_spsign(scr,nidx)

      deallocate(scr)
      
      return
      end
