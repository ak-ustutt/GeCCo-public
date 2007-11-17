      subroutine set_strmapdim_c(nstr1str2,nstr1,nstr2,
     &                         nblk12,
     &                         nstr1r,nstr2r,map12)

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     nstr1str2(*), nstr1(*), nstr2(*)
      integer, intent(in) ::
     &     nblk12,
     &     nstr1r(*), nstr2r(*), map12(*)
            
      integer ::
     &     idxmap, iblk12, iblk, len1, len2, nblk

      idxmap = 0
      do iblk12 = 1, nblk12
c dbg
c        print *,'iblk12: ',iblk12
c        print *,'map12: ',map12(idxmap+1:idxmap+4)
c dbg
        idxmap = idxmap+1
        nblk = map12(idxmap) ! already adapted for arbitrary nblk>1
                             ! (not sure whether this works ...)
        len1 = 1
        do iblk = 1, nblk
          idxmap = idxmap+1
          len1 = len1*nstr1r(map12(idxmap))
        end do
        idxmap = idxmap+1
        len2 = 1
        nblk = map12(idxmap)
        do iblk = 1, nblk
          idxmap = idxmap+1
          len2 = len2*nstr2r(map12(idxmap))
        end do
        nstr1str2(iblk12) = len1*len2
        nstr1(iblk12) = len1
        nstr2(iblk12) = len2
c dbg
c        print *,'len1, len2: ',len1,len2
c dbg
      end do

      return
      end
