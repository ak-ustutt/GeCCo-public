*----------------------------------------------------------------------*
      subroutine set_strmapdim_c2(nstr1str2,nstr1,nstr2,
     &                         ireo1,ireo2,
     &                         nblk12,
     &                         nstr1r,nstr2r,map12)
*----------------------------------------------------------------------*
*     set up length of strings and mappings for 
*      str12->str1,str2 resoultion for each block of str12 
*     (there are nblk12 blocks)
*     ireo1,ireo2 contain the info to which block of str12 the
*     given str1 or str2 block contributes
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'

      integer, intent(out) ::
     &     nstr1str2(*), nstr1(*), nstr2(*),
     &     ireo1(*), ireo2(*)
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
          ireo1(map12(idxmap)) = iblk12
        end do
        idxmap = idxmap+1
        len2 = 1
        nblk = map12(idxmap)
        do iblk = 1, nblk
          idxmap = idxmap+1
          len2 = len2*nstr2r(map12(idxmap))
          ireo2(map12(idxmap)) = iblk12
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
