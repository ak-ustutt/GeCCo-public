*----------------------------------------------------------------------*
      integer function get_lenmap(lstr1,lstr2,map12,nblk12)
*----------------------------------------------------------------------*
*     get length of str1-str2 mapping
*     lstr1(*), lstr2(*) contain the lengthes of non-zero blocks
*     map12 contains information how the blocks of 1 and 2 are related
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     nblk12, map12(*), lstr1(*), lstr2(*)

      integer ::
     &     idxmap, iblk, nblk, len, iblk12

      idxmap = 0
      get_lenmap = 0
      do iblk12 = 1, nblk12
        idxmap = idxmap+1
        nblk = map12(idxmap) ! already adapted for arbitrary nblk>1
        len = 1
        do iblk = 1, nblk
          idxmap = idxmap+1
          len = len*lstr1(map12(idxmap))
        end do
        idxmap = idxmap+1
        nblk = map12(idxmap)
        do iblk = 1, nblk
          idxmap = idxmap+1
          len = len*lstr2(map12(idxmap))
        end do
        get_lenmap = get_lenmap+len
      end do

      return
      end
