*----------------------------------------------------------------------*
      subroutine update_reo(ireo,nel,
     &     idx1,idx2,idx_merge,imvleft,nmvleft)
*----------------------------------------------------------------------*
*     auxiliary routine: update reordering array [ireo(old) => new]
*     idx1 and idx2 are removed, instead idx_merge is inserted
*     imvleft(1:nmvleft) contains indices which must be left of 
*     idx_merge in the new ordering
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, idx1, idx2, idx_merge, nmvleft,
     &     imvleft(nmvleft)
      integer, intent(inout) ::
     &     ireo(nel)

      integer ::
     &     idx, ireo2(nel), idx_off, imv, inum, isub
      integer, external ::
     &     idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is update_reo')
        write(luout,'(x,a,10i5)') 'on input: ',ireo(1:nel)
        write(luout,*) 'idx1, idx2: ', idx1, idx2
        write(luout,*) 'idx_merge:  ',idx_merge
        write(luout,*) 'imvleft:    ',imvleft(1:nmvleft)
      end if

      if (idx2.le.idx1 .or. idx_merge.lt.idx1 .or. idx_merge.gt.idx2)
     &     call quit(1,'update_reo','error 1')

      ! set up new permutation vector
      do idx = 1, nel
        ireo2(idx) = idx
      end do

      do idx = idx1, idx_merge-1
        ireo2(idx) = ireo2(idx+1)
      end do
      ireo2(idx_merge) = idx1
      ireo2(idx2) = 0
c      do idx = idx2+1, nel
c        ireo2(idx) = ireo2(idx)-1
c      end do

      if (ntest.ge.100) then
        write(luout,*) 'ireo2: ',ireo2(1:nel)
      end if

      ! apply to ireo
      call perm_mult(ireo,ireo,ireo2,nel)

      if (ntest.ge.100) then
        write(luout,*) 'new ireo: ',ireo(1:nel)
      end if

      ! ensure 1-increment:
      ! e.g.  1 3 0 4  -->  1 2 0 3
      inum = 0
      isub = 0
      do while (inum.lt.nel)
        inum = inum+1
        idx = idxlist(inum,ireo,nel,1)
        if (idx.gt.0) then
          ireo(idx) = inum-isub
        else
          isub = isub+1
        end if
      end do

      if (ntest.ge.100) then
        write(luout,*) 'mod. ireo: ',ireo(1:nel)
      end if

      ! reset permutation vector ...
      do idx = 1, nel
        ireo2(idx) = idx
      end do

      ! an set up permutation for imvleft
      idx_off = idx_merge
      do imv = 1, nmvleft        
        do idx = idx_off+imv+1, imvleft(imv)
          ireo2(idx) = ireo2(idx-1)
        end do
        ireo2(idx_off+imv) = imvleft(imv)
      end do

      if (ntest.ge.100) then
        write(luout,*) 'ireo2: ',ireo2(1:nel)
      end if

      call perm_mult(ireo,ireo,ireo2,nel)

      if (ntest.ge.100) then
        write(luout,*) 'final ireo: ',ireo(1:nel)
      end if

      return
      end
