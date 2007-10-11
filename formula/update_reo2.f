*----------------------------------------------------------------------*
      subroutine update_reo2(ireo,iloweq,nel,
     &     idx1,idx2,idx_merge,imvleft,nmvleft)
*----------------------------------------------------------------------*
*     auxiliary routine: update reordering array [ireo(old) => new]
*     idx1 and idx2 are removed, instead idx_merge is inserted
*     iloweq(1:nel) holds the lowest equivalent vertex (old numbering)
*     i.e. the lowest vertex number contributing to a merged vertex
*
*     UNUSED, NECCESSARY?:
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
     &     ireo(nel), iloweq(nel)

      integer ::
     &     idx, ireo2(nel), idx1m, idx2m, idxmm,
     &     idxmv, imv
      integer, external ::
     &     idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is update_reo')
        write(luout,'(x,a,10i5)') 'on input: ',ireo(1:nel)
        write(luout,*) 'iloweq: ', iloweq(1:nel)
        write(luout,*) 'idx1, idx2: ', idx1, idx2
        write(luout,*) 'idx_merge:  ',idx_merge
        write(luout,*) 'imvleft:    ',imvleft(1:nmvleft)

        call perm_inv(ireo2,ireo,nel)
        write(luout,*) 'ireo reorders'
        write(luout,'(x,a,10i5)') 'this: ',(idx,idx=1,nel)
        write(luout,'(x,a,10i5)') 'to  : ',ireo2(1:nel)
        
      end if

c      if (nmvleft.gt.0) call quit(1,'update_reo','mvleft needed?')

      if (idx2.le.idx1 .or. idx_merge.lt.idx1 .or. idx_merge.gt.idx2)
     &     call quit(1,'update_reo','error 1')

      ! -------------------------------
      !  set up new permutation vector
      ! -------------------------------

      ! get lowest equivalent vertex (old ordering)
      idx1m = iloweq(idx1)
      idx2m = iloweq(idx2)
      idxmm = iloweq(idx_merge)
c dbg
c      print *,'a) idx1m, idx2m, idxmm: ',idx1m, idx2m, idxmm
c dbg      

      ! old ordering to latest ordering (def'd by input ireo)
c      idx1m = idxlist(idx1m,ireo,nel,1)
c      idx2m = idxlist(idx2m,ireo,nel,1)
c      idxmm = idxlist(idxmm,ireo,nel,1)
      idx1m = ireo(idx1m)
      idx2m = ireo(idx2m)
      idxmm = ireo(idxmm)
      if (idx1m.le.0.or.idx2m.le.0.or.idxmm.le.0)
     &     call quit(1,'update_reo','inconsistent ireo (input)')
c dbg
c      print *,'b) idx1m, idx2m, idxmm: ',idx1m, idx2m, idxmm
c dbg      

      ! init with unity:
      do idx = 1, nel
        ireo2(idx) = idx
      end do
      ! make permutation:
      do idx = idx1m, idxmm-1
        ireo2(idx) = ireo2(idx+1)
      end do
      ireo2(idxmm) = idx1m
      do idx = max(idx2m,idxmm), nel-1
        ireo2(idx) = ireo2(idx+1)
      end do
      ireo2(nel) = idx2m

      if (ntest.ge.100) then
        write(luout,*) 'ireo2 (ori): ',ireo2(1:nel)
      end if

      ! actually, we need the inverse, so:
      call perm_inv(ireo2,ireo2,nel)

      if (ntest.ge.100) then
        write(luout,*) 'ireo2(inv): ',ireo2(1:nel)
      end if

      ! apply to ireo
      ! I assume this is better:
      call perm_mult(ireo,ireo2,ireo,nel)
c     not this
c      call perm_mult(ireo,ireo,ireo2,nel)

      if (nmvleft.gt.0) then
        if (ntest.ge.100) then
          write(luout,*) 'ireo before mvleft: ',ireo(1:nel)
        end if

        do idx = 1, nel
          ireo2(idx) = idx
        end do        

        ! an set up permutation for imvleft
        !idx_off = idxmm !idx_merge
        do imv = 1, nmvleft        
          idxmv = ireo(iloweq(imvleft(imv)))
c dbg
c          print *,'imv, idxmv, idxmm: ',imv, idxmv, idxmm
c dbg
          do idx = idxmm+imv, idxmv
            ireo2(idx) = ireo2(idx-1)
          end do
          ireo2(idxmm+imv-1) = idxmv
        end do

        if (ntest.ge.100) then
          write(luout,*) 'ireo2: ',ireo2(1:nel)
        end if

        call perm_mult(ireo,ireo2,ireo,nel)

      end if

      ! update iloweq
      iloweq(idx1) = min(iloweq(idx1),idx1,idx2)
      iloweq(idx2) = min(iloweq(idx2),idx1,idx2)

      if (ntest.ge.100) then
        write(luout,*) 'final ireo: ',ireo(1:nel)
        write(luout,*) 'final iloweq: ',iloweq(1:nel)

        call perm_inv(ireo2,ireo,nel)
        write(luout,*) 'ireo reorders'
        write(luout,'(x,a,10i5)') 'this: ',(idx,idx=1,nel)
        write(luout,'(x,a,10i5)') 'to  : ',ireo2(1:nel)

      end if

      return
      end
