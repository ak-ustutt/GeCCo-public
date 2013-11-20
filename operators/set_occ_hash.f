*----------------------------------------------------------------------*
      subroutine set_occ_hash(occ_hash,op)
*----------------------------------------------------------------------*
*     set up hash table for finding operator block for given occupation
*     currently for 4 indices, works only for CCAA case.
*     needed in import routines to avoid excessive calls to iocc_blk
*     which is much too slow
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_operator.h'

      type(operator), intent(in) ::
     &     op
      integer, intent(out) ::
     &     occ_hash((ngastp*(ngastp+1)/2)**2)

      integer ::
     &     nj, nblk, iblk, iblkoff, idx,
     &     ica, hpvx, ij, nocc, ll, ihash,
     &     list(4)
      integer, pointer ::
     &     occ(:,:,:)

      nblk = op%n_occ_cls
      nj   = op%njoined

      occ => op%ihpvca_occ

      occ_hash(1:) = -1

      ! loop over operator occupations
      do iblk = 1, nblk

        if (op%ica_occ(1,iblk).ne.2) cycle
        if (op%formal_blk(iblk)) cycle

        iblkoff = (iblk-1)*nj 

        ! get hpvx of first 4 indices
        idx = 0
        outer : do ica = 1, 2
          do hpvx = 1, ngastp
            do ij = 1, nj
              nocc = occ(hpvx,ica,iblkoff+ij)
              if (nocc.eq.0) cycle
              do ll = 1, nocc
                idx = idx+1
                list(idx) = hpvx
                if (idx.eq.4) exit outer
              end do
            end do
          end do
        end do outer

        ihash = ((list(2)-1)*list(2)/2+list(1)-1)*(ngastp*(ngastp+1)/2)+
     &           (list(4)-1)*list(4)/2+list(3)
        occ_hash(ihash) = iblk

      end do 
        
      if (ntest.ge.100) then
        write(lulog,*) 'hash table:'
        call iwrtma(occ_hash,ngastp*(ngastp+1)/2,ngastp*(ngastp+1)/2,
     &                       ngastp*(ngastp+1)/2,ngastp*(ngastp+1)/2)

      end if

      end
