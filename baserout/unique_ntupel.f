      subroutine unique_ntupel(ntupel,nel,nlist,nunique)
*
*     on input: list of N-tupels of integers ntupel(nel,nlist)
*
*     on output: removed all identical N-tupels, ntupel(nel,nunique)
*              contains contiguous list of unique N-tupels
*
      implicit none

      integer, intent(in) ::
     &       nel, nlist
      integer, intent(inout) ::
     &       ntupel(nel,nlist)
      integer, intent(out) ::
     &       nunique

      integer ::
     &       idx, jdx, kdx, iel
      logical ::
     &       same, unique


      kdx = 1
      do idx = 2, nlist
        
c dbg
c        print *,'present: ',idx,'->',ntupel(1:nel,idx)
c dbg
        unique = .true.
        do jdx = 1, kdx
c dbg
c        print *,'cmp to:  ',jdx,'->',ntupel(1:nel,jdx)
c dbg
          same = .true.
          do iel = 1, nel
            same = same.and.ntupel(iel,idx).eq.ntupel(iel,jdx) 
          end do
c dbg
c       print *,'same? ',same
c dbg
          unique = unique.and..not.same
        end do
c dbg
c       print *,'unique? ',unique
c dbg

        if (unique) then
          kdx = kdx+1
          if (kdx.ne.idx) then
            ntupel(1:nel,kdx) = ntupel(1:nel,idx)
          end if
        end if

      end do

      nunique = kdx

      end

