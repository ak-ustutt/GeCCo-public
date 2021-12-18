*----------------------------------------------------------------------*
      integer function ieqfac(itlist,icidxlist,nel)
*----------------------------------------------------------------------*
*     obtain (inverse of) prefactor for equivalent contractions
*     i.e. contraction index and operator index (itlist) are equal
*     itlist is supposed to be ordered, and icidxlist is ordered
*     for equivalent operators
*     VERY OLD AND OBSOLETE ROUTINE
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nel, itlist(nel), icidxlist(nel)

      integer ::
     &     itlast, iclast, idx, neq
      integer, external ::
     &     ifac

      if (ntest.ge.100) then
        write(lulog,*) '----------------'
        write(lulog,*) ' ieqfac speaks:'
        write(lulog,*) '----------------'
        write(lulog,*) ' itlist    : ',itlist(1:nel)
        write(lulog,*) ' icidxlist : ',icidxlist(1:nel)
      end if

      itlast = itlist(1)
      iclast = icidxlist(1)
      neq = 1
      ieqfac = 1
      do idx = 2, nel
        if (itlist(idx).eq.itlast) then
          if (icidxlist(idx).eq.iclast) then
            neq = neq+1
          else
            if (neq.gt.1) ieqfac = ieqfac*ifac(neq)
            neq = 1
            iclast = icidxlist(idx)
          end if
        end if
        if (itlist(idx).ne.itlast.or.idx.eq.nel) then
          if (neq.gt.1) ieqfac = ieqfac*ifac(neq)
          neq = 1
          itlast = itlist(idx)
          iclast = icidxlist(idx)
        end if
      end do

      if (ntest.ge.100) then
        write(lulog,*) ' result: ',ieqfac
      end if

      return
      end
