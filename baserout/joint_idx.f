*----------------------------------------------------------------------*
      integer function joint_idx(idx1,idx2,ibase,sort_it)
*----------------------------------------------------------------------*
*     pack two packed index lists (base ibase) into one
*----------------------------------------------------------------------*
      implicit none

      integer, intent(in) ::
     &     idx1, idx2, ibase
      logical, intent(in) ::
     &     sort_it

      integer ::
     &     nmax, itest, itest2, nexp1, nexp2
      integer, allocatable ::
     &     idx_exp(:)

      integer, external ::
     &     int_expand, int_pack

      ! test possible max. packing
      itest = ibase
      nmax = 1
      do
        itest2 = itest*ibase
        if (itest2.lt.itest) exit
        itest = itest2
        nmax = nmax+1
      end do

      allocate(idx_exp(2*nmax))

      nexp1 = int_expand(idx1,ibase,idx_exp)
      nexp2 = int_expand(idx2,ibase,idx_exp(nexp1+1))

      if (nexp1+nexp2.gt.nmax)
     &     call quit(1,'joint_idx','cannot join indices')
      ! sort new vertex list
      if (sort_it) call isort(idx_exp,nexp1+nexp2,+1)
      ! pack again
      joint_idx = int_pack(idx_exp,nexp1+nexp2,ibase)

      deallocate(idx_exp)

      return
      end
