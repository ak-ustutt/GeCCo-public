*----------------------------------------------------------------------*
      subroutine idx_reform(idxprqs,igmprqs,ispprqs,
     &                    idx_in,nidx,nx_max,del_tri,
     &                    reord,igamorb,igasorb,hpvx_gas,norb)
*----------------------------------------------------------------------*
*     idx_in contains indices in the order of the external integral
*     program, in mulliken order (pq|rs)
*     reorder to GeCCo ordering (using reord) and store in Dirac
*     sequence <pr|qs>
*     also, provide symmetry and subspace number of orbital
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'

      logical, intent(in) ::
     &     del_tri
      integer, intent(in) ::
     &     nidx, nx_max, norb
      integer(2), intent(out) ::
     &     idxprqs(4,nidx), igmprqs(4,nidx), ispprqs(4,nidx)
      integer(2), intent(in) ::
     &     idx_in(4,nidx)
      integer, intent(in) ::
     &     reord(norb), igamorb(norb), igasorb(norb), hpvx_gas(ngastp)


      integer ::
     &     iidx, nx

*======================================================================*
      integer ::
     &     idxpq, ld
      integer(2) ::
     &     p,q

      idxpq(p,q,ld) = (min(p,q)-1)*ld+max(p,q)
*======================================================================*

      do iidx = 1, nidx
        ! dirac         mulliken
        idxprqs(1,iidx) = reord(idx_in(1,iidx))
        idxprqs(3,iidx) = reord(idx_in(2,iidx))
        idxprqs(2,iidx) = reord(idx_in(3,iidx))
        idxprqs(4,iidx) = reord(idx_in(4,iidx))
      end do

      do iidx = 1, nidx
        igmprqs(1,iidx) = igamorb(idxprqs(1,iidx))
        igmprqs(2,iidx) = igamorb(idxprqs(2,iidx))
        igmprqs(3,iidx) = igamorb(idxprqs(3,iidx))
        igmprqs(4,iidx) = igamorb(idxprqs(4,iidx))
      end do

      do iidx = 1, nidx
        ispprqs(1,iidx) = igasorb(idxprqs(1,iidx))
        ispprqs(2,iidx) = igasorb(idxprqs(2,iidx))
        ispprqs(3,iidx) = igasorb(idxprqs(3,iidx))
        ispprqs(4,iidx) = igasorb(idxprqs(4,iidx))
      end do

      ! remove "upper" triangle, i.e. pr > qs
      if (del_tri) then
        do iidx = 1, nidx
          if (idxpq(idxprqs(1,iidx),idxprqs(2,iidx),norb).gt.
     &         idxpq(idxprqs(3,iidx),idxprqs(4,iidx),norb)) 
     &         idxprqs(1,iidx) = 0
        end do
      end if

      if (nx_max.ge.4) return

      ! rise flag for integrals to be ignored ...
      do iidx = 1, nidx
        nx = 0
        if (hpvx_gas(ispprqs(1,iidx)).eq.IEXTR) nx = nx+1
        if (hpvx_gas(ispprqs(2,iidx)).eq.IEXTR) nx = nx+1
        if (hpvx_gas(ispprqs(3,iidx)).eq.IEXTR) nx = nx+1
        if (hpvx_gas(ispprqs(4,iidx)).eq.IEXTR) nx = nx+1
        if (nx.gt.nx_max) idxprqs(1,iidx) = 0
      end do

      return
      end
