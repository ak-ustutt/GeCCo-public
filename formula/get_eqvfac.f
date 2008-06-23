      integer function get_eqvfac(neqv,fix_vtx,nvtx)

      integer, intent(in) ::
     &     nvtx,
     &     neqv(nvtx)
      logical, intent(in) ::
     &     fix_vtx(nvtx)

      integer ::
     &     ieqvfac, ivtx

      integer, external ::
     &     ifac

      ieqvfac = 1
      do ivtx = 1, nvtx
        if (fix_vtx(ivtx).or.neqv(ivtx).lt.0) cycle
        ieqvfac = ieqvfac*ifac(neqv(ivtx))
      end do
      
      get_eqvfac = ieqvfac

      return
      end
