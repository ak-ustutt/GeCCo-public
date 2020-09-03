      integer ::
     &     ci_iter,     ! Number of iterations for the ci problem, set in solve_evp
     &     n_sv         ! Number of singular values retained
      real(8) ::
     &     sv_min,      ! Smallest singular value included
     &     sv_max       ! Largest singular value excluded
      logical ::
     &     lmol,        ! True if using molpro interface
     &     no_print     ! True if inside MRCC main loop (don't print out the ci iterations)
      common/molpro_out/ ci_iter,n_sv,sv_min,sv_max,lmol,no_print

