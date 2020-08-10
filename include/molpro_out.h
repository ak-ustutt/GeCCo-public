      integer ::
     &     ci_iter,     ! Number of iterations for the ci problem, set in solve_evp
     &     n_sv         ! Number of singular values retained
      real(8) ::
     &     sv_min,      ! Smallest singular value included
     &     sv_max       ! Largest singular value excluded
      logical ::
     &     lmol         ! True if using molpro interface
      common/molpro_out/ ci_iter, lmol

