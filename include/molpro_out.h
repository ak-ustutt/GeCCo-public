      integer ::
     &     ci_iter,     ! Number of iterations for the ci problem, set in solve_evp
     &     n_sv,        ! Number of singular values retained
     &     n_sv_um      ! Number of singular values retained, from update_matrix.f
      real(8) ::
     &     sv_min,      ! Smallest singular value included
     &     sv_max,      ! Largest singular value excluded
     &     sv_min_um,   ! Smallest singular value included, from update_matrix.f
     &     sv_max_um    ! Largest singular value excluded, from update_matrix.f
      logical ::
     &     lmol         ! True if using molpro interface
      common/molpro_out/ ci_iter,n_sv,n_sv_um,sv_min,sv_max,sv_min_um,
     &                   sv_max_um,lmol

