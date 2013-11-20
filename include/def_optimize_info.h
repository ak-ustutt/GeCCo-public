      integer, parameter ::
     &       optinf_prc_direct = 0, ! future: diagonal not on file
     &       optinf_prc_file   = 1, ! the usual case
     &       optinf_prc_blocked = 2,! as used for R12
     &       optinf_prc_mixed = 3,  ! both 1 and 2
     &       optinf_prc_traf = 4,   ! via transformation to orth. basis 
     &       optinf_prc_norm = 5,   ! like usual, but normal. to 1
     &       optinf_prc_spinp = 6,  ! like 1, but with spin projection
     &       optinf_prc_invH0 = 7,  ! invert non-diagonal H0
     &       optinf_prc_prj = 8     ! apply (projection) formula

      ! input variables to control optimization
      type optimize_info
        logical ::
     &       variational,     ! variational energy
     &       linear,          ! linear
     &       skip_resx,       ! skip evaluation of not needed residuals
     &       resume           ! use last eigenvec. as initial guess
        integer ::
     &       nopt,            ! sets of operators descr. variables
     &                        ! to be optimized
     &       nroot,           ! number of roots to be solved for
     &       norder,          ! order of optimization
     &       mode_evp,
     &       mode_leq,
     &       mode_nleq,
     &       maxmacit, maxmicit, micifac,
     &       maxsbsp, max_incore,
     &       update_prc,      ! allow updating of preconditioner
     &       optref           ! optimization route for relaxed reference
        real(8) ::
     &       trini,
     &       shift,           ! shift linear equations by value
     &       mic_ahead        ! control weaker conv. thr. for micro it.
        integer, pointer ::   ! # of wave-function parameters (wfp's)
     &       nwfpar(:),                     ! dimension: nopt
     &       typ_prc(:),      ! type of preconditioner
     &       nsec(:),         ! # of sections in wfp list
     &       nwfpsec(:),      ! # of wfp's in section, dim=sum(nsec)
     &       idstsec(:)       ! first idx in section, dim=sum(nsec)
        real(8), pointer ::
     &       thrgrd(:),       ! dimension: nopt
     &       signsec(:)       ! contraction sign of section, dim=sum(nsec)
      end type optimize_info

      integer, parameter ::
     &     mode_leq_conjg = 0,
     &     mode_leq_subsp = 1
      integer, parameter ::
     &     mode_nleq_pert = 0,
     &     mode_nleq_diis = 1,
     &     mode_nleq_assj = 2,
     &     mode_nleq_2ndo = 3
          
