      integer, parameter ::
     &       optinf_prc_direct = 0, ! future: diagonal not on file
     &       optinf_prc_file   = 1, ! the usual case
     &       optinf_prc_blocked = 2,! as used for R12
     &       optinf_prc_mixed = 3   ! both 1 and 2

      ! input variables to control optimization
      type optimize_info
        logical ::
     &       variational,     ! variational energy
     &       linear           ! linear
        integer ::
     &       nopt,            ! sets of operators descr. variables
     &                        ! to be optimized
     &       nroot,           ! number of roots to be solved for
     &       norder,          ! order of optimization
     &       mode_evp,
     &       mode_leq,
     &       mode_nleq,
     &       maxmacit, maxmicit, micifac,
     &       maxsbsp, max_incore
        real(8) ::
     &       trini,
     &       shift            ! shift linear equations by value
        integer, pointer ::   ! # of wave-function parameters (wfp's)
     &       nwfpar(:),                     ! dimension: nopt
     &       typ_prc(:)       ! type of preconditioner
        real(8), pointer ::
     &       thrgrd(:)        ! dimension: nopt
      end type optimize_info

      integer, parameter ::
     &     mode_leq_conjg = 0,
     &     mode_leq_subsp = 1
      integer, parameter ::
     &     mode_nleq_pert = 0,
     &     mode_nleq_diis = 1,
     &     mode_nleq_assj = 2,
     &     mode_nleq_2ndo = 3
          
