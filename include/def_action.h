      integer, parameter ::
     &     iaction_import     = 1,
     &     iaction_evaluate   = 2,
     &     iaction_solve_leq  = 3,
     &     iaction_solve_nleq = 4,
     &     iaction_solve_evp  = 5,
     &     iaction_solve_gevp = 6

      type action
        integer ::
     &       action_type
        integer ::
     &       nop_in, nop_out
        integer, allocatable ::
     &       idxopdef_in(:), idxopdef_out(:),
     &       idxopfile_in(:,:), idxopfile_out(:,:)
        integer ::
     &       idx_formula
        type(filinf) ::
     &       fform_opt
      end type action
