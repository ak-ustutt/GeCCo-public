      integer, parameter ::
     &     iaction_import     = 1,
     &     iaction_evaluate   = 2,
     &     iaction_solve_leq  = 3,
     &     iaction_solve_nleq = 4,
     &     iaction_solve_evp  = 5,
     &     iaction_solve_gevp = 6,
     &     iaction_setup_prc  = 7,
     &     iaction_prop_eval  = 8,
     &     iaction_symmetrise = 9,
c dbg
     &     iaction_diagonal   = 10
c dbg
      type action
        integer ::
     &       action_type
        integer ::
     &       nop_in, nop_out, nop_opt
        integer, allocatable ::
     &       idxopdef_in(:), idxopdef_out(:)
        integer ::
     &       nform
        integer, allocatable ::
     &       idx_formula(:)
        type(filinf) ::
     &       fform_opt
      end type action
