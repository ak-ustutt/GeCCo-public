      integer, parameter ::
     &     len_target_name = 32,
     &     len_command_name = 32,
     &     len_command_par  = 4096,!2048,!1024,
     &     len_str_batch = 32

      integer, parameter ::
     &     ttype_op  = 1,
     &     ttype_frm = 2,
     &     ttype_frmopt = 3,
     &     ttype_opme = 4,
     &     ttype_gen = 0

      integer, parameter ::
     &     aatype_label = 0,
     &     aatype_log = 1,
     &     aatype_int = 2,
     &     aatype_occ = 3,
     &     aatype_restr = 4,
     &     aatype_rl8 = 5,
     &     aatype_str = 6

      type action_arg
        integer ::
     &       type
        integer ::
     &       arg_dim, n_str_batch
        character(len_target_name) ::
     &       arg_label
        character(len_target_name), pointer ::
     &       val_label(:)
        logical, pointer ::
     &       val_log(:)
        integer, pointer ::
     &       val_int(:)
        integer, pointer ::
     &       val_occ(:,:,:)
        integer, pointer ::
     &       val_restr(:,:,:,:,:,:)
        real(8), pointer ::
     &       val_rl8(:)
        character(len=len_str_batch), pointer ::
     &       val_str(:)
        logical ::
     &       required, def_provided

      end type action_arg

      type action

        integer ::
     &       type      ! type of targets modified by this action
        character(len_command_name) ::
     &       command   ! the command to execute
        logical ::
     &       new       ! new behaviour?
        ! needed for old behaviour
        integer ::
     &       n_update,  ! the first n_update targets will get updated 
     &       n_labels, n_parameter_strings        
        character(len_target_name), pointer ::
     &       labels(:)      ! list of target labels
        character(len_command_par), pointer ::
     &       parameters(:)  ! list of parameters as formatted string
        ! needed for new behaviour
        integer ::
     &       n_arguments ! number of arguments
        type(action_arg), pointer ::
     &       arg(:)

      end type action

      type target
      
        ! THIS TARGET:
        integer ::
     &       type     ! type of this target
        logical ::
     &       required ! .true.: make, even if not needed by other targets
        integer ::
     &       my_idx     ! index on array in target info
        integer(8) ::
     &       last_mod   ! last modification (should be identical to
                        ! value in target_info%last_mod
        character*(len_target_name) ::
     &       name
                        ! name of this target, must be valid name
                        ! that can be found on operator/formula list        
        ! DEPENDENCIES:
        integer ::
     &       n_joined_with, n_depends_on ! dimensions for arrays below:
        character*(len_target_name), pointer ::
     &       joined_with(:), ! targets, which are simultaneously evaluated;
     &       depends_on(:)    ! targets, on which this target depends;
        integer, pointer ::  ! the same info resolved to indices in
     &       idx_joined_with(:), ! target list
     &       idx_depends_on(:)  
        ! RULES:
        integer ::
     &       n_rules
        type(action), pointer ::
     &       rules(:)

      end type target
