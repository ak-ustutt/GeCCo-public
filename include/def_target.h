      integer, parameter ::
     &     len_target_name = 32,
     &     len_command_name = 32,
     &     len_command_par  = 512

      integer, parameter ::
     &     ttype_op  = 1,
     &     ttype_frm = 2,
     &     ttype_frmopt = 3,
     &     ttype_opme = 4,
     &     ttype_gen = 0

      type action

        integer ::
     &       type      ! type of targets modified by this action
        character(len_command_name) ::
     &       command   ! the command to execute
        integer ::
     &       n_update,  ! the first n_update targets will get updated 
     &       n_labels, n_parameter_strings
        character(len_target_name), pointer ::
     &       labels(:)      ! list of target labels
        character(len_command_par), pointer ::
     &       parameters(:)  ! list of parameters as formatted string
        
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
