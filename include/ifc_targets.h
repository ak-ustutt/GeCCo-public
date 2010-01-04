      interface
        subroutine set_arg(name_target,command,
     &                                      arg_label,arg_dim,tgt_info,
     &     val_label,val_log,val_int,val_occ,val_restr,val_rl8,val_str)
        import
        implicit none
        type(target_info), intent(inout), target ::
     &       tgt_info
        character(len=*), intent(in) ::
     &       name_target, command, arg_label
        integer, intent(in) ::
     &       arg_dim
        character(len=*), intent(in), optional ::
     &       val_label(arg_dim)
        logical, intent(in), optional ::
     &       val_log(arg_dim)
        integer, intent(in), optional ::
     &       val_int(arg_dim)
        integer, intent(in), optional ::
     &       val_occ(ngastp,2,arg_dim)
        integer, intent(in), optional ::
     &       val_restr(:,:,:,:,:,:)
        real(8), intent(in), optional ::
     &       val_rl8(arg_dim)
        character(len=*), intent(in), optional ::
     &       val_str
        end subroutine
	  
        subroutine get_arg(arg_label,rule,tgt_info,
     &     val_label,val_label_list,val_log,val_log_list,
     &     val_int,val_int_list,val_occ,val_restr,
     &     val_rl8,val_rl8_list,val_str,ndim)
        import
        type(action), intent(in), target ::
     &       rule
        type(target_info), intent(in) ::
     &       tgt_info
        character(len=*), intent(in) ::
     &       arg_label
        character(len=*), intent(out), optional ::
     &       val_label
        character(len=*), intent(out), optional ::
     &       val_label_list(:)
        logical, intent(out), optional ::
     &       val_log
        logical, intent(out), optional ::
     &       val_log_list(:)
        integer, intent(out), optional ::
     &       val_int
        integer, intent(out), optional ::
     &       val_int_list(:)
        integer, intent(out), optional ::
     &       val_occ(:,:,:)
        integer, intent(out), optional ::
     &       val_restr(:,:,:,:,:,:)
        real(8), intent(out), optional ::
     &       val_rl8
        real(8), intent(out), optional ::
     &       val_rl8_list(:)
        character(len=*), intent(out), optional ::
     &       val_str
        integer, intent(out), optional ::
     &       ndim
        end subroutine

      end interface
