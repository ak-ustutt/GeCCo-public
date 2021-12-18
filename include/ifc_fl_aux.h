      interface
        function find_fl_item(fl_start,command,command_list,nlist,
     &                        label,label_res,backward)
        import
        implicit none
        type(formula_item), pointer :: find_fl_item
        type(formula_item), pointer :: fl_start
        integer, intent(in), optional :: command
        integer, intent(in), optional :: command_list(*), nlist
        character(len=len_opname), intent(in), optional :: label
        character(len=len_opname), intent(in), optional :: label_res
        logical, intent(in), optional :: backward
        end function
      end interface
