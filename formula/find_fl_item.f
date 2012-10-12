      function find_fl_item(fl_start,command,command_list,nlist,
     &                      label,label_res,backward)

      implicit none

      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_operator.h'
      include 'def_formula_item.h'

      type(formula_item), pointer :: find_fl_item
      type(formula_item), pointer :: fl_start
      integer, intent(in), optional :: command
      integer, intent(in), optional :: command_list(*), nlist
      character(len=len_opname), intent(in), optional :: label
      character(len=len_opname), intent(in), optional :: label_res
      logical, intent(in), optional :: backward

      logical ::
     &    back, ok
      integer, pointer ::
     &    comm(:)
      integer ::
     &    icom, jcom, ncom
      type(formula_item), pointer ::
     &       fl_pnt

      back = .false.
      if (present(backward)) back = backward
      ncom = 1
      if (present(nlist)) ncom = nlist
      allocate(comm(ncom))
      comm(1) = command_end_of_formula
      if (present(command)) comm(1) = command
      if (present(command_list)) comm(1:ncom) = command_list(1:ncom)

      fl_pnt => fl_start

      do
        if (.not.back) fl_pnt => fl_pnt%next
        if (     back) fl_pnt => fl_pnt%prev
        if (.not.associated(fl_pnt)) exit
       
        ok = .false.
        do icom = 1, ncom
          if (fl_pnt%command.eq.comm(icom)) jcom = icom
          ok = ok.or.(fl_pnt%command.eq.comm(icom))
        end do
        if (.not.ok) cycle

        if (present(label)) then
          ok = .false.
          select case (comm(jcom))
          case (command_bc,command_add_bc,
     &          command_bc_reo,command_add_bc_reo)
            ok = (trim(fl_pnt%bcontr%label_op1).eq.trim(label).or.
     &            trim(fl_pnt%bcontr%label_op2).eq.trim(label).or.
     &            trim(fl_pnt%bcontr%label_res).eq.trim(label))
          case (command_add_intm,command_cp_intm)
            ok = (trim(fl_pnt%bcontr%label_op1).eq.trim(label).or.
     &            trim(fl_pnt%bcontr%label_res).eq.trim(label))
          case (command_new_intermediate)
            ok = (trim(fl_pnt%interm%name).eq.trim(label))
          case (command_del_intermediate)
            ok = (trim(fl_pnt%label).eq.trim(label))
          case (command_reorder) 
            ok = (trim(fl_pnt%reo%label_in ).eq.trim(label).or.
     &            trim(fl_pnt%reo%label_out).eq.trim(label))
          end select
          if (ok) exit
        else if (present(label_res)) then
          ok = .false.
          select case (comm(jcom))
          case (command_bc,command_add_bc,
     &          command_bc_reo,command_add_bc_reo)
            ok = trim(fl_pnt%bcontr%label_res).eq.trim(label_res)
          case (command_add_intm,command_cp_intm)
            ok = trim(fl_pnt%bcontr%label_res).eq.trim(label_res)
          case (command_reorder) 
            ok = trim(fl_pnt%reo%label_out).eq.trim(label_res)
          end select
          if (ok) exit
        else
          exit
        end if

      end do    

      find_fl_item => fl_pnt

      deallocate(comm)

      return
      end function

