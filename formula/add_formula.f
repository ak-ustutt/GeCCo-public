*----------------------------------------------------------------------*
      subroutine add_formula(form_info,label)
*----------------------------------------------------------------------*
*     add new entry to formula info
*----------------------------------------------------------------------*
      implicit none

      include 'def_filinf.h'
      include 'mdef_formula_info.h'

      type(formula_info), intent(inout) ::
     &     form_info
      character ::
     &     label*(*)

      character(form_maxlen_label*2) ::
     &     name

      type(formula_list), pointer ::
     &     list_pnt

      ! advance to end of formula list
      list_pnt => form_info%form_list
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%form)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt%next%next => null()
        list_pnt => list_pnt%next
      else if (form_info%nform.gt.0) then
        ! should only happen the first time:
        call quit(1,'add_formula','Suspicious unused entry !')
      end if
      allocate (list_pnt%form)
      
      ! assign label such that we can search for it
      list_pnt%form%label = trim(label)
c      list_pnt%form%fhand => null()
      write(name,'(a,".fml")') trim(list_pnt%form%label)
      call file_init(list_pnt%form%fhand,name,ftyp_sq_unf,0)

      form_info%nform = form_info%nform+1

      call update_form_arr(form_info)

      end
