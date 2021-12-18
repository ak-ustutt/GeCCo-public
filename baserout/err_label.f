      subroutine err_label(where,what_type,what_label)
      
      implicit none
      include 'stdunit.h'

      character(len=*), intent(in) ::
     &     where,what_type, what_label

      call quit(1,'err_label for:'//trim(where),
     &     trim(what_type)//
     &     ' label not on list: "'//
     &     trim(what_label)//'"')

      end
