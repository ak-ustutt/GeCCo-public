      subroutine clean_sort_buffer(sbuffer)

      implicit none

      include 'def_sort_buffer.h'
      include 'ifc_memman.h'

      type(sort_buffer), intent(inout) ::
     &     sbuffer
      integer ::
     &     ifree

      ! will deallocate the complete memory section associated with 
      ! the buffer:
      ifree = mem_flushmark('sbuffer')

      ! also deallocate the pointer array
      deallocate(sbuffer%bin)

      return
      end
