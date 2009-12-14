      subroutine clean_occ_list(olist)

      implicit none

      include 'opdim.h'
      include 'def_occ_list.h'

      type(occ_list), intent(inout) ::
     &     olist
      
      olist%n_vtx_inf = 0
      olist%n_occ = 0
      
      deallocate(olist%vtx_inf)
      deallocate(olist%occ)

      olist%vtx_inf  => null()
      olist%occ => null()

      olist%max_vtx_inf = 0
      olist%max_occ = 0
      
      return
      end
