      subroutine init_occ_list(olist)

      implicit none

      include 'opdim.h'
      include 'def_occ_list.h'

      type(occ_list), intent(inout) ::
     &     olist
      
      olist%n_vtx_inf = 0
      olist%n_occ = 0
      
      allocate(olist%vtx_inf(4,olist_vtx_inf_inc))
      allocate(olist%occ(ngastp,2,olist_occ_inc))
      olist%max_vtx_inf = olist_vtx_inf_inc
      olist%max_occ = olist_occ_inc
      
      return
      end
