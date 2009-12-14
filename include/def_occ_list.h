      integer, parameter ::
     &     olist_vtx_inf_inc = 10,
     &     olist_occ_inc = 30
      

      type occ_list
      integer ::
     &     n_vtx_inf, n_occ, max_vtx_inf, max_occ
      integer, pointer ::
     &     vtx_inf(:,:),
     &     occ(:,:,:)
      end type occ_list
