! include 'opdim.h'

      type reorder_list

      logical ::
     &     is_bc_result
      integer ::
     &     idxsuper, to, from
      integer ::
     &     occ_shift(ngastp,2)

      end type reorder_list

      type reorder_info
      
      ! raw info
      integer ::
     &     nreo
      type(reorder_list), pointer ::
     &     reo(:)
      ! processed info
      integer ::
     &     n_op_reo
      integer, pointer ::
     &     merge_stp1(:),
     &     merge_stp2(:),
     &     iocc_reo(:,:,:),
     &     iocc_opreo0(:,:,:)
      integer ::
     &     ncblk_reo0, nablk_reo0, ncblk_reo, nablk_reo
      integer, pointer ::
     &     map_reo1c(:), map_reo1a(:),
     &     map_reo2c(:), map_reo2a(:),
     &     cinfo_opreo0c(:,:), cinfo_opreo0a(:,:),
     &     cinfo_reo_c(:,:), cinfo_reo_a(:,:)

      end type reorder_info
