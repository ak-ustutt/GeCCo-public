*----------------------------------------------------------------------*
      subroutine init_reo_info(reo_info)
*----------------------------------------------------------------------*
*     initialize all pointers 
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'def_reorder_info.h'

      type(reorder_info), intent(inout) ::
     &     reo_info

      reo_info%nreo = 0
      reo_info%nreo_i0 = 0
      reo_info%n_op_reo = 0

      reo_info%nca_vtx => null()

      reo_info%reo => null()
      
      reo_info%merge_stp1 => null()
      reo_info%merge_stp2 => null()
      reo_info%iocc_reo  => null()
      reo_info%iocc_opreo0 => null()
      reo_info%from_to => null()

      reo_info%map_reo1c => null()
      reo_info%map_reo1a => null()
      reo_info%map_reo2c => null()
      reo_info%map_reo2a => null()

      reo_info%cinfo_opreo0c => null()
      reo_info%cinfo_opreo0a => null()
      reo_info%cinfo_reo_c => null()
      reo_info%cinfo_reo_a => null()

      return
      end
