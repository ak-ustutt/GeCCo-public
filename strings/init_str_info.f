*----------------------------------------------------------------------*
      subroutine init_str_info(str_info)
*----------------------------------------------------------------------*
*     initialize the str_info array
*----------------------------------------------------------------------*

      implicit none

      include 'def_graph.h'
      include 'def_strinf.h'

      type(strinf), intent(inout) ::
     &     str_info
      
      str_info%ngraph = 0
      str_info%max_igtyp = 0
      str_info%max_idxms = 0
      str_info%gtab => null()
      str_info%ispc_typ => null()
      str_info%ispc_occ => null()
      str_info%igas_restr => null()
      str_info%g => null()

      return
      end
