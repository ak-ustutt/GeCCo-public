! need to include before: opdim.h      
      type cntr_arc

         integer ::
     &      link(2)   ! indices of operators in contraction list
         integer ::
     &      occ_cnt(ngastp,2)
                      ! occupation describing the contraction
                      ! convention: V1(C)V2(C^+)
      end type cntr_arc

      type cntr_vtx
 
        integer ::  ! type of operator defining vertex
     &     idx_op,  !  index
     &     iblk_op  !  block
!     &     ncntr    ! number of contractions

      end type cntr_vtx

      type contraction

        integer ::     ! type of result:
     &     idx_res,    ! index of operator type (0 for scalar) 
     &     iblk_res    ! block of operator type (0 for scalar)

        real(8) ::     
     &       fac       ! prefactor

        integer ::
     &       nvtx,     ! number of vertices (=operators)
     &       narc,     ! number of arcs (=interaction lines/raw contractions)
     &       nfac,     ! number of factors (factorization info)
     &       mxvtx,    ! current sizes of subarrays
     &       mxarc,    ! current sizes of subarrays
     &       mxfac     ! current sizes of subarrays

        type(cntr_vtx), pointer ::
     &       vertex(:) ! description of vertices
        type(cntr_arc), pointer ::
     &       arc(:)    ! description of arcs
        integer, pointer ::
     &       inffac(:,:) ! factorization info (4,nfac)
*----------------------------------------------------------------------*
*	factorization info organized as:
*         ((vertex1,vertex2,intermediate1),
*          (vertex3/intermediate1,vertex4,intermediate2/result))
*       vertex1,...        number of vertex (1...n)
*	intermediate1, ... further numbers (n+1,...)
*	result             0
*----------------------------------------------------------------------*
        
      end type contraction
