! need to include before: opdim.h      
      integer, parameter ::
     &      ld_inffac = 5

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
     &     iblk_op  !  block (super vertices: compound vertex/block index)
        logical ::
     &     dagger   ! the operators enters as its adjoint

      end type cntr_vtx

      type contraction

        integer ::     ! type of result:
     &     idx_res,    ! index of operator type (0 for scalar) 
     &     iblk_res    ! block of operator type (0 for scalar)
        logical ::
     &     dagger      ! the result must be transposed 
                       ! (intended for formal purposes)

        real(8) ::     
     &       fac       ! prefactor

        integer ::
     &       nsupvtx,  ! number of super-vertices (joined vertices)
     &       nvtx,     ! number of vertices (=operators)
     &       narc,     ! number of arcs (=interaction lines/raw contractions)
     &       nxarc,    ! number of ext. arcs (=open lines[*])
     &       nfac,     ! number of factors (factorization info)
     &       mxvtx,    ! current sizes of subarrays
     &       mxarc,    ! current sizes of subarrays
     &       mxxarc,   ! current sizes of subarrays
     &       mxfac     ! current sizes of subarrays
        ! [*] only needed for multi-vertex results

        type(cntr_vtx), pointer ::
     &       vertex(:) ! description of vertices
        type(cntr_arc), pointer ::
     &       arc(:)    ! description of arcs
        type(cntr_arc), pointer ::
     &       xarc(:)   ! description of external arcs (optional [*])
        integer, pointer :: !  super-vertex info:
     &       joined(:,:), !  vertices per supervertex (nvtx,nsupvtx)
     &       svertex(:)   !  supervertex to which vertex belongs (nvtx)
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
