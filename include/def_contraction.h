! need to include before: opdim.h
      integer, parameter ::
     &     ld_inffac = 5,
     &     max_vtx_group = 10

      type cntr_arc

         integer ::
     &      link(2)   ! indices of operators in contraction list
         integer ::
     &      occ_cnt(ngastp,2)
                      ! occupation describing the contraction
                      ! convention: V1(C)V2(C^+)
      end type cntr_arc

c      type cntr_proto_arc
c
c         integer ::
c     &      link_group(max_vtx_group,2)
c                       ! list of verticex groups to be contracted
c         integer ::
c     &      occ_cnt_min(ngastp,2),  minimum
c     &      occ_cnt_max(ngastp,2)
c      end type cntr_arc
c
      type cntr_vtx

        integer ::  ! type of operator defining vertex
     &     idx_op,  !  index
     &     iblk_op  !  block (super vertices: compound vertex/block index)
        logical ::
     &     dagger   ! the operators enters as its adjoint

      end type cntr_vtx

      ! a type for storing explicit index information
      type string_element
        integer :: vtx     ! vertex to which current index belongs
        integer :: ca      ! creation or annihilation
        integer :: hpvx    ! orbital space
        integer :: cnt     ! contraction (ext=F) or external vertex (ext=T) for index
        integer :: idx     ! index number (within orbital space
        logical :: ext     ! ext=F: cnt is contraction number; ext=T: cnt is external vertex
        logical :: del     ! marker for internal purposes
      end type string_element

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
     &       joined(:,:), !  vertices per supervertex (0:nvtx,nsupvtx)
     &                    !  0: number of vertices, 1-nvtx: vertex numbers
     &       svertex(:)   !  supervertex to which vertex belongs (nvtx)
        integer, pointer ::
     &       inffac(:,:)  ! factorization info (4,nfac)
        logical ::
     &       unique_set   ! unique representation (topo etc.) defined?
        integer(8), pointer ::
     &       vtx(:), topo(:,:), xlines(:,:)
        integer(8) ::
     &       hash         ! for quick comparison
        logical ::
     &       index_info         ! index information defined? (see below)
        integer ::
     &       nidx, nxidx
        type(string_element), pointer ::
     &       contr_string(:),   ! index information for contraction (= diagram labels)
     &       result_string(:)   ! index information for final result (= pairing for open lines)
        integer ::
     &       total_sign         ! sign for explicitly labelled diagram
        real(8) ::
     &       eqvl_fact          ! equivalent line factor for external post processing
*----------------------------------------------------------------------*
*	factorization info organized as:
*         ((vertex1,vertex2,intermediate1),
*          (vertex3/intermediate1,vertex4,intermediate2/result))
*       vertex1,...        number of vertex (1...n)
*	intermediate1, ... further numbers (n+1,...)
*	result             0
*----------------------------------------------------------------------*
      end type contraction

      integer, parameter ::
     &     maxlen_bc_label = 32

      type binary_contr
      integer ::
     &     n_operands, n_cnt
      real(8) ::
     &     fact, fact_itf
!                fact_itf is corrected for different merge convention
!                (should give the "standard" sign from diagram rules)
      character(len=maxlen_bc_label) ::
     &     label_res, label_op1, label_op2
      integer ::
     &     iblk_res, iblk_op1, iblk_op2,
     &     nj_res, nj_op1, nj_op2,
     &     ngas, nspin
      logical ::
     &     tra_res, tra_op1, tra_op2,
     &     perm(ngastp,2)
      integer, pointer ::
     &     occ_res(:,:,:),
     &     occ_op1(:,:,:),
     &     occ_op2(:,:,:),
     &     rst_res(:,:,:,:,:,:),
     &     rst_op1(:,:,:,:,:,:),
     &     rst_op2(:,:,:,:,:,:),
     &     occ_ex1(:,:,:),
     &     occ_ex2(:,:,:),
     &     occ_cnt(:,:,:),
     &     rst_ex1(:,:,:,:,:,:),
     &     rst_ex2(:,:,:,:,:,:),
     &     rst_cnt(:,:,:,:,:,:),
     &     merge_op1(:), merge_op2(:),
     &     merge_op1op2(:), merge_op2op1(:),
     &     itf_index_info(:)
      end type binary_contr

      type reorder
      character(len=maxlen_bc_label) ::
     &     label_out, label_in
      integer ::
     &     iblk_out, iblk_in,
     &     nj_out, nj_in
      logical ::
     &     tra_out, tra_in
      integer ::
     &     nreo,nreo_i0,sign,
     &     ngas, nspin
      integer, pointer ::
     &     from_to(:,:),
     &     occ_shift(:,:,:),
     &     occ_opout(:,:,:),
     &     rst_opout(:,:,:,:,:,:),
     &     occ_opin(:,:,:),
     &     rst_opin(:,:,:,:,:,:),
     &     occ_op0(:,:,:),
     &     merge_stp1(:),merge_stp1inv(:),
     &     merge_stp2(:),merge_stp2inv(:)
      end type reorder
