      integer, parameter ::
     &     INDEX_LEN = 8,       ! Length of index string
     &     MAXINT = 8,          ! Maximum number of intermediates that contribute to a result
     &     MAX_SPIN_CASES = 10,  ! Maximum number of intermediates that contribute to a result
     &     MAXVTX = 10  ! Maximum number of intermediates that contribute to a result

*----------------------------------------------------------------------*
      type spin_cases
*----------------------------------------------------------------------*
!     Hold information about a tensor and the spin summed cases which need to be
!     printed
*----------------------------------------------------------------------*

      character(len=MAXLEN_BC_LABEL) ::
     &     name                 ! Name of tensor
      integer, dimension(INDEX_LEN,MAX_SPIN_CASES) ::
     &     cases = 0                ! Matrix containing possible spin cases each row is a different case
      integer ::
     &     ncase = 0,           ! Number of different spin cases
     &     itype(INDEX_LEN)
      logical ::
     &     symm_res = .false.    ! True if intermediate contributes to a symmetric residual

      end type spin_cases



*----------------------------------------------------------------------*
      type twodarray
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      integer ::
     &     elements(2)

      end type twodarray


*----------------------------------------------------------------------*
      type spin_info
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      integer, pointer ::
     &     spin(:,:) => null()

      end type spin_info


*----------------------------------------------------------------------*
      type spin_info2
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      integer ::
     &     spin(2,INDEX_LEN/2)

      end type spin_info2


*----------------------------------------------------------------------*
      type spin_cases2
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*
      type(spin_info2) ::
     &      t_spin(3)

      end type spin_cases2


*----------------------------------------------------------------------*
      type itf_contr
*----------------------------------------------------------------------*
!     Object that holds information needed to define a line in an ITF algo file
!            Res[abij] += T1[acik] T2[cbkj]
!     Label: 3            1        2
*----------------------------------------------------------------------*

      character(len=MAXLEN_BC_LABEl) ::
     &     label_t1,            ! Name of first tensor in the contraction
     &     label_t2,            ! Name of second tensor in the contraction
     &     label_res,           ! Name of tensor result
     &     old_name             ! If symmetrising result, the previous label it had
      character(len=INDEX_LEN) ::
     &     idx1,                ! First tensor index string
     &     idx2,                ! Second tensor index string
     &     idx3,                 ! Result tensor index string
     &     vidx1(MAXVTX),                ! First tensor index string
     &     vidx2(MAXVTX),                ! Second tensor index string
     &     vidx3(MAXVTX)                 ! Result tensor index string
      integer ::
     &     rank1,               ! Rank of first tensor
     &     rank2,               ! Rank of second tensor
     &     rank3,               ! Rank of result tensor
     &     out,                 ! File to print to
     &     tout,                ! ITF Task file
     &     permute,             ! 0,1,2,3: permutation cases
     &     command,             ! Type of contraction, borrowed from formula_item
     &     contri,              ! Number of contraction indices in a line
     &     nops1(ngastp),
     &     nops2(ngastp),
     &     nops3(ngastp),
     &     c(4,2),              ! Operator numbers of contraction index
     &     e1(4,2),             ! Operator numbers of external index 1
     &     e2(4,2),             ! Operator numbers of external index 2
     &     spin_cases,          ! Number of printed spin cases
     &     cntr(4),             ! Keep count of various quantaties inbetween calls to command_to_itf()
     &     nj_op1,
     &     nj_op2,
     &     nj_res
      logical ::
     &     inter(3) = .false.,  ! True if tensor is an intermediate
     &     int(3) = .false.,    ! True if tensor is an integral
     &     den(3) = .false.,    ! True if tensor is a denisty matrix
     &     symm_res = .false.,   ! True if intermediate contributes to a symmetric residual
     &     binary = .true.,     ! True if a binary contraction
     &     product = .false., ! True if the line is a tensor product
     &     perm_case(4),
     &     j_int = .false.,     ! True if integral is special
     &     symmetric = .false.,
     &     nosym = .false.,     ! True if R[apiq] (residual has no symmetry between indicies)
     &     abba_line = .false., ! True if this line is the R[apiq] abba spin case
     &     k4e_line = .false., ! True if a K4E line
     &     intpp = .false.,
     &     tasks = .false.
      real(8) ::
     &     fact                 ! Factor

      integer ::
     &     ninter = 0,                  ! Number of intermediates used in a contraction
     &     itype(MAXINT,INDEX_LEN),
     &     itype2(MAXVTX,MAXINT,INDEX_LEN)
      type(spin_info) ::
     &     t_spin(3)                    ! Spin info of result, t1 and t2 tensors
      integer, pointer ::
     &  vertex(:) => null(),
     &  v1(:,:,:) => null(),
     &  v2(:,:,:) => null(),
     &  v3(:,:,:) => null(),
     &  vc1(:,:,:) => null(),
     &  vc2(:,:,:) => null(),
     &  vc3(:,:,:) => null(),
     &  vnops1(:) => null(),
     &  vnops2(:) => null(),
     &  vnops3(:) => null()

      integer ::
     &      nspin_cases
      type(spin_cases2)
     &     all_spins(20)

      end type itf_contr


*----------------------------------------------------------------------*
      type pair_list
*----------------------------------------------------------------------*
!     List of index paris
*----------------------------------------------------------------------*

      type(pair), pointer ::
     &     plist(:) => null()

      end type pair_list


*----------------------------------------------------------------------*
      type pair
*----------------------------------------------------------------------*
!     Paired index
*----------------------------------------------------------------------*

      character(len=1) ::
     &     pindex(2),           ! Holds index pair, creation op in (1), annihilation in (2)
     &     link                 ! If external indices appear on different tensors, need a contraction index to link them

      integer ::
     &     ops(2),              ! Information on which tensor the index belogns
     &     nval(3)              ! Numerical value of index (1=e, 2=a, 3=c, 4=x)

      logical ::
     &     linked = .false.     ! Doesn't need a linking index

      end type pair


*----------------------------------------------------------------------*
      type index_str
*----------------------------------------------------------------------*
!
*----------------------------------------------------------------------*

      character(len=1), pointer ::
     &     str(:) => null()

      integer, pointer ::
     &     cnt_poss(:) => null(),
     &     itype(:) => null()

      end type index_str
