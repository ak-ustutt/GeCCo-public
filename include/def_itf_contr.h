      integer, parameter ::
     &     INDEX_LEN = 8,       ! Length of index string
     &     MAXINT = 8,          ! Maximum number of intermediates that contribute to a result
     &     MAX_SPIN_CASES = 10  ! Maximum number of intermediates that contribute to a result

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
     &     label_res            ! Name of tensor result
      character(len=INDEX_LEN) ::
     &     idx1,                ! First tensor index string
     &     idx2,                ! Second tensor index string
     &     idx3,                ! Result tensor index string
     &     inter1,              ! Intermediate spin name (ie. abab or aaaa)
     &     inter2               ! Intermediate spin name (ie. abab or aaaa)
      integer ::
     &     rank1,               ! Rank of first tensor
     &     rank2,               ! Rank of second tensor
     &     rank3,               ! Rank of result tensor
     &     logfile,             ! File to print to
     &     permute,             ! 0,1,2,3: permutation cases
     &     command,             ! Type of contraction, borrowed from formula_item
     &     contri,              ! Number of contraction indices in a line
     &     nops1(ngastp),
     &     nops2(ngastp),
     &     nops3(ngastp),
     &     c(4,2),              ! Operator numbers of contraction index
     &     e1(4,2),             ! Operator numbers of external index 1
     &     e2(4,2),             ! Operator numbers of external index 2
     &     spin_cases           ! Number of printed spin cases
      logical ::
     &     inter(3) = .false.,  ! True if tensor is an intermediate
     &     int(3) = .false.,    ! True if tensor is an integral
     &     den(3) = .false.,    ! True if tensor is a denisty matrix
     &     symm_res = .false.,   ! True if intermediate contributes to a symmetric residual
     &     binary = .true.,     ! True if a binary contraction
     &     permutation = .false., ! True if the line is a result of the (1+P_ij^ab) permutation
     &     product = .false., ! True if the line is a tensor product
     &     perm_case(4),
     &     j_int = .false.      ! True if integral is special
      real(8) ::
     &     fact                 ! Factor

      ! Objects needed in the search for intermediates
      type(spin_cases), pointer ::
     &     inter_spins(:) => null()     ! Array of intermediates + spin cases
      integer ::
     &     ninter = 0,                  ! Number of intermediates used in a contraction
     &     itype(INDEX_LEN)
      logical ::
     &     print_line = .true.          ! Should this line be printed
      type(spin_info) ::
     &     t_spin(3),                   ! Spin info of result, t1 and t2 tensors
     &     i_spin                       ! Spin info of an intermediate

      end type itf_contr


*----------------------------------------------------------------------*
      type itf_contr2
*----------------------------------------------------------------------*
!     Object that holds information needed to define a line in an ITF algo file
!            Res[abij] += T1[acik] T2[cbkj]
!     Label: 3            1        2
*----------------------------------------------------------------------*

      character(len=MAXLEN_BC_LABEl) ::
     &     label_t1,            ! Name of first tensor in the contraction
     &     label_t2,            ! Name of second tensor in the contraction
     &     label_res            ! Name of tensor result
      character(len=INDEX_LEN) ::
     &     idx1,                ! First tensor index string
     &     idx2,                ! Second tensor index string
     &     idx3,                ! Result tensor index string
     &     inter1,              ! Intermediate spin name (ie. abab or aaaa)
     &     inter2               ! Intermediate spin name (ie. abab or aaaa)
      integer ::
     &     rank1,               ! Rank of first tensor
     &     rank2,               ! Rank of second tensor
     &     rank3,               ! Rank of result tensor
     &     logfile,             ! File to print to
     &     permute,             ! 0,1,2,3: permutation cases
     &     command,             ! Type of contraction, borrowed from formula_item
     &     contri,              ! Number of contraction indices in a line
     &     nops1(ngastp),
     &     nops2(ngastp),
     &     nops3(ngastp),
     &     c(4,2),              ! Operator numbers of contraction index
     &     e1(4,2),             ! Operator numbers of external index 1
     &     e2(4,2),             ! Operator numbers of external index 2
     &     spin_cases           ! Number of printed spin cases
      logical ::
     &     inter(3) = .false.,  ! True if tensor is an intermediate
     &     int(3) = .false.,    ! True if tensor is an integral
     &     den(3) = .false.,    ! True if tensor is a denisty matrix
     &     symm_res = .false.,   ! True if intermediate contributes to a symmetric residual
     &     binary = .true.,     ! True if a binary contraction
     &     permutation = .false., ! True if the line is a result of the (1+P_ij^ab) permutation
     &     product = .false., ! True if the line is a tensor product
     &     perm_case(4),
     &     j_int = .false.      ! True if integral is special
      real(8) ::
     &     fact                 ! Factor

      ! Objects needed in the search for intermediates
      type(spin_cases), pointer ::
     &     inter_spins(:) => null()     ! Array of intermediates + spin cases
      integer ::
     &     ninter = 0,                  ! Number of intermediates used in a contraction
     &     itype(INDEX_LEN)
      logical ::
     &     print_line = .true.          ! Should this line be printed
      type(spin_info) ::
     &     t_spin(3),                   ! Spin info of result, t1 and t2 tensors
     &     i_spin                       ! Spin info of an intermediate

      integer ::
     &      nspin_cases
      type(spin_cases2)
     &     all_spins(10)

      end type itf_contr2


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
