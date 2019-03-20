      integer, parameter ::
     &     index_len = 8,       ! Length of index string
     &     maxint = 4           ! Maximum number of intermediates that contribute to a result

*----------------------------------------------------------------------*
      type spin_cases
*----------------------------------------------------------------------*
!     Hold information about a tensor and the spin summed cases which need to be
!     printed
*----------------------------------------------------------------------*

      character(len=maxlen_bc_label) ::
     &     name                 ! Name of tensor
      integer, dimension(4,6) ::
     &     cases = 0                ! Matrix containing possible spin cases each row is a different case
      integer ::
     &     ncase = 0            ! Number of different spin cases

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
      type itf_contr
*----------------------------------------------------------------------*
!     Object that holds information needed to define a line in an ITF algo file
*----------------------------------------------------------------------*

      character(len=maxlen_bc_label) ::
     &     label_t1,            ! Name of first tensor in the contraction
     &     label_t2,            ! Name of second tensor in the contraction
     &     label_res            ! Name of tensor result
      character(len=index_len) ::
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
     &     contri               ! Number of contraction indices in a line
      logical ::
     &     inter(3) = .false.,  ! True if tensor is an intermediate
     &     int(3) = .false.,    ! True if tensor is an integral
     &     symm = .false.       ! True if going to symmetrise result
      real(8) ::
     &     fact                 ! Factor
      integer ::
     &     spin_case(4) ! Spin case of result, ie. [1111] = all alpha

      ! Objects needed in the search for intermediates
      type(spin_cases), pointer ::
     &     inter_spins(:) => null()     ! Array of intermediates + spin cases
      integer ::
     &     ninter = 0                   ! Number of intermediates used in a contraction
      logical ::
     &     print_line = .true.          ! Should this line be printed
      type(spin_info) ::
     &     t_spin(3)

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
     &     ops(2)               ! Information on which tensor the index belogns

      logical ::
     &     linked = .false.     ! Doesn't need a linking index

      end type pair
