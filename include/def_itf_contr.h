      integer, parameter ::
     &     index_len = 8        ! Length of index string

      type spin_cases

      character(len=maxlen_bc_label) ::
     &     name    
      integer, dimension(4,6) ::
     &     cases                ! matrix containing possible spin cases each row is a different case
      integer ::
     &     ncase = 0

      end type spin_cases


      type itf_contr
      ! Object that holds information needed to define a line in an ITF algo file

      character(len=maxlen_bc_label) ::
     &     label_t1,            ! Name of first tensor in the contraction
     &     label_t2,            ! Name of second tensor in the contraction
     &     label_res            ! Name of tensor result
      character(len=index_len) ::
     &     idx1,                ! First tensor index string
     &     idx2,                ! Second tensor index string
     &     idx3                 ! Result tensor index string
      integer ::
     &     rank1,               ! Rank of first tensor
     &     rank2,               ! Rank of second tensor
     &     rank3,               ! Rank of result tensor
     &     logfile,             ! File to print to
     &     permute,             ! 0,1,2,3: permutation cases
     &     spin_idx=0,          ! Groups spin summed lines together
     &     command              ! Type of contraction, borrowed from formula_item
      logical ::
     &     inter(3)=.false.     ! True if result is intermediate
      real(8) ::
     &     fact                 ! Factor
      integer ::
     &     spin_case(4) ! Spin case of result, ie. [1111] = all alpha

      type(spin_cases), pointer ::
     &     inter_spins(:) => null()
      integer ::
     &     ninter = 0
      logical ::
     &     print_line = .true.

      end type itf_contr






      ! TODO: Get rid of below...

      type itf_spin_parts

      type(itf_contr), pointer ::
     &     spin_parts(:) => null()

      end type itf_spin_parts


      type itf_intermediate_spin
      ! Container of possible spin cases of intermediates used in an itf_contr

      type(itf_spin_parts), pointer ::
     &     spin_case(:) => null()               ! Spin summed binary contraction
      character(len=maxlen_bc_label) ::
     &     name                         ! Name of intermediate

      end type itf_intermediate_spin


      type itf_intermediate
      ! Container of intermediates used in an itf_contr, ie. I1, I2, I3...

      type(itf_intermediate_spin), pointer ::
     &     interm(:) => null()             ! Array of spin summed cases of a specific intermediate
      integer ::
     &     size

      end type itf_intermediate
