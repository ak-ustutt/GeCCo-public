      integer, parameter ::
     &     index_len = 8        ! Length of index string

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
      end type itf_contr


      type itf_intermediate_spin
      ! Link list of possible spin cases of intermediates used in an itf_contr

      type(itf_contr) ::
     &     inter_contr                  ! Spin summed binary contraction
      integer ::
     &     spin_case(index_len)         ! Spin case of result, ie. [1111] = all alpha

      type(itf_intermediate_spin), pointer ::
     &     next_spin_inter => null()    ! Pointer to the next spin case

      end type itf_intermediate_spin


      type itf_intermediate
      ! Link list of intermediates used in an itf_contr, ie. I1, I2, I3...

      type(itf_intermediate_spin) ::
     &     inter_contr_spin             ! Linked list of spin summed cases of a specific intermediate

      type(itf_intermediate), pointer ::
     &     next_inter => null()         ! Pointer to next intermediate

      end type itf_intermediate
