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
