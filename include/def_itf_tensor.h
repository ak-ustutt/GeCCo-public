      integer, parameter ::
     &     index_len = 8

      type itf_tensor
      character(len=maxlen_bc_label) ::
     &     name                 ! Name of tensor
      character(len=index_len) ::
     &     idx                  ! Index string, for now only rank 4
      character(len=index_len) ::
     &     c_ops,               ! Array of creation operators
     &     a_ops                ! Array of annihilation operators
      integer ::
     &     rank,                ! Tensor rank
     &     spin_case=0          ! Spin case (see below)
      real(8) ::
     &     fact                 ! Factor
! Index convention, integer....
      end type itf_tensor

! Spin case:
!  1 = 4T - T^+
