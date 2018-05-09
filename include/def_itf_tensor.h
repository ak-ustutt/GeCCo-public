      integer, parameter ::
     &     index_len = 8

      type itf_tensor
      character(len=maxlen_bc_label) ::
     &     name                 ! Name of tensor
      character(len=index_len) ::
     &     idx,                 ! Index string, for now only rank 2 and 4
     &     next_idx,            ! Next result index string
     &     prev_idx             ! Previous result index string
      character(len=index_len) ::
     &     c_ops,               ! Array of creation operators
     &     a_ops                ! Array of annihilation operators
      integer ::
     &     rank,                ! Tensor rank
     &     idx_convention,      ! Index convention ie. integrals or amplitudes
     &     t_numb               ! t1=1, t2=2, res=3
      integer, dimension(4) ::
     &     ncops,               ! Number of h/p/v/x creation operators
     &     naops                ! Number of h/p/v/x annihilation operators
      real(8) ::
     &     fact                 ! Factor
      end type itf_tensor
