      type optimize_status
        type(file_array), pointer ::
     &       ffrsbsp(:), ffvsbsp(:), ffssbsp(:), ffscr(:), ffext(:)
        integer ::
     &       ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &       mxdim_sbsp, nadd, ndel, mode_2nd
        real(8) ::
     &       energy_last, xngrd_last, crate_last, trrad, gamma
        real(8), pointer ::
     &       sbspmat(:)
        integer, pointer ::
     &       iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:)
        type(optimize_status), pointer ::
     &       next_state
        integer ::
     &       i_state ! Not realy necessary
      end type optimize_status
