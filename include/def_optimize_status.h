      type optimize_status
        type(file_array), pointer ::
     &       ffrsbsp(:), ffvsbsp(:), ffssbsp(:)
        integer ::
     &       ndim_rsbsp, ndim_vsbsp, ndim_ssbsp,
     &       mxdim_sbsp, nadd, ndel
        real(8) ::
     &       energy_last, xngrd_last, crate_last, trrad
        real(8), pointer ::
     &       sbspmat(:)
        integer, pointer ::
     &       iord_rsbsp(:), iord_vsbsp(:), iord_ssbsp(:)
      end type optimize_status
