      type optimize_status
        type(file_array), pointer ::
     &       ffrsbsp(:), ffvsbsp(:)
        integer ::
     &       ndim_rsbsp, ndim_vsbsp, mxdim_sbsp
        real(8) ::
     &       energy_last, xngrd_last, crate_last, trrad
        real(8), pointer ::
     &       sbspmat(:)
        integer, pointer ::
     &       iord_rsbsp(:), iord_vsbsp(:)
      end type optimize_status
