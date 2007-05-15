      type strmap_offsets
        integer, pointer ::
     &     msms(:), msmsgmgm(:)
      end type strmap_offsets

      type strmapinf
        integer, pointer ::
     &     idx_strmap(:)        ! for each pair of graphs def'd above:
                                ! start index of strmap on file
                                ! or -1 if none available (yet)
        integer ::
     &     idx_last
        type(filinf) ::
     &     ffstrmap
        type(strmap_offsets), pointer ::
     &       offsets(:)
      end type strmapinf
