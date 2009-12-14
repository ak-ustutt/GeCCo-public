      type strmap_offsets
        integer, pointer ::
     &     msms(:), msmsgmgm(:)
      end type strmap_offsets

      type flpmap_offsets
        integer, pointer ::
     &     ms(:), msgm(:)
      end type flpmap_offsets

      type strmapinf
        integer ::
     &     mxgraph
        integer, pointer ::
     &     idx_strmap(:)        ! for each pair of graphs def'd above:
                                ! start index of strmap on file
                                ! or -1 if none available (yet)
        integer ::
     &     idx_last
        type(filinf) ::
     &     ffstrmap
        integer, pointer ::
     &     maxlen_blk(:)
        type(strmap_offsets), pointer ::
     &       offsets(:)

        ! info for frozen-core-maps
        integer, pointer ::
     &     idx_fcmap(:)
        integer, pointer ::
     &     maxlen_blk_fc(:)
        type(flpmap_offsets), pointer ::
     &     offsets_fc(:)

        ! info for flip-maps
        integer, pointer ::
     &     idx_flipmap(:)
        integer, pointer ::
     &     maxlen_blk_flip(:)
        type(flpmap_offsets), pointer ::
     &     offsets_flip(:)

      end type strmapinf
