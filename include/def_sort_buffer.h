      type sort_bin
      
      integer ::
     &     length
      real(8), pointer ::
     &     value(:)
      integer, pointer ::
     &     index(:)
      integer, pointer ::
     &     chain(:)

      end type sort_bin

      type sort_buffer
      integer ::
     &     nbins, max_length, max_chain
      type(sort_bin), pointer ::
     &     bin(:)
      end type sort_buffer
