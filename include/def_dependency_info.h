      type dependency_info
      
        integer ::
     &       ntargets, ndepend
        integer, pointer ::
     &       idxlist(:),
     &       depends_on_idxlist(:,:)

      end type dependency_info
