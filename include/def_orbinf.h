      type orbinf
        integer ::
     &     nsym, ngas
        integer ::
     &     ntoob
        integer, allocatable ::
     &     igassh(:,:),
     &     ntoobs(:), ireots(:), ireost(:),
     &     igamorb(:), igasorb(:),
     &     mostnd(:,:,:),
     &     iad_gas(:),
     &     ihpvgas(:), ngas_hpv(:), nactt_hpv(:),
     &     idx_gas(:), ioff_gas(:)

      end type orbinf
