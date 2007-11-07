      type orbinf
        integer ::
     &     nsym, ngas
        integer ::
     &      ntoob,caborb,nbast
        integer, allocatable ::
     &     igassh(:,:),
     &     nbas(:), ntoobs(:), ireots(:), ireost(:),
     &     igamorb(:), igasorb(:),
     &     mostnd(:,:,:),
     &     iad_gas(:), gas_reo(:),
     &     ihpvgas(:), ngas_hpv(:), nactt_hpv(:),
     &     idx_gas(:), ioff_gas(:), cab_orb(:),
     &     xreosym(:), xreotyp(:)

      end type orbinf
