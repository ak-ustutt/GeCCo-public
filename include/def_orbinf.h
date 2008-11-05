      type orbinf
        integer ::
     &     nsym, ngas, nspin
        integer ::
     &      ntoob,caborb,nbast,nxbast
        integer, allocatable ::
     &     igassh(:,:),
     &     nbas(:), ntoobs(:), ireots(:), ireost(:),
     &     igamorb(:), igasorb(:),
     &     mostnd(:,:,:),
     &     iad_gas(:), gas_reo(:),
     &     ihpvgas(:,:), ngas_hpv(:), nactt_hpv(:), norb_hpv(:,:),
     &     idx_gas(:), ioff_gas(:), cab_orb(:), nxbas(:),
     &     xreosym(:), xreotyp(:)

      end type orbinf
