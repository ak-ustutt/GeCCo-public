      type orbinf
        integer ::
     &     nsym, ngas, nspin
        integer ::
     &      ntoob,caborb,nbast,nxbast,nactel,nactorb,lsym,imult,ims,
     &      ncore_mpro,mem_mpro
        character(len=256) ::
     &      name_intfile_mpro
        integer, pointer ::
     &     igassh(:,:),
     &     nbas(:), ntoobs(:), ireots(:), ireost(:),
     &     igamorb(:), igasorb(:),
     &     mostnd(:,:,:),
     &     iad_gas(:), gas_reo(:),
     &     ihpvgas(:,:), ngas_hpv(:), nactt_hpv(:), norb_hpv(:,:),
     &     idx_gas(:), ioff_gas(:), cab_orb(:), nxbas(:),
     &     xreosym(:), xreotyp(:)
        integer ::
     &     n_bound_orbs, n_freeze_rcmd
        integer, pointer ::
     &     isym_bound_orbs(:)

      end type orbinf
