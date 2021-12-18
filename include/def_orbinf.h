      type orbinf
        integer ::
     &     nsym, ngas, nspin
        integer ::
     &      ntoob,caborb,nbast,nxbast,nactel,nactorb,lsym,imult,ims,
     &      ncore_ext,mem_ext,bufflen,ibufflen
        character(len=256) ::
     &      name_intfile_ext
        integer, pointer ::
     &     igassh(:,:),
     &     nbas(:), ntoobs(:), ireots(:), ireost(:),
     &     igamorb(:), igasorb(:),
     &     mostnd(:,:,:),
     &     iad_gas(:), gas_reo(:),
     &     ihpvgas(:,:), ngas_hpv(:), nactt_hpv(:), norb_hpv(:,:),
     &     idx_gas(:), ioff_gas(:), cab_orb(:), nxbas(:),
     &     xreosym(:), xreotyp(:), ext_gamorb(:), ext_fcreo(:)
        integer ::
     &     n_bound_orbs, n_freeze_rcmd
        integer, pointer ::
     &     isym_bound_orbs(:)

      end type orbinf
