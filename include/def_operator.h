*----------------------------------------------------------------------*
*     parameters
*----------------------------------------------------------------------*
      ! operator types
      integer, parameter ::
     &     optyp_operator = 1,      ! single vertex operator
     &     optyp_density  = 2,      ! DX/EX 2-vertex operator
     &     optyp_intermediate = 3   ! general n-vertex operator

      ! vertex types as returned by vtx_type()
      integer, parameter ::
     &     vtxtyp_scalar= 0,   ! scalar
     &     vtxtyp_ph    = 1,   ! P/H spaces only (incl. X)
     &     vtxtyp_ph_ex = 2,   !  typ 1, only excitations 
     &     vtxtyp_ph_dx = 3,   !  typ 1, only deexcitations 
     &     vtxtyp_val   = 4    ! V spaces as well

      integer, parameter ::
     &     len_opname = 32
*----------------------------------------------------------------------*
*     operator definition
*----------------------------------------------------------------------*
      type operator
        integer ::
     &     id   ! might be obsolete
        character(len_opname) ::
     &     name
        character(2*len_opname+2) ::  
     &     assoc_list     ! name of the currently associated ME-list
*----------------------------------------------------------------------*
*
*       dagger == true: Operator defined as adjungate, 
*                       in particular: operator elements
*                       are saved in same sequence as 
*                       in the non-daggered case
*                       
*----------------------------------------------------------------------*
	logical ::
     &     dagger               ! daggered operator: OBSOLETE!!!
                                ! C <-> A are to be interchanged
        integer ::
     &     type,                ! 1: operator, 2: density, 3: intermed.
     &     njoined,             ! for intermediate only: number of joined 
     &                          !      vertices
     &     n_occ_cls,           ! number of occupation classes
     &     ngas,nspin           ! info from orb_info<-for convenience 

        integer ::
     &     hermitian            ! 1(-1): (Anti-)Hermitian, 0: not Hermitian

        integer ::
     &     order,               ! perturbation order
     &     species              ! 1: t-amplitude, 2: lagr. multipl., 3: other
        integer, pointer ::
     &     ifreq(:)             ! frequency index

        integer, pointer ::
     &     blk_version(:)       ! block version (in case of similar blocks)

        logical ::
     &       formal                   ! formal operator only?
        logical, pointer ::
     &       formal_blk(:)            ! is this block formal only?
        integer, pointer ::
     &       ihpvca_occ(:,:,:)        ! actual hole/particle/valence (HPV)
                                      ! occupations per creation/annihilation
                                      ! string (CA)
        integer, pointer ::
     &       ica_occ(:,:)             ! occupations summed up for CA
        integer, pointer :: 
     &       igasca_restr(:,:,:,:,:,:)    ! associated subspace restrictions
           ! igasca_restr(2,      ngas, 2  , 2         , nspin, n_occ_cls)
           !              min/max, spc, C/A, restr/mask, spcas, occ.class
           !              restr: defining allowed subspace occupations
           !              mask-restr: set of restrictions defining subspc.
           !                     occupations to be skipped
	                      ! hpv occupation + subspace restriction =
			      ! occupation class
      end type operator	      
