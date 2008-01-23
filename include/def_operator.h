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
     &     len_opname = 8
*----------------------------------------------------------------------*
*     operator definition
*----------------------------------------------------------------------*
      type operator
        integer ::
     &     id   ! might be obsolete
        character(len_opname) ::
     &     name
        character(2*len_opname) ::  
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
     &     dagger               ! daggered operator:
                                ! C <-> A are to be interchanged
        integer ::
     &     type,                ! 1: operator, 2: density, 3: intermed.
     &     njoined,             ! for intermediate only: number of joined 
     &                          !      vertices
     &     n_occ_cls,           ! number of occupation classes
     &     ngas                 ! info from orb_info<-for convenience 

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
     &       igasca_restr(:,:,:,:,:)    ! associated subspace restrictions
           ! igasca_restr(2,      ngas, 2  , 2               , n_occ_cls)
           !              min/max, spc, C/A, restr/mask-restr, occ.class
           !              restr: defining allowed subspace occupations
           !              mask-restr: set of restrictions defining subspc.
           !                     occupations to be skipped
	                      ! hpv occupation + subspace restriction =
			      ! occupation class
      end type operator	      
