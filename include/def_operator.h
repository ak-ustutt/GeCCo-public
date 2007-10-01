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
*----------------------------------------------------------------------*
*     auxiliary types
*----------------------------------------------------------------------*
      type leninf
        integer, pointer ::
     &     gam_ms(:,:)
      end type leninf
      type leninfx
        integer ::
     &     maxd
        integer, pointer ::
     &     ndis(:,:)
        integer, pointer ::
     &     did(:,:,:)
        integer, pointer ::
     &     d_gam_ms(:,:,:)
      end type leninfx
      type leninfx2
        integer, pointer ::
     &     d_gam_ms(:,:,:)
      end type leninfx2

      integer, parameter ::
     &     len_opname = 8
*----------------------------------------------------------------------*
*     operator definition
*----------------------------------------------------------------------*
      type operator
        integer ::
     &     id
        character ::
     &     name*(len_opname)
*----------------------------------------------------------------------*
*     basic definitions:
*
*       dagger == true: Operator defined as adjungate, 
*                       in particular: operator elements
*                       are saved in same sequence as 
*                       in the non-daggered case
*                       
*       absym, casym imply that only triangles are saved
*       convention: MS >= 0, index(alpha) >= index(beta)
*                   index(A) >= index(C)
*----------------------------------------------------------------------*
	logical ::
     &     dagger               ! daggered operator:
                                ! C <-> A are to be interchanged
        integer ::
     &     type,                ! 1: operator, 2: density, 3: intermed.
     &     njoined              ! for intermediate only: number of joined 
                                !      vertices
        integer ::
     &     absym,               ! symmetry on interchange of alpha/beta
     &                          ! (time reversal sym.) values: 0/+1/-1
     &     casym                ! symmetry on interchange C<->A
                                ! (adjungation) values: 0/+1/-1
        
        integer ::
     &       gamt,                    ! total symmetry
     &       s2,mst,                  ! total spin and ms
     &       n_occ_cls                ! number of occupation classes

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
	! dimension info for operators
	integer :: 
     &       len_op	      ! total length
	integer, pointer ::
     &       off_op_occ(:),   ! offset for occupation class
     &       len_op_occ(:)    ! length per occupation class
        type(leninf), pointer ::
     &       off_op_gmo(:),    ! offset for IRREP, MS and occupation class
     &       len_op_gmo(:)     ! length per IRREP, MS and occupation class
        type(leninfx), pointer ::
     &       off_op_gmox(:)    ! offset per IRREP, MS and occupation class,
                               ! and IRREP and MS distr. over HPV per C/A
        type(leninfx2), pointer ::
     &       len_op_gmox(:)    ! length per IRREP, MS and occupation class,
        integer, pointer ::
     &       idx_graph(:,:,:)     ! graph needed for addressing of
                              ! the corresp. PHV/CA part
      end type operator	      
