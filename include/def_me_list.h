*----------------------------------------------------------------------*
*     matrix-element lists (type me_list)
*----------------------------------------------------------------------*
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
     &     mxlen_melabel = len_opname*2
*----------------------------------------------------------------------*
*     me_list definition
*----------------------------------------------------------------------*
      type me_list
      
         ! my label
         character(mxlen_melabel) ::
     &        label

         ! the assigned operator
         type(operator), pointer ::
     &        op

         ! the file handle with the actual matrix elements
         ! (may be incore)
         type(filinf), pointer ::
     &        fhand

         ! assigned frequency
         real(8) ::
     &        frequency

*----------------------------------------------------------------------*
*     symmetry and storage info:
*----------------------------------------------------------------------*
        integer ::
     &       gamt,                    ! total symmetry
     &       s2,mst                   ! total spin and ms

        integer ::
     &       nsym                     ! for convenience: the nsym
                                      ! value used for dimensioning
                                      ! the arrays below (should be
                                      ! the same as orbinf%nsym)
                                      ! set in init_me_list()
                                      ! do not change elsewhere!

        integer ::
     &       absym,             ! symmetry on interchange of alpha/beta
     &                          ! (time reversal sym.) values: 0/+1/-1
     &       casym              ! symmetry on interchange C<->A
                                ! (adjungation) values: 0/+1/-1

        logical ::
     &       fix_vertex_ms      ! Flag as to whether the ms value for each
                                ! block should be fixed (usually false, but
                                ! set up to be true for the non-antisym
                                ! Hamiltonian elements needed for the F12 Z
                                ! intermediate). May need later sophistication.
*----------------------------------------------------------------------*
*       absym, casym imply that only triangles are saved
*       convention: MS >= 0, index(alpha) >= index(beta)
*                   index(A) >= index(C)
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
*     dimension info for me_list:
*----------------------------------------------------------------------*
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
        type(leninfx2), pointer ::
     &       ld_op_gmox(:)     ! leading dim. per IRREP, MS and occ. class,
        integer, pointer ::
     &       idx_graph(:,:,:)  ! graph needed for addressing of
                               ! the corresp. PHVX/CA part
      end type me_list
