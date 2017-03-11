*----------------------------------------------------------------------*
*     parameters
*----------------------------------------------------------------------*
      integer, parameter ::
     &     maxpop = 6, maxcmp = 99, maximum_order = 9

      character(len=maxpop), parameter ::
     &     pert_ops = 'drpvam'          ! d: dipole     r: position
                                        ! p: momentum   v: dipole velocity
                                        ! a: angular momentum
                                        ! m: magnetic moment
      character(len=6*maxpop), parameter ::
     &     dalton_int_names = 'DIPLENDIPLENDIPVELDIPVELANGMOMANGMOM'
      integer, parameter ::
     &     pert_op_sign(maxpop) = (/2,1,4,3,3,8/)

*----------------------------------------------------------------------*
*     pert_op_info and pert_component_info definitions
*----------------------------------------------------------------------*
      type pert_op_info
        character(len=1) ::
     &       name,           ! one of pert_ops
     &       comp            ! directional component
        character(len=8) ::
     &       int_name        ! name of dalton integral list
        integer ::
     &       isym,           ! irrep
     &       sign            ! 1=+, 2=-, 3=i, 4=-i,5=1/2,...,8=-i/2
      end type pert_op_info

      type pert_component_info
        integer ::
     &       redun,          ! refers to highest idx of similar component
     &       pop_idx,         ! idx of assigned pert. op.
     &       order           ! order of the response parameters associated with this component
        real(8) ::
     &       freq            ! assigned frequency
      end type pert_component_info
