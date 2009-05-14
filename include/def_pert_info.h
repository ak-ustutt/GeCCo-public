*----------------------------------------------------------------------*
*     parameters
*----------------------------------------------------------------------*
      integer, parameter ::
     &     maxpop = 4, maxcmp = 99, maximum_order = 9

      character(len=maxpop), parameter ::
     &     pert_ops = 'urpl'
      character(len=6*maxpop), parameter ::
     &     dalton_int_names = 'DIPLENDIPLENDIPVELANGMOM'
      integer, parameter ::
     &     pert_op_sign(maxpop) = (/2,1,3,3/)

*----------------------------------------------------------------------*
*     pert_op_info and pert_component_info definitions
*----------------------------------------------------------------------*
      type pert_op_info
        character(len=1) ::
     &       name,           ! e.g. r, p (momentum), l (ang. mom.)
     &       comp            ! directional component
        character(len=8) ::
     &       int_name        ! name of dalton integral list
        integer ::
     &       isym,           ! irrep
     &       sign            ! 1=+, 2=-, 3=i, 4=-i
      end type pert_op_info

      type pert_component_info
        integer ::
     &       redun,          ! refers to highest idx of similar component
     &       pop_idx         ! idx of assigned pert. op.
        real(8) ::
     &       freq            ! assigned frequency
      end type pert_component_info
