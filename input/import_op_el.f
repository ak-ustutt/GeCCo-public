*----------------------------------------------------------------------*
      subroutine import_op_el(idxop,idxffop,ffops,ops,nops,
     &     env_type,str_info,orb_info)
*----------------------------------------------------------------------*
*     import matrix elements for operator idxop from environment
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'

      integer, intent(in) ::
     &     idxop, idxffop, nops
      type(operator), intent(in) ::
     &     ops(nops)
      type(filinf), intent(in) ::
     &     ffops(nops)
      character, intent(in) ::
     &     env_type*(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      select case(trim(env_type))
      case ('dalton','DALTON')
        ! what to import?
        select case(trim(ops(idxop)%name))
        case ('H')
          call import_hamint_dalton(ops(idxop),ffops(idxffop),
     &         str_info,orb_info)
        case default
          call quit(1,'import_op_el','DALTON: cannot handle operator '
     &         //trim(ops(idxop)%name))
        end select
      case ('intern','INTERN')
        call quit(1,'import_op_el','type INTERN not implemented')
      case ('aces2','ACES2')
        call quit(1,'import_op_el','type ACES2 not implemented')
      case ('tmole','TMOLE')
        call quit(1,'import_op_el','type TMOLE not implemented')
      case default
        call quit(1,'import_op_el','unknown type '//trim(env_type))
      end select

      return
      end
