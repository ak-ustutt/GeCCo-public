*----------------------------------------------------------------------*
      subroutine import_op_el(idxop,idxffop,
     &     op_info,
     &     env_type,str_info,orb_info)
*----------------------------------------------------------------------*
*     import matrix elements for operator idxop from environment
*     note: idxffop is obsolete now
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 1000

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'
      include 'explicit.h'

      integer, intent(in) ::
     &     idxop, idxffop
      type(operator_info), intent(in) ::
     &     op_info
      character, intent(in) ::
     &     env_type*(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ipri,mode
      type(operator), pointer ::
     &     op_target
      type(filinf), pointer ::
     &     opfil_target

      op_target => op_info%op_arr(idxop)%op
      ! init file for operator with standard name
      !                opfil_target is dummy
      call assign_file_to_op(idxop,.true.,opfil_target,
     &                       1,1,1,
     &                       1,op_info)

      opfil_target => op_info%opfil_arr(idxop)%fhand

      select case(trim(env_type))
      case ('dalton','DALTON')
        ! what to import?
        select case(trim(op_target%name))
        case (op_ham)
          call import_hamint_dalton(op_target,opfil_target,
     &         str_info,orb_info)
          if(explicit)then
            mode=0
            call import_r12_dalton(op_target,opfil_target,'MO_G',
     &           mode,str_info,orb_info)
          endif

        ! Get other integrals needed for R12 calculations.
        case(op_rint)
          if(.not.op_target%formal)then
            mode=1
            call import_r12_dalton(op_target,opfil_target,'MO_R',
     &           mode,str_info,orb_info) 
          else
            write(luout,*)'R12 operator is purely formal.'
          endif

        case(op_rinba)
          if(.not.op_target%formal)then
            mode=1
            call import_r12_dalton(op_target,opfil_target,'MO_R',
     &           mode,str_info,orb_info) 
          else
            write(luout,*)'R12-DAG operator is purely formal.'
          endif

        case(op_ttr)
          mode=2
          call import_r12_dalton(op_target,opfil_target,'MO_TTR',
     &         mode,str_info,orb_info) 
          
c        case(op_del_inter)
c          mode=1
c          call import_r12_dalton(op_target,opfil_target,'DELTA',
c     &         mode,str_info,orb_info)

        case default
          call quit(1,'import_op_el','DALTON: cannot handle operator '
     &         //trim(op_target%name))
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

      if (ntest.ge.10.and.(.not.op_target%formal)) then
        write(luout,*)
        write(luout,*) 'imported operator: ',trim(op_target%name)
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_op_file(luout,ipri,opfil_target,op_target,
     &       1,op_target%n_occ_cls,
     &       str_info,orb_info)
      end if

      return
      end
