*----------------------------------------------------------------------*
      subroutine import_op_el(label_mel,
     &     op_info,
     &     env_type,str_info,orb_info)
*----------------------------------------------------------------------*
*     import matrix elements from environment for ME-list with 
*     label "label_mel" 
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 10

      include 'stdunit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'

      character(*), intent(in) ::
     &     label_mel
      type(operator_info), intent(in) ::
     &     op_info
      character, intent(in) ::
     &     env_type*(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer ::
     &     ipri, mode, idx_mel
      type(me_list), pointer ::
     &     mel_target

      integer, external ::
     &     idx_mel_list

      idx_mel = idx_mel_list(label_mel,op_info)
      if (idx_mel.lt.0)
     &     call quit(1,'import_op_el','Label not on list: "'//
     &     trim(label_mel)//'"')

      mel_target => op_info%mel_arr(idx_mel)%mel

      select case(trim(env_type))
      case ('dalton_special','DALTON_SPECIAL')
        ! what to import?
        select case(trim(mel_target%op%name))
        case (op_ham)
          mode=0
          call import_r12_dalton(mel_target,'MO_G',
     &           mode,str_info,orb_info)
          ! call after 2int import, as the above routine
          ! zeroes all blocks, including E0 and F:
          call import_fock_dalton(mel_target,str_info,orb_info)
        ! Get other integrals needed for R12 calculations.
        case(op_rint)
          mode=1
          call import_r12_dalton(mel_target,'MO_R',
     &           mode,str_info,orb_info) 

        case(op_rinba)
          mode=1
          call import_r12_dalton(mel_target,'MO_R',
     &           mode,str_info,orb_info) 

        case(op_f2)
          mode=1
          call import_r12_dalton(mel_target,'MO_F2',
     &         mode,str_info,orb_info)

        case(op_ttr)
          mode=2
          call import_r12_dalton(mel_target,'MO_TTR',
     &         mode,str_info,orb_info) 

        case(op_exchange)
c          call import_exchange_dalton(mel_target,'MO_K',
c     &                                str_info,orb_info)
c dbg
          print *,'Hello'
c dbg

        case default
          call quit(1,'import_op_el',
     &         'DALTON_SPECIAL: cannot handle operator '
     &         //trim(mel_target%op%name))
        end select
          
      case ('dalton','DALTON')
        ! what to import?
        select case(trim(mel_target%op%name))
        case (op_ham)
          call import_hamint_dalton(mel_target,str_info,orb_info)

        case default
          call quit(1,'import_op_el','DALTON: cannot handle operator '
     &         //trim(mel_target%op%name))
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

      if (ntest.ge.10.and.(.not.mel_target%op%formal)) then
        write(luout,*)
        write(luout,*) 'imported list: ',trim(mel_target%label)
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_mel_file(luout,ipri,mel_target,
     &       1,mel_target%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      return
      end
