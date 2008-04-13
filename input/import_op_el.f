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
     &     ipri, mode, scaling, idx_mel
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
        case (op_ham,op_g_x)
          mode=1
          scaling=0
          call import_2el_dalton(mel_target,'MO_G',
     &           mode,scaling,str_info,orb_info)
          ! call after 2int import, as the above routine
          ! zeroes all blocks, including E0 and F:
          if (trim(mel_target%op%name).eq.op_ham)
     &       call import_fock_dalton(mel_target,str_info,orb_info,
     &         .true.,'MO_F')
        ! Get other integrals needed for R12 calculations.
        case(op_rint)
          mode=1
          scaling=1
          call import_2el_dalton(mel_target,'MO_F12',
     &           mode,scaling,str_info,orb_info) 

        case(op_rintbar)
          mode=3
          scaling=1
          call import_2el_dalton(mel_target,'MO_F12BAR',
     &           mode,scaling,str_info,orb_info) 

        case(op_rdagbar)
          mode=3
          scaling=1
          call import_2el_dalton(mel_target,'MO_FDGBAR',
     &           mode,scaling,str_info,orb_info) 

        case(op_rinttilde)
          mode=3
          scaling=1
          call import_2el_dalton(mel_target,'MO_F12TLD',
     &           mode,scaling,str_info,orb_info) 

        case(op_rintbreve)
          mode=3
          scaling=1
          call import_2el_dalton(mel_target,'MO_F12BRV',
     &           mode,scaling,str_info,orb_info) 

c        case(op_rintc)
c          mode=3
c          call import_2el_dalton(mel_target,'MO_F12C',
c     &           mode,scaling,str_info,orb_info) 

        case(op_rinba)
          call quit(1,'import_op_el',
     &         'import of F12^+ is obsolete')
          mode=1
          scaling=1
          call import_2el_dalton(mel_target,'MO_F12',
     &           mode,scaling,str_info,orb_info) 

        case(op_ff)
          mode=1
          scaling=2
          call import_2el_dalton(mel_target,'MO_FF',
     &         mode,scaling,str_info,orb_info)

        case(op_ffbar)
          mode=3
          scaling=2
          call import_2el_dalton(mel_target,'MO_FFBAR',
     &         mode,scaling,str_info,orb_info)

        case(op_ttr)
          mode=2
          scaling=1
          call import_2el_dalton(mel_target,'MO_TTF',
     &           mode,scaling,str_info,orb_info) 

        case(op_rttr)
          mode=1
          scaling=2
          call import_2el_dalton(mel_target,'MO_FTF',
     &           mode,scaling,str_info,orb_info) 

        case(op_gr)
          mode=1
          scaling=1
          call import_2el_dalton(mel_target,'MO_FG',
     &           mode,scaling,str_info,orb_info) 

        case(op_exchange)
          call quit(1,'import_op_el','K is not yet ready')
c          call import_exchange_dalton(mel_target,'MO_K',
c     &                                str_info,orb_info)

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
c dbg
c      stop 'test import'
c dbg

      return
      end
