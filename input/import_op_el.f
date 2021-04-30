*----------------------------------------------------------------------*
      subroutine import_op_el(label_mel,op_label,env_type,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     import matrix elements from environment to ME-list with
*     label "label_mel" 
*
*     op_label contains the type of operator to be imported,
*     env_type is the environment type (Molpro, Dalton, so on)
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000
      logical, parameter ::
     &     hard_restart = .false. ! currently for testing only

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'
      include 'def_strmapinf.h'

      character(len=*), intent(in) ::
     &     label_mel
      type(operator_info), intent(in) ::
     &     op_info
      character(len=*), intent(in) ::
     &     env_type, op_label
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info

      integer, parameter ::
     &     use_scaling = 10  ! use 0 to turn off scaling

      integer ::
     &     ipri, mode, scaling, idx_mel
      real(8) ::
     &     xdum
      logical ::
     &     anti, list_exists, act_orbs, watch
      type(me_list), pointer ::
     &     mel_target

      integer, external ::
     &     idx_mel_list

      if (use_scaling.lt.2.and.use_scaling.ne.0)
     &     call quit(1,'import_op_el',
     &     'use use_scaling=0 to switch off scaling')

      idx_mel = idx_mel_list(label_mel,op_info)
      if (idx_mel.lt.0)
     &     call quit(1,'import_op_el','Label not on list: "'//
     &     trim(label_mel)//'"')

      mel_target => op_info%mel_arr(idx_mel)%mel

      anti = .true.
      act_orbs = orb_info%nactt_hpv(IVALE).gt.0

      call touch_file_rec(mel_target%fhand)

      ! hard restart
      inquire(file=trim(mel_target%fhand%name),exist=list_exists)
      if (list_exists.and.hard_restart) then
        write(lulog,*) 'found: ',trim(mel_target%fhand%name)
        write(lulog,*) 'trying a hard restart ... watch out!'
        return
      end if

      select case(trim(env_type))
      case ('dalton_special','DALTON_SPECIAL')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT','G_INT','F_INT')
          if(trim(mel_target%op%name).eq.op_g_z) anti = .false.
          mode=1
          scaling=min(use_scaling,0)
          if (trim(op_label).eq.'G_INT'.or.
     &        trim(op_label).eq.'H_INT') then
            if (anti) then
              call import_2el_dalton_new(mel_target,'MO_G',
     &             mode,scaling,anti,str_info,orb_info)
            else ! G_Z still works with old import only:
              call import_2el_dalton(mel_target,'MO_G',
     &             mode,scaling,anti,str_info,orb_info)
            end if
          end if
          ! call after 2int import, as the above routine
          ! zeroes all blocks, including E0 and F:
          if (trim(op_label).eq.'FI_INT'.or.
     &        (trim(op_label).eq.'H_INT'.and.act_orbs))
     &       call import_fock_dalton(mel_target,str_info,orb_info,
     &         .true.,'MO_FI')
          if (trim(op_label).eq.'F_INT'.or.
     &        (trim(op_label).eq.'H_INT'.and..not.act_orbs))
     &       call import_fock_dalton(mel_target,str_info,orb_info,
     &         .true.,'MO_F')
        ! Integrals needed for R12 calculations.
        case('F12_INT')
          mode=1
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_F12',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12BAR_INT')
          mode=3
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_F12BAR',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FDGBAR_INT')
          mode=3
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_FDGBAR',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12TLD_INT')
          mode=3
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_F12TLD',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12BRV_INT')
          mode=3
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_F12BRV',
     &           mode,scaling,anti,str_info,orb_info) 
        case('FF_INT')
          mode=1
          scaling=min(use_scaling,2)
          call import_2el_dalton_new(mel_target,'MO_FF',
     &         mode,scaling,anti,str_info,orb_info)

        case('FFBAR_INT')
          mode=3
          scaling=min(use_scaling,2)
          call import_2el_dalton_new(mel_target,'MO_FFBAR',
     &         mode,scaling,anti,str_info,orb_info)

        case('FFG_INT')
          mode=1
          scaling=min(use_scaling,2)
          call import_2el_dalton_new(mel_target,'MO_FFG',
     &         mode,scaling,anti,str_info,orb_info)

        case('TTF_INT')
          mode=2
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_TTF',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FTF_INT')
          mode=1
          scaling=min(use_scaling,2)
          call import_2el_dalton_new(mel_target,'MO_FTF',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FG_INT')
          mode=1
          scaling=min(use_scaling,1)
          call import_2el_dalton_new(mel_target,'MO_FG',
     &           mode,scaling,anti,str_info,orb_info) 

        case('K_INT')
          call quit(1,'import_op_el','K is not yet ready')
c          call import_exchange_dalton(mel_target,'MO_K',
c     &                                str_info,orb_info)

        case('Z_LIST','P_LIST')
          call import_intm_fc(mel_target,mel_target%op%name,
     &         str_info,orb_info)

        case default
        ! Integrals for general one-electron operators
          call import_propint(mel_target,op_label,env_type,
     &         str_info,orb_info)

        end select
          
      case ('dalton','DALTON')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          call import_hamint_dalton(mel_target,str_info,orb_info)

        case default
        ! Integrals for general one-electron operators
          call import_propint(mel_target,op_label,env_type,
     &         str_info,orb_info)

        end select
      case ('dalton64','DALTON64')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          call import_hamint_dalton64(mel_target,str_info,orb_info)

        case default
        ! Integrals for general one-electron operators
          call import_propint(mel_target,op_label,env_type,
     &         str_info,orb_info)

        end select

      case ('gamess','GAMESS')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          call import_hamint_gamess(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el','GAMESS: cannot handle op_label "'
     &         //trim(op_label)//'"')
        end select

      case ('molpro_dump','MOLPRO_DUMP')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          call import_hamint_molpro_dump(mel_target,str_info,orb_info)

        case default
        ! Integrals for general one-electron operators
          call import_propint(mel_target,op_label,env_type,
     &         str_info,orb_info)

        end select

      case ('molpro_ifc','MOLPRO_IFC')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          anti = .true.
          call import_hamint_molpro(mel_target,anti,str_info,orb_info)
        case ('Hnox_INT')
          anti = .false.
          call import_hamint_molpro(mel_target,anti,str_info,orb_info)

        case default
        ! Integrals for general one-electron operators
          call import_propint(mel_target,op_label,env_type,
     &         str_info,orb_info)

        end select

      case ('cfour','CFOUR')
        ! what to import?
        select case(trim(op_label))
        case ('H_INT')
          call import_hamint_cfour(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el',
     &         'CFOUR: cannot handle op_label "'
     &         //trim(op_label)//'"')
        end select

      case ('intern','INTERN')
        call quit(1,'import_op_el','type INTERN not implemented')
      case ('aces2','ACES2')
        call quit(1,'import_op_el','type ACES2 not implemented')
      case ('tmole','TMOLE')
        call quit(1,'import_op_el','type TMOLE not implemented')
      case default
        call quit(1,'import_op_el',
     &                    'unknown type "'//trim(env_type)//'"')
      end select

      !debug: set here a specific target to get more output
      watch = .false.
      !watch = trim(op_label).eq.'Hnox_INT'
      if (ntest.ge.10.and.(.not.mel_target%op%formal).or.watch) then
        write(lulog,*)
        write(lulog,*) 'imported list: ',trim(mel_target%label)
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000.or.watch) ipri = 5
        call wrt_mel_file(lulog,ipri,mel_target,
     &       1,mel_target%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      return
      end
