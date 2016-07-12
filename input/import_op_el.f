*----------------------------------------------------------------------*
      subroutine import_op_el(label_mel,
     &     list_type,env_type,trplt,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
*     import matrix elements from environment for ME-list with 
*     label "label_mel" 
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
     &     env_type, list_type
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      logical, intent(in) ::
     &     trplt

      integer, parameter ::
     &     use_scaling = 10  ! use 0 to turn off scaling

      integer ::
     &     ipri, mode, scaling, idx_mel
      real(8) ::
     &     xdum
      logical ::
     &     anti, list_exists, act_orbs
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
        select case(trim(list_type))
        case ('H_INT','G_INT','F_INT')
          if(trim(mel_target%op%name).eq.op_g_z) anti = .false.
          mode=1
          scaling=min(use_scaling,0)
          if (trim(list_type).eq.'G_INT'.or.
     &        trim(list_type).eq.'H_INT') then
c     &         call import_2el_dalton(mel_target,'MO_G',
c     &           mode,scaling,anti,str_info,orb_info)
c     &         call import_2el_dalton(mel_target,'MO_G',
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
          if (trim(list_type).eq.'FI_INT'.or.
     &        (trim(list_type).eq.'H_INT'.and.act_orbs))
     &       call import_fock_dalton(mel_target,str_info,orb_info,
     &         .true.,'MO_FI')
          if (trim(list_type).eq.'F_INT'.or.
     &        (trim(list_type).eq.'H_INT'.and..not.act_orbs))
     &       call import_fock_dalton(mel_target,str_info,orb_info,
     &         .true.,'MO_F')
        ! Get other integrals needed for R12 calculations.
        case('F12_INT')
          mode=1
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_F12',
          call import_2el_dalton_new(mel_target,'MO_F12',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12BAR_INT')
          mode=3
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_F12BAR',
          call import_2el_dalton_new(mel_target,'MO_F12BAR',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FDGBAR_INT')
          mode=3
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_FDGBAR',
          call import_2el_dalton_new(mel_target,'MO_FDGBAR',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12TLD_INT')
          mode=3
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_F12TLD',
          call import_2el_dalton_new(mel_target,'MO_F12TLD',
     &           mode,scaling,anti,str_info,orb_info) 

        case('F12BRV_INT')
          mode=3
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_F12BRV',
          call import_2el_dalton_new(mel_target,'MO_F12BRV',
     &           mode,scaling,anti,str_info,orb_info) 
        case('FF_INT')
          mode=1
          scaling=min(use_scaling,2)
c          call import_2el_dalton(mel_target,'MO_FF',
          call import_2el_dalton_new(mel_target,'MO_FF',
     &         mode,scaling,anti,str_info,orb_info)

        case('FFBAR_INT')
          mode=3
          scaling=min(use_scaling,2)
c          call import_2el_dalton(mel_target,'MO_FFBAR',
          call import_2el_dalton_new(mel_target,'MO_FFBAR',
     &         mode,scaling,anti,str_info,orb_info)

        case('FFG_INT')
          mode=1
          scaling=min(use_scaling,2)
c          call import_2el_dalton(mel_target,'MO_FFG',
          call import_2el_dalton_new(mel_target,'MO_FFG',
     &         mode,scaling,anti,str_info,orb_info)

        case('TTF_INT')
          mode=2
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_TTF',
          call import_2el_dalton_new(mel_target,'MO_TTF',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FTF_INT')
          mode=1
          scaling=min(use_scaling,2)
c          call import_2el_dalton(mel_target,'MO_FTF',
          call import_2el_dalton_new(mel_target,'MO_FTF',
     &           mode,scaling,anti,str_info,orb_info) 

        case('FG_INT')
          mode=1
          scaling=min(use_scaling,1)
c          call import_2el_dalton(mel_target,'MO_FG',
          call import_2el_dalton_new(mel_target,'MO_FG',
     &           mode,scaling,anti,str_info,orb_info) 

        case('K_INT')
          call quit(1,'import_op_el','K is not yet ready')
c          call import_exchange_dalton(mel_target,'MO_K',
c     &                                str_info,orb_info)

        case('Z_LIST','P_LIST')
          call import_intm_fc(mel_target,mel_target%op%name,
     &         str_info,orb_info)

        case ('XDIPLEN','YDIPLEN','ZDIPLEN')
          call import_propint_dalton(mel_target,list_type,1,trplt,
     &         str_info,orb_info)

        case ('XDIPVEL','YDIPVEL','ZDIPVEL',
     &        'XANGMOM','YANGMOM','ZANGMOM')
          call import_propint_dalton(mel_target,list_type,-1,trplt,
     &         str_info,orb_info)

        case default
          call quit(1,'import_op_el',
     &         'DALTON_SPECIAL: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
          
      case ('dalton','DALTON')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_dalton(mel_target,str_info,orb_info)
        case ('XDIPLEN','YDIPLEN','ZDIPLEN')
          call import_propint_dalton(mel_target,list_type,1,trplt,
     &         str_info,orb_info)
        case ('XDIPVEL','YDIPVEL','ZDIPVEL',
     &        'XANGMOM','YANGMOM','ZANGMOM')
          call import_propint_dalton(mel_target,list_type,-1,trplt,
     &         str_info,orb_info)

        case default
          call quit(1,'import_op_el','DALTON: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
      case ('dalton64','DALTON64')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_dalton64(mel_target,str_info,orb_info)
        case ('XDIPLEN','YDIPLEN','ZDIPLEN')
          call import_propint_dalton(mel_target,list_type,1,trplt,
     &         str_info,orb_info)
        case ('XDIPVEL','YDIPVEL','ZDIPVEL',
     &        'XANGMOM','YANGMOM','ZANGMOM')
          call import_propint_dalton(mel_target,list_type,-1,trplt,
     &         str_info,orb_info)

        case default
          call quit(1,'import_op_el','DALTON: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
      case ('gamess','GAMESS')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_gamess(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el','GAMESS: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
      case ('molpro_dump','MOLPRO_DUMP')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_molpro_dump(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el',
     &         'MOLPRO_DUMP: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
      case ('molpro_ifc','MOLPRO_IFC')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_molpro(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el',
     &         'MOLPRO_IFC: cannot handle list_type "'
     &         //trim(list_type)//'"')
        end select
      case ('cfour','CFOUR')
        ! what to import?
        select case(trim(list_type))
        case ('H_INT')
          call import_hamint_cfour(mel_target,str_info,orb_info)
        case default
          call quit(1,'import_op_el',
     &         'CFOUR: cannot handle list_type "'
     &         //trim(list_type)//'"')
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

      if (ntest.ge.10.and.(.not.mel_target%op%formal)) then
        write(lulog,*)
        write(lulog,*) 'imported list: ',trim(mel_target%label)
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
c dbg
c        if (trim(list_type).eq.'FF_INT') ipri = 5
c dbg
        call wrt_mel_file(lulog,ipri,mel_target,
     &       1,mel_target%op%n_occ_cls,
     &       str_info,orb_info)
c dbg
c       if (trim(list_type).eq.'H_INT') stop 'H stop'
c dbg
      end if

      return
      end
