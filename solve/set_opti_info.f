*----------------------------------------------------------------------*
      subroutine set_opti_info(opti_info,
     &     mode,nopt,nroot,mel_opt,prc_str)
*----------------------------------------------------------------------*
      implicit none
      
      include 'def_optimize_info.h'
      include 'ifc_input.h'
      include 'ifc_memman.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_me_list_array.h'

      type(optimize_info), intent(inout) ::
     &     opti_info
      integer, intent(in) ::
     &     nopt, mode, nroot
      type(me_list_array), intent(in) ::
     &     mel_opt(nopt)
      character(*), intent(in) ::
     &     prc_str

      character ::
     &     str*16
      integer ::
     &     ifree, iopt, ioff, n_states, ciroot
      logical ::
     &     multistate

      integer, external ::
     &     ndisblk_mel

      opti_info%shift = 0d0
      opti_info%variational = .false.
      if (mode.eq.1) opti_info%linear = .false.
      if (mode.eq.2) opti_info%linear = .true.
      if (mode.eq.3) opti_info%linear = .true.

      ! get global values
      call get_argument_value('calculate.solve','maxiter',
     &     ival=opti_info%maxmacit)
      call get_argument_value('calculate.solve','maxsub',
     &     ival=opti_info%maxsbsp)
      call get_argument_value('calculate.solve','check_incore',
     &     ival=opti_info%max_incore)

      ifree = mem_alloc_int (opti_info%nwfpar,nopt,'nwfpar')
      ifree = mem_alloc_int (opti_info%typ_prc,nopt,'typprc')
      ifree = mem_alloc_real(opti_info%thrgrd,nopt,'thrgrd')

      call get_argument_value('calculate.solve','conv',
     &     xval=opti_info%thrgrd(1))
      if (nopt.gt.1)
     &     opti_info%thrgrd(2:nopt) = opti_info%thrgrd(1)

      opti_info%nopt = nopt
      do iopt = 1, nopt
        opti_info%nwfpar(iopt) = mel_opt(iopt)%mel%len_op
      end do

      opti_info%nroot = nroot

      ! specials:
      if (mode.eq.1) then

        if (nroot.gt.1)
     &       call quit(1,'set_opti_info',
     &       'nroot>1 not yet considered for non-linear equations!')

        if(is_argument_set('calculate.solve.non_linear','maxiter').gt.0)
     &       call get_argument_value('calculate.solve.non_linear',
     &       'maxiter',
     &       ival=opti_info%maxmacit)
        if (is_argument_set('calculate.solve.non_linear','maxsub').gt.0)
     &       call get_argument_value('calculate.solve.non_linear',
     &       'maxsub',
     &       ival=opti_info%maxsbsp)

        call get_argument_value('calculate.solve.non_linear','maxmic',
     &       ival=opti_info%micifac)
        opti_info%maxmicit = opti_info%maxmacit*opti_info%micifac

        if (is_argument_set('calculate.solve.non_linear','conv').gt.0) 
     &    then
          call get_argument_value('calculate.solve.non_linear','conv',
     &         xval=opti_info%thrgrd(1))
          if (nopt.gt.1)
     &         opti_info%thrgrd(2:nopt) = opti_info%thrgrd(1)
        end if

        str(1:16) = ' '
        call get_argument_value('calculate.solve.non_linear','method',
     &       str=str)
        call uppcas(str)
        select case(trim(str))
        case('PERT') 
          opti_info%mode_nleq = mode_nleq_pert
          opti_info%norder = 1
        case('DIIS') 
          opti_info%mode_nleq = mode_nleq_diis
          opti_info%norder = 1
        case('ASSJ','RLE') 
          opti_info%mode_nleq = mode_nleq_assj
          opti_info%norder = 1
        case default
          call quit(0,'set_opti_info','invalid method: '//trim(str))
        end select

        call get_argument_value('method.MR','ciroot',
     &       ival=ciroot)
        call get_argument_value('method.MR','multistate',
     &       lval=multistate)
        if(multistate)then
         n_states = ciroot
         if(n_states.lt.2)
     &        call quit(1,'set_opti_info',
     &        'Please, do not use multistate=T for one state!')
        else
         n_states = 1
        end if

        call get_argument_value('calculate.solve.non_linear','tr_ini',
     &       xval=opti_info%trini)

        call get_argument_value('calculate.solve.non_linear','optref',
     &       ival=opti_info%optref)
        call get_argument_value('calculate.solve.non_linear',
     &       'update_prc',
     &       ival=opti_info%update_prc)
        call get_argument_value('calculate.solve.non_linear',
     &       'preopt',
     &       lval=opti_info%skip_resx)
        opti_info%skip_resx = opti_info%skip_resx.and.nopt.eq.n_states

        call get_argument_value('calculate.solve.non_linear',
     &       'mic_ahead',
     &       xval=opti_info%mic_ahead)

      else if (mode.eq.2) then

        if (is_argument_set('calculate.solve.linear','maxiter').gt.0)
     &       call get_argument_value('calculate.solve.linear',
     &       'maxiter',
     &       ival=opti_info%maxmacit)
        if (is_argument_set('calculate.solve.linear','maxsub').gt.0)
     &       call get_argument_value('calculate.solve.linear',
     &       'maxsub',
     &       ival=opti_info%maxsbsp)

        ! here the interpretation is: maxsbsp = sbsp per root, so:
        opti_info%maxsbsp = opti_info%maxsbsp*nroot

        if (is_argument_set('calculate.solve.linear','conv').gt.0) then
          call get_argument_value('calculate.solve.linear','conv',
     &         xval=opti_info%thrgrd(1))
          if (nopt.gt.1)
     &         opti_info%thrgrd(2:nopt) = opti_info%thrgrd(1)
        end if

        call get_argument_value('calculate.solve.linear','method',
     &       str=str)
        call uppcas(str)
        select case(trim(str(1:8)))
        case('CONJGRAD') 
          opti_info%mode_leq = mode_leq_conjg
          opti_info%norder = 1
        case('SUBSPACE') 
          opti_info%mode_leq = mode_leq_subsp
          opti_info%norder = 1
        case default
          call quit(0,'set_opti_info','invalid method: '//trim(str))
        end select
     
        ! set shift
        opti_info%shift = mel_opt(1)%mel%frequency
        do iopt = 2,nopt
          if (mel_opt(iopt)%mel%frequency.ne.opti_info%shift)
     &      call quit(0,'set_opti_info'
     &                 ,'shift associated with solution '//
     &                'vectors must be identical for all of them')
        end do
 
      else if (mode.eq.3) then

        if (is_argument_set('calculate.solve.eigen','maxiter').gt.0)
     &       call get_argument_value('calculate.solve.eigen',
     &       'maxiter',
     &       ival=opti_info%maxmacit)
        if (is_argument_set('calculate.solve.eigen','maxsub').gt.0)
     &       call get_argument_value('calculate.solve.eigen',
     &       'maxsub',
     &       ival=opti_info%maxsbsp)

        ! here the interpretation is: maxsbsp = sbsp per root, so:
        opti_info%maxsbsp = opti_info%maxsbsp*nroot

        if (is_argument_set('calculate.solve.eigen','conv').gt.0) then
          call get_argument_value('calculate.solve.eigen','conv',
     &         xval=opti_info%thrgrd(1))
          if (nopt.gt.1)
     &         opti_info%thrgrd(2:nopt) = opti_info%thrgrd(1)
        end if

        call get_argument_value('calculate.solve.eigen','method',
     &       str=str)
        call uppcas(str)
        select case(trim(str(1:8)))
        case('DAVIDSON') 
          opti_info%mode_evp = mode_leq_subsp
          opti_info%norder = 1
        case default
          call quit(0,'set_opti_info','invalid method: '//trim(str))
        end select

        call get_argument_value('calculate.solve.eigen','resume',
     &       lval=opti_info%resume)
      
      else
        call quit(1,'set_opti_info','illegal value of mode')
      end if

      ! set type of preconditioner
      ioff = 0
      do iopt = 1, nopt
        select case(prc_str(ioff+1:ioff+3))
        case('DIA') 
          opti_info%typ_prc(iopt) = optinf_prc_file
        case('BLK')
          opti_info%typ_prc(iopt) = optinf_prc_blocked
        case('MIX')
          opti_info%typ_prc(iopt) = optinf_prc_mixed
        case('TRF') 
          opti_info%typ_prc(iopt) = optinf_prc_traf
        case('TR0') 
          opti_info%typ_prc(iopt) = optinf_prc_traf_spc
        case('IH0') 
          opti_info%typ_prc(iopt) = optinf_prc_invH0
        case('NRM')
          opti_info%typ_prc(iopt) = optinf_prc_norm
        case('SPP')
          opti_info%typ_prc(iopt) = optinf_prc_spinp
        case('PRJ')
          opti_info%typ_prc(iopt) = optinf_prc_prj
        case('SRP')
          opti_info%typ_prc(iopt) = optinf_prc_spinrefp
        case default
          call quit(1,'set_opti_info','cannot interpret string: '//
     &         trim(prc_str))
        end select
        ioff = ioff+4
      end do

      return
      end
