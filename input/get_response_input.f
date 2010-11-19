*----------------------------------------------------------------------*
      subroutine get_response_input(ncnt,maxord,ncmp,cmp,npop,pop,
     &                              orb_info)
*----------------------------------------------------------------------*
*     process response functions input
*
*     matthias, 2009
*----------------------------------------------------------------------*
      implicit none

      include 'def_target.h'
      include 'def_orbinf.h'
      include 'ifc_input.h'
      include 'def_pert_info.h'

      integer, intent(in) ::
     &     ncnt
      type(orbinf), intent(in) ::
     &     orb_info

      integer, intent(out) ::
     &     maxord(ncnt), ncmp, npop
      type(pert_op_info), intent(out) ::
     &     pop(3*maxpop)
      type(pert_component_info), intent(out) ::
     &     cmp(maxcmp)
      
      integer ::
     &     icnt, idx, jdx, ipop, pos, sign
      logical ::
     &     skip
      character(len=1) ::
     &     pert_ord
      character(len=6) ::
     &     int_name
      character(len_command_par) ::
     &     pert, pertop

      integer, external ::
     &     pert_sym

      do icnt = 1,ncnt
        call get_argument_value('method.response','order',
     &       keycount=icnt,ival=maxord(icnt))
      end do
      ncmp = ncnt*maxval(maxord)
      cmp(:)%redun = 0
      cmp(:)%pop_idx = 0
      cmp(:)%freq = 0d0
      if (maxval(maxord).gt.maximum_order) then
        write(pert_ord,'(i1)') maximum_order
        call quit(1,'get_response_input',
     &        'maxord must not exceed '//pert_ord)
      end if
      pert(1:len_command_par) = ' '
      pertop(1:len_command_par) = ' '

      npop = 0
      do icnt = 1,ncnt
        pos = (icnt-1)*maxval(maxord) + 1
        call get_argument_value('method.response','comp',
     &       keycount=icnt,str=pert(pos:len_command_par))
        call get_argument_value('method.response','pert',
     &       keycount=icnt,str=pertop(pos:len_command_par))

        call get_argument_value('method.response','freq',
     &       keycount=icnt,xarr=cmp(pos+1:ncmp)%freq)
        cmp(pos+maxord(icnt):ncmp)%freq = 0d0
        if (maxord(icnt).gt.0)
     &       cmp(pos)%freq = -sum(cmp(pos+1:pos+maxord(icnt)-1)%freq)

        ! duplicate values for pert if necessary and not specified
        do idx = pos+1,pos+maxord(icnt)-1
          if (pert(idx:idx).eq.' ')
     &        pert(idx:idx) = pert(idx-1:idx-1)
          if (pertop(idx:idx).eq.' ')
     &        pertop(idx:idx) = pertop(idx-1:idx-1)
        end do

        do idx = pos,pos+maxord(icnt)-1
          ! check if perturbation operator input is ok
          if (pert(idx:idx).ne.'X' .and.
     &          pert(idx:idx).ne.'Y' .and.
     &          pert(idx:idx).ne.'Z')
     &          call quit(1,'get_response_input',
     &          'comp must contain X,Y,Z')
          skip = .false.
          do ipop = 1,len(pert_ops)
            if (pertop(idx:idx).eq.pert_ops(ipop:ipop)) then
              skip = .true.
              int_name = dalton_int_names(6*ipop-5:6*ipop)
              sign = pert_op_sign(ipop)
              exit
            end if
          end do
          if (.not.skip) call quit(1,'get_response_input',
     &          'perturbation operator "'//pertop(idx:idx)//
     &          '" is currently not allowed.')
          skip = .false.
          do ipop = 1,npop
            if (pop(ipop)%comp.eq.pert(idx:idx).and.
     &              pop(ipop)%name.eq.pertop(idx:idx)) then
              skip = .true.
              cmp(idx)%pop_idx = ipop
              exit
            end if
          end do

          if (.not.skip) then
            npop = npop + 1
            pop(npop)%comp = pert(idx:idx)
            pop(npop)%name = pertop(idx:idx)
            pop(npop)%int_name = pert(idx:idx)//int_name//' '
            pop(npop)%sign = sign
            pop(npop)%isym = pert_sym(pop(npop)%int_name,orb_info)
            cmp(idx)%pop_idx = npop
          end if

          ! determine redundancies
          do jdx = 1,idx-1
            if (pert(idx:idx).eq.pert(jdx:jdx) .and.
     &          pertop(idx:idx).eq.pertop(jdx:jdx) .and.
     &          abs(cmp(idx)%freq - cmp(jdx)%freq).lt.1d-12)
     &         cmp(idx)%redun = jdx
          end do
          if (cmp(idx)%redun.eq.0) cmp(idx)%redun = idx
        end do
      end do

      return
      end
