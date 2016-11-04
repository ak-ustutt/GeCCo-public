*----------------------------------------------------------------------*
      subroutine get_mrcc_response_input(orb_info)
*----------------------------------------------------------------------*
*     process response functions input for ic-MRCC and write everything 
*     in a file.
*     Adapted from get_response_input and this now works only for one
*     perturbation, though it is written for any number of perturbtions.
*
*     Pradipta, 2016
*----------------------------------------------------------------------*
      implicit none

      include 'def_target.h'
      include 'def_orbinf.h'
      include 'ifc_input.h'
      include 'def_pert_info.h'
      include 'def_filinf.h'

      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffpropinf
      type(pert_op_info) ::
     &     pop(3*maxpop)
      type(pert_component_info) ::
     &     cmp(maxcmp)
      integer ::
     &     ncnt, icnt, idx, jdx, ipop, pos, sign
      integer ::
     &     ncmp, npop, luprop, idum, maxord, new_pos
      integer, allocatable :: 
     &     ord(:)
      integer, allocatable ::
     &     prop_comp(:,:), conj_comp(:,:), conj_prop(:,:)
      real ::
     &     order
      logical ::
     &     skip, trplt
      character(len=1) ::
     &     pert_ord
      character(len=6) ::
     &     int_name
      character(len_command_par) ::
     &     pert, pertop
      character(20),parameter::
     &     name_propinf="prop_info.gecco"

      integer, external ::
     &     pert_sym

      ncnt = is_keyword_set('method.MRCC.response')
      if (ncnt.eq.0) then
        call quit(1,'get_mrcc_response_input',
     &            'response keyword is not set')
      else if(ncnt.gt.1) then
!       call quit(1,'get_mrcc_response_input',
!    &            'response keyword can be set only once, for now')
      end if

      allocate(ord(ncnt))

      do icnt = 1,ncnt
        call get_argument_value('method.MRCC.response','order',
     &       keycount=icnt,ival=ord(icnt))
      end do

      call get_argument_value('calculate.properties','triplet',
     &     lval=trplt)

      maxord = maxval(ord(:))
      allocate(prop_comp(ncnt,maxord))
      allocate(conj_comp(ncnt,maxord))
      allocate(conj_prop(ncnt,maxord))
      prop_comp=0
      conj_comp=0
      conj_prop=0

      ncmp = ncnt*maxord
      cmp(:)%redun = 0
      cmp(:)%pop_idx = 0
      cmp(:)%freq = 0d0
      cmp(:)%order=0
      if (maxord.gt.maximum_order) then
        call quit(1,'get_mrcc_response_input',
     &        'ord must not exceed '//pert_ord)
      end if
      pert(1:len_command_par) = ' '
      pertop(1:len_command_par) = ' '

      npop = 0
      do icnt = 1,ncnt
        pos = (icnt-1)*maxord + 1
        call get_argument_value('method.MRCC.response','comp',
     &       keycount=icnt,str=pert(pos:len_command_par))
        call get_argument_value('method.MRCC.response','pert',
     &       keycount=icnt,str=pertop(pos:len_command_par))
        ! Getting the frequencies of all the perturbations and storing them
        ! starting from the second place in cmp%freq
        call get_argument_value('method.MRCC.response','freq',
     &       keycount=icnt,xarr=cmp(pos+1:ncmp)%freq)
        cmp(pos+ord(icnt):ncmp)%freq = 0d0
        ! Then putting the frequency at the first position of cmp
        ! This is the negative of the sum of all the frequencies
        if (ord(icnt).gt.0)
     &       cmp(pos)%freq = -sum(cmp(pos+1:pos+ord(icnt)-1)%freq)

        ! duplicate values for pert if necessary and not specified
        do idx = pos+1,pos+ord(icnt)-1
          if (pert(idx:idx).eq.' ')
     &        pert(idx:idx) = pert(idx-1:idx-1)
          if (pertop(idx:idx).eq.' ')
     &        pertop(idx:idx) = pertop(idx-1:idx-1)
        end do

        do idx = pos,pos+ord(icnt)-1
          ! check if perturbation operator input is ok
          if (pert(idx:idx).ne.'X' .and.
     &          pert(idx:idx).ne.'Y' .and.
     &          pert(idx:idx).ne.'Z')
     &          call quit(1,'get_mrcc_response_input',
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
          if (.not.skip) call quit(1,'get_mrcc_response_input',
     &          'perturbation operator "'//pertop(idx:idx)//
     &          '" is currently not allowed.')
          skip = .false.
          ! get the order corresponding to the component 
          order=((real(ord(icnt))-1.0d0)/2.0d0)
          cmp(idx)%order=ceiling(order)

          if(cmp(idx)%order.gt.0) then
            do ipop = 1,npop
              if (pop(ipop)%comp.eq.pert(idx:idx).and.
     &                pop(ipop)%name.eq.pertop(idx:idx)) then
                skip = .true.
                cmp(idx)%pop_idx = ipop
                exit
              end if
            end do
          end if

          prop_comp(icnt,idx-(icnt-1)*maxord) = idx

          if (.not.skip) then
            npop = npop + 1
            pop(npop)%comp = pert(idx:idx)
            pop(npop)%name = pertop(idx:idx)
            pop(npop)%int_name = pert(idx:idx)//int_name//' '
            pop(npop)%sign = sign
            if(trplt) then
                pop(npop)%isym = 1  ! If the perturbation is triplet,
                                    !then we can ignore the isym and set it to 1
            else
                pop(npop)%isym = pert_sym(pop(npop)%int_name,orb_info)
            end if
            cmp(idx)%pop_idx = npop
          end if

          ! determine redundancies
          do jdx = 1,idx-1
            if (pert(idx:idx).eq.pert(jdx:jdx) .and.
     &          pertop(idx:idx).eq.pertop(jdx:jdx) .and.
     &          cmp(idx)%order.eq.cmp(jdx)%order.and.
     &          abs(cmp(idx)%freq - cmp(jdx)%freq).lt.1d-12)
     &         cmp(idx)%redun = jdx
          end do
          if (cmp(idx)%redun.eq.0) cmp(idx)%redun = idx
        end do
      end do
      
      ! Here we get and store the information about each call  
      ! for the response calculation. 
      do icnt = 1, ncnt
        pos = (icnt-1)*maxord + 1
          if (ord(icnt).eq.1)then
!             prop_comp(icnt,1)=pos
              conj_comp(icnt,1)=pos
              conj_prop(icnt,1)=cmp(pos)%pop_idx
          elseif(ord(icnt).eq.2)then
            if (cmp(pos)%redun.eq.pos) then
              new_pos=cmp(pos+1)%redun
              conj_comp(icnt,1)=new_pos
              conj_prop(icnt,1)=cmp(new_pos)%pop_idx
            else
              new_pos=cmp(pos)%redun
              prop_comp(icnt,1)=new_pos
              conj_comp(icnt,1)=new_pos
              conj_prop(icnt,1)=cmp(new_pos)%pop_idx
            end if
            if (cmp(pos+1)%redun.eq.pos+1) then
              new_pos=cmp(pos)%redun
              conj_comp(icnt,2)=new_pos
              conj_prop(icnt,2)=cmp(new_pos)%pop_idx
            else
              new_pos=cmp(pos+1)%redun
              prop_comp(icnt,2)=new_pos
              conj_comp(icnt,2)=new_pos
              conj_prop(icnt,2)=cmp(new_pos)%pop_idx
            end if
          end if
      end do

      call file_init(ffpropinf,trim(name_propinf),ftyp_sq_frm,idum)
      call file_open(ffpropinf)
      luprop = ffpropinf%unit

      write(luprop,*) 'npop', npop
      do ipop=1,npop
        write(luprop,*) 'pop(npop)%comp:  ', pop(ipop)%comp
        write(luprop,*) 'pop(npop)%name:  ', pop(ipop)%name
        write(luprop,*) 'pop(npop)%int_name:  ', pop(ipop)%int_name
        write(luprop,*) 'pop(npop)%sign:  ', pop(ipop)%sign
        write(luprop,*) 'pop(npop)%isym:  ', pop(ipop)%isym
      enddo
      write(luprop,*) ncnt, maxord
      do idx=1,ncnt*maxord
        write(luprop,*) 'cmp(idx)%pop_idx', cmp(idx)%pop_idx
        write(luprop,*) 'cmp(idx)%freq', cmp(idx)%freq
        write(luprop,*) 'cmp(idx)%redun', cmp(idx)%redun
        write(luprop,*) 'cmp(idx)%order', cmp(idx)%order
      end do
      write(luprop,*) 'order: ', ord(1:ncnt)
      do icnt=1,ncnt
        write(luprop,*) 'prop_comp:', prop_comp(icnt,:)
        write(luprop,*) 'conj_comp:', conj_comp(icnt,:)
        write(luprop,*) 'conj_prop:', conj_prop(icnt,:)
      end do

      deallocate(ord)

      return
      end
