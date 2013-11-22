      subroutine scale_copy_op(label_res,label_inp,fac,nfac,mode,nspc,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00
      
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'
      include 'def_optimize_info.h'

      integer, intent(in) ::
     &     nfac, nspc
      real(8), intent(in) ::
     &     fac(nfac)
      character(*), intent(in) ::
     &     label_res, label_inp(1)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info
      character(len=*), intent(in) ::
     &     mode

      type(me_list), pointer ::
     &     me_res, me_inp
      type(me_list_array), pointer ::
     &     me_vec(:), me_shape(:)
      type(filinf), pointer ::
     &     ffop_src, ffop_tgt
      type(optimize_info) ::
     &     opti_info

      integer ::
     &     idx_res, idx_inp, idx, idxnd_src, idxnd_tgt, len_op, nbuff,
     &     ipri, idxst_tgt, idxst_src, idisc_off_src, idisc_off_tgt,
     &     ifac, nblkmax, ifree, nblk, idx_shape, isec, nsec,
     &     smapre_num, negpre_num
      logical ::
     &     open_close_res, open_close_inp,
     &     same, warning
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      real(8), pointer ::
     &     buffer(:), buf_in(:), signsec(:)
      integer, pointer ::
     &     lensec(:), idstsec(:)

      integer, external ::
     &     idx_mel_list

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('scale_copy_op')

      if(ntest.ge.100)then
        write(lulog,*) '===================='
        write(lulog,*) ' scale & copy       '
        write(lulog,*) '===================='
        write(lulog,*) 'Result: ',trim(label_res)
        write(lulog,*) 'Input:  ',trim(label_inp(1))
        write(lulog,*) 'mode:   ',trim(mode)
        write(lulog,*) 'The factors (applied periodically): '
        do idx = 1, nfac
          write(lulog,'(3x,f12.6)') fac(idx)
        end do
      endif

      if (trim(mode).eq.'precond'.and.nfac.ne.1)
     &    call quit(1,'scale_copy_op',
     &    'for mode precond, only a single factor must be given')

      idx_res = idx_mel_list(label_res,op_info)
      idx_inp = idx_mel_list(label_inp(1),op_info)
      idx_shape = -1
      if (nspc.gt.0) idx_shape = idx_mel_list(label_inp(2),op_info)

      if (idx_res.lt.0) then
        write(lulog,*) '"',trim(label_res),'"'
        write(lulog,*) idx_res
        call quit(1,'scale_copy_op','label not on list (1)')
      end if
      if (idx_inp.lt.0) then
        write(lulog,*) '"',trim(label_inp(1)),'"'
        write(lulog,*) idx_inp
        call quit(1,'scale_copy_op','label not on list (2)')
      end if

      same = idx_res.eq.idx_inp
c      if (same) call quit(1,'scale_copy_op',
c     &       'copy means to have a second list...')

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      ffop_tgt => me_res%fhand
      if (.not.associated(ffop_tgt))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = ffop_tgt%unit.le.0

      if(open_close_res)then
        call file_open(ffop_tgt)
      endif

      me_inp => op_info%mel_arr(idx_inp)%mel
      ffop_src => me_inp%fhand
      if (.not.same) then
        if (.not.associated(ffop_src))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_inp%label))
        open_close_inp = ffop_src%unit.le.0
        if (open_close_inp) then
          call file_open(ffop_src)
        endif
      else
        open_close_inp = .false.
      end if

      ! if one list is shorter, we will just end copying process there
      len_op = min(me_inp%len_op,me_res%len_op)

      ! offset on file (if more than one instance of list exists)
      idisc_off_src =
     &     ffop_src%length_of_record*(ffop_src%current_record-1)
      idisc_off_tgt =
     &     ffop_tgt%length_of_record*(ffop_tgt%current_record-1)

      ! is there a sign correction (due to formal contraction)?
      if (idx_shape.ge.0) then
        allocate(me_vec(1),me_shape(1))
        me_vec(1)%mel => me_res
        me_shape(1)%mel => op_info%mel_arr(idx_shape)%mel
        ! put sign corrections on opti_info
        call set_opti_info_signs(opti_info,1,1,
     &            me_vec,me_shape,me_shape,me_shape,.false.)
        deallocate(me_vec,me_shape)
      else
        ifree = mem_alloc_int(opti_info%nsec,1,'nsec')
        ifree = mem_alloc_int(opti_info%nwfpsec,1,'nwfpsec')
        ifree = mem_alloc_int(opti_info%idstsec,1,'idstsec')
        ifree = mem_alloc_real(opti_info%signsec,1,'signsec')
        opti_info%nsec(1) = 1
        opti_info%nwfpsec(1) = len_op
        opti_info%idstsec(1) = 1
        opti_info%signsec(1) = 1d0
      end if
      nsec = opti_info%nsec(1)
      lensec => opti_info%nwfpsec(1:nsec)
      idstsec => opti_info%idstsec(1:nsec)
      signsec => opti_info%signsec(1:nsec)
      warning = .false.

      ! record length hopefully the same
      if (ffop_tgt%reclen.ne.ffop_src%reclen)
     &   call quit(1,'scale_copy_op',
     &             'not prepared for different reclen''s')
      nblkmax = ifree/ffop_src%reclen
      if (trim(mode).eq.'mult'.or.trim(mode).eq.'precond')
     &      nblkmax = nblkmax/2
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop_src%reclen
        call quit(1,'scale_copy_op',
     &            'not even 1 record fits into memory?')
      end if

      ! loop over sections
      do isec = 1, nsec
      ifree = mem_setmark('scale_copy_section')
      len_op = lensec(isec)

      if (.not.ffop_src%buffered.and.
     &    .not.ffop_tgt%buffered) then

        nblk = min((len_op-1)/ffop_src%reclen + 1,nblkmax)

        nbuff = min(len_op,nblk*ffop_src%reclen)

        ifree = mem_alloc_real(buffer,nbuff,'buffer')
        if (trim(mode).eq.'mult'.or.trim(mode).eq.'precond') 
     &          ifree = mem_alloc_real(buf_in,nbuff,'buf_in')

        ifac = 1
        idxst_src = idisc_off_src+idstsec(isec)
        idxst_tgt = idisc_off_tgt+idstsec(isec)
        negpre_num = 0
        smapre_num = 0
        do while(idxst_src.le.idisc_off_src+idstsec(isec)-1+len_op)
          idxnd_src = min(idisc_off_src+idstsec(isec)-1+len_op,
     &                    idxst_src-1+nbuff)
          idxnd_tgt = min(idisc_off_tgt+idstsec(isec)-1+len_op,
     &                    idxst_tgt-1+nbuff)
          call get_vec(ffop_src,buffer,idxst_src,idxnd_src)

          select case(trim(mode))
          case('square')
            ! take square and apply scaling factors
            do idx = 1, idxnd_src-idxst_src+1
              buffer(idx) = signsec(isec)*fac(ifac)*(buffer(idx)**2)
              ifac = ifac + 1
              if (ifac.gt.nfac) ifac = 1
            end do
          case('mult')
            ! multiply both lists element-wise
            call get_vec(ffop_tgt,buf_in,idxst_tgt,idxnd_tgt)
            do idx = 1, idxnd_src-idxst_src+1
              buffer(idx) = signsec(isec)*
     &                      fac(ifac)*buffer(idx)*buf_in(idx)
              ifac = ifac + 1
              if (ifac.gt.nfac) ifac = 1
            end do
          case('precond')
            ! divide lists element-wise
            ! here, buf_in contains nominator !!!
            call get_vec(ffop_tgt,buf_in,idxst_tgt,idxnd_tgt)
            call diavc(buffer,buf_in,
     &           signsec(isec)*fac(1),buffer,
     &           0d0,idxnd_src-idxst_src+1)
c          case('atleast')
c            ! elements should be at least the given number(s)
c            do idx = 1, idxnd_src-idxst_src+1
c              buffer(idx) = max(buffer(idx),fac(ifac))
c              ifac = ifac + 1
c              if (ifac.gt.nfac) ifac = 1
c            end do
c          case('atleastwarn')
c            ! warn if any element is below the given number(s)
c            do idx = 1, idxnd_src-idxst_src+1
c              if (buffer(idx).lt.fac(ifac)) then
c                warning = .true.
c                write(lulog,'(a,i14,a,f20.12,a,f20.12)')
c     &               'Element',idx,' with value',buffer(idx),
c     &               ' is below',fac(ifac)
c              end if
c              ifac = ifac + 1
c              if (ifac.gt.nfac) ifac = 1
c            end do
          case('prc_thresh')
            ! elements should be at least the given number(s)
            ! trigger a warning if any element was negative
            do idx = 1, idxnd_src-idxst_src+1
              if (buffer(idx).lt.fac(ifac)) then
                smapre_num = smapre_num + 1
                if (buffer(idx).lt.0d0) negpre_num = negpre_num + 1
                buffer(idx) = fac(ifac)
              end if
              ifac = ifac + 1
              if (ifac.gt.nfac) ifac = 1
            end do
          case default
            ! apply scaling factors (periodically)
            do idx = 1, idxnd_src-idxst_src+1
              buffer(idx) = signsec(isec)*fac(ifac)*buffer(idx)
c dbg
c              if (buffer(idx).lt.-1d-14) then
c                write(lulog,*) 'changing sign for el.# ',idx
c                buffer(idx) = abs(buffer(idx))
c              end if
c dbgend
              ifac = ifac + 1
              if (ifac.gt.nfac) ifac = 1
            end do
          end select

          call put_vec(ffop_tgt,buffer,idxst_tgt,idxnd_tgt)
          idxst_src = idxnd_src+1
          idxst_tgt = idxnd_tgt+1
        end do

        if (trim(mode).eq.'prc_thresh') then
          if (smapre_num.gt.0) then
            write(lulog,'(1x,a,i9,a,g9.2)')
     &           'number of small preconditioner elements: ',
     &           smapre_num, '; set to ',fac(1:nfac)
          end if
          if (negpre_num.gt.0) then
            write(lulog,'(1x,a,i9)')
     &           'number of negative preconditioner elements: ',
     &           negpre_num
            call warn('scale_copy_op',
     &                'negative preconditioner elements!')
          end if
        end if

      else

        call quit(1,'scale_copy_op','adapt for buffering')

      end if

      ifree = mem_flushmark()
      end do

      call touch_file_rec(ffop_tgt)

      ! needed: close loop over active records

      if (open_close_res)
     &     call file_close_keep(ffop_tgt)

      if (open_close_inp)
     &     call file_close_keep(ffop_src)

      if (ntest.ge.10) then
        write(lulog,*) 'dump of scaled list:'
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_mel_file(lulog,ipri,me_res,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      if (warning) call warn('scale_copy_op',
     &          'At least one element is below the recommended value.')

      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'time for scaling',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
