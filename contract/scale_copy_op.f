!> Subroutine to scale matrix element lists
!!
!!This suroutine has mutiple uses depending on mode. 
!!+ "":  scales all matrix elements by fac.
!!+ "mult": multiplies the MELs elementwise and ap
!!+ "square": squares all matrix elements and multiplies them by fac.
!!+ "precond"": divides label_res elementwise by label_inp and applies fac(1).
!!+ "prc_thresh": sets all elements to at least fac.
!!+ any other:  scales all matrix elements by fac.
!!
!!\param[in] mode one of: "square", "mult", "precond", "prc_thresh" or empty
!!\param[in] fac list of factors which are applied (repeatedly)
!!\param[in] op_info
!!\param[out] label_res
      subroutine scale_copy_op(label_res,label_inp,fac,nfac,mode,nspc,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00
      character(len=13), parameter ::
     &     i_am = 'scale_copy_op'
      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'
      include 'def_optimize_info.h'
      include 'par_scale_copy_modes.h'


      
      interface
         subroutine mel_scale_copy(
     &     me_inp, me_res,
     &     buf1,buf2, lbuf, nbuf,
     &     fac, nfac,
     &     mode,
     &        opti_info)
         import optimize_info,me_list
         integer,intent(in)::
     &        lbuf, nbuf, nfac,
     &        mode
         type(me_list), intent(inout)::
     &        me_inp, me_res
         type(optimize_info),intent(in)::
     &        opti_info
         real(8), intent(inout) ::
     &        buf1(:), buf2(:)
         real(8),intent(in)::
     &        fac(nfac)
         end subroutine
      end interface

      
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
     &     idx_res, idx_inp, idx, idxnd_src, idxnd_tgt,  nbuff, lbuff,
     &     ipri, idxst_tgt, idxst_src,
     &     ifac, nblkmax, ifree, nblk, idx_shape,
     &     imode, len_op, isec
      logical ::
     &     open_close_res, open_close_inp,
     &     same, warning
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      real(8), pointer ::
     &     buffer(:), buf_in(:), signsec(:)

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

      
      select case(trim(mode))
      case('square')
         imode=MOD_SQUARE
      case('mult')
         imode=MOD_MULT
      case('precond')
         imode=MOD_PRECOND
         if (nfac.ne.1)
     &        call quit(1,i_am,
     &        'for mode precond, only a single factor must be given')
      case('prc_thresh')
         imode=MOD_PRC_THRESH
      case default
         imode = MOD_SCALE
      end select
      
      idx_res = idx_mel_list(label_res,op_info)
      idx_inp = idx_mel_list(label_inp(1),op_info)
      idx_shape = -1
      if (nspc.gt.0) idx_shape = idx_mel_list(label_inp(2),op_info)

      if (idx_res.lt.0) then
        write(lulog,*) '"',trim(label_res),'"'
        write(lulog,*) idx_res
        call quit(1,i_am,'result label not on list (1)')
      end if
      if (idx_inp.lt.0) then
        write(lulog,*) '"',trim(label_inp(1)),'"'
        write(lulog,*) idx_inp
        call quit(1,i_am,'input label not on list')
      end if




      ! Point to the relevant operators and their associated files.
      same = idx_res.eq.idx_inp
! result List
      
      me_res => op_info%mel_arr(idx_res)%mel
      ffop_tgt => me_res%fhand
      if (.not.associated(ffop_tgt))
     &     call quit(1,i_am,'no file handle defined for '//
     &                  trim(me_res%label))

      open_close_res = ffop_tgt%unit.le.0
      if(open_close_res)then
        call file_open(ffop_tgt)
      endif

      ! input List
      me_inp => op_info%mel_arr(idx_inp)%mel
      ffop_src => me_inp%fhand
      if (.not.same) then
        if (.not.associated(ffop_src))
     &     call quit(1,i_am,'no file handle defined for '//
     &                  trim(me_inp%label))
        open_close_inp = ffop_src%unit.le.0
        if (open_close_inp) then
          call file_open(ffop_src)
        endif
      else
        open_close_inp = .false.
      end if

      


! if one list is shorter, we will just end copying process there
      
      ! is there a sign correction (due to formal contraction)?
      if (idx_shape.ge.0) then
        allocate(me_vec(1),me_shape(1))
        me_vec(1)%mel => me_res
        me_shape(1)%mel => op_info%mel_arr(idx_shape)%mel
        ! put sign corrections on opti_info
        call set_opti_info_signs(opti_info,1,1,
     &            me_vec,me_shape,me_shape,me_shape,.false.)
        deallocate(me_vec,me_shape)
        print *,"left opti_info signs"
      else
        ifree = mem_alloc_int(opti_info%nsec,1,'nsec')
        ifree = mem_alloc_int(opti_info%nwfpsec,1,'nwfpsec')
        ifree = mem_alloc_int(opti_info%idstsec,1,'idstsec')
        ifree = mem_alloc_real(opti_info%signsec,1,'signsec')
        opti_info%nsec(1) = 1
        opti_info%nwfpsec(1) = min(me_inp%len_op,me_res%len_op)
        opti_info%idstsec(1) = 1
        opti_info%signsec(1) = 1d0
      end if


      ! record length hopefully the same
      if (ffop_tgt%reclen.ne.ffop_src%reclen)
     &   call quit(1,i_am,
     &     'not prepared for different reclen''s')
      
      nblkmax = ifree/ffop_src%reclen
      if (trim(mode).eq.'mult'.or.trim(mode).eq.'precond')
     &      nblkmax = nblkmax/2
      if (nblkmax.le.0) then
        write(lulog,*) 'free memory (words):  ',ifree
        write(lulog,*) 'block length (words): ',ffop_src%reclen
        call quit(1,i_am,
     &            'not even 1 record fits into memory?')
      end if

      
      ifree = mem_setmark('scale_copy_section')
      
      len_op =0
      do isec = 1,  opti_info%nsec(1)
         len_op = max(len_op,opti_info%nwfpsec(isec))
      end do

      if (.not.ffop_src%buffered .and.
     &     .not.ffop_tgt%buffered) then

         nblk = min((len_op-1)/ffop_src%reclen + 1,nblkmax)

         lbuff = min(len_op,nblk*ffop_src%reclen)

         ifree = mem_alloc_real(buffer,lbuff,'buffer')
         
         if (trim(mode).eq.'mult'.or.trim(mode).eq.'precond')then
            ifree = mem_alloc_real(buf_in,lbuff,'buf_in')
            nbuff=2
         else
            buf_in => null()
            nbuff=1
         end if


         
          call mel_scale_copy(
     &     me_inp, me_res,
     &     buffer,buf_in, lbuff, nbuff ,
     &     fac, nfac,
     &     imode,
     &     opti_info)
         
        
         ifree = mem_flushmark()
      else

         call quit(1,i_am,'adapt for buffering')
         
      end if
        
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
      end subroutine
      

