*----------------------------------------------------------------------*
      subroutine set_prc4op(label_prc,mode_str,
     &                      label_inp,nlabel_inp,
     &                      op_info,
     &                      str_info,orb_info)
*----------------------------------------------------------------------*
*     driver routine for setting up a diagonal preconditioner
*
*     mode_str: '' or 'dia-F'  -- the usual choice: diagonal using F
*               'dia-R12'      -- diagonal precond. for R12
*               'dia-R12noX'   -- R12: but do not use X
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 000

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, intent(in) ::
     &     nlabel_inp
      character(*), intent(in) ::
     &     mode_str,
     &     label_prc, label_inp(nlabel_inp)
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     use_h, use_b1, use_x1, use_b, use_x, too_few,
     &     open_close_b,
     &     open_close_x,
     &     open_close_ham,
     &     open_close_prc
      integer ::
     &     ifree, nops,
     &     idxprc, idxham, idx_b, idx_x, x2nblk, b2nblk,
     &     occ_test(ngastp,2)
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0

      type(me_list), pointer ::
     &     me_ham, me_prc, me_b, me_x

      real(8), target ::
     &     xdummy(1)

      integer, pointer ::
     &     b2off(:), x2off(:)
      real(8), pointer ::
     &     h1dia(:), b1dia(:), x1dia(:), b2dia(:), x2dia(:)

      integer, external ::
     &     idx_mel_list, ndisblk_mel, iblk_occ

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('prc4op')

      too_few = .false.
      select case(trim(mode_str))
      case('','dia-F')
        use_h = .true.
        use_b = .false.
        use_x = .false.
        if (nlabel_inp.lt.1) too_few = .true.
      case('dia-R12')
        use_h = .true.
        use_b = .true.
        use_x = .true.
        if (nlabel_inp.lt.3) too_few = .true.
      case('dia-R12noX')
        use_h = .true.
        use_b = .true.
        use_x = .false.
        if (nlabel_inp.lt.2) too_few = .true.
      case default
        call quit(1,'set_prc4op','unknown mode: "'//trim(mode_str)//'"')
      end select
      if (too_few)
     &     call quit(1,'set_prc4op',
     &     'Too few labels for mode: "'//trim(mode_str)//'"')

      idxprc = idx_mel_list(label_prc,op_info)
      if (idxprc.lt.0) then
        call quit(1,'set_prc4op','label not on list: '//trim(label_prc))
      end if
      me_prc => op_info%mel_arr(idxprc)%mel
      if (.not.associated(me_prc%fhand))
     &     call quit(1,'set_prc4op','no file handle defined for '//
     &                  trim(me_prc%label))
      open_close_prc = me_prc%fhand%unit.le.0

      if (use_h) then
        idxham = idx_mel_list(label_inp(1),op_info)
        if (idxham.lt.0) then
          call quit(1,'set_prc4op',
     &         'label not on list: '//trim(label_inp(1)))
        end if
        me_ham => op_info%mel_arr(idxham)%mel
        if (.not.associated(me_ham%fhand))
     &       call quit(1,'set_prc4op','no file handle defined for '//
     &       trim(me_ham%label))
        open_close_ham = me_ham%fhand%unit.le.0
      end if

      if (use_b) then
        idx_b = idx_mel_list(label_inp(2),op_info)
        if (idx_b.lt.0) then
          call quit(1,'set_prc4op',
     &         'label not on list: '//trim(label_inp(2)))
        end if
        me_b => op_info%mel_arr(idx_b)%mel
        if (.not.associated(me_b%fhand))
     &       call quit(1,'set_prc4op','no file handle defined for '//
     &       trim(me_b%label))
        open_close_b = me_ham%fhand%unit.le.0
        occ_test = 0
        occ_test(IPART,1:2) = 1
        use_b1 = iblk_occ(occ_test,.false.,me_b%op).gt.0
      end if

      if (use_x) then
        idx_x = idx_mel_list(label_inp(3),op_info)
        if (idx_x.lt.0) then
          call quit(1,'set_prc4op',
     &         'label not on list: '//trim(label_inp(3)))
        end if
        me_x => op_info%mel_arr(idx_x)%mel
        if (.not.associated(me_x%fhand))
     &       call quit(1,'set_prc4op','no file handle defined for '//
     &       trim(me_x%label))
        open_close_x = me_ham%fhand%unit.le.0
        occ_test = 0
        occ_test(IPART,1:2) = 1
        use_x1 = iblk_occ(occ_test,.false.,me_x%op).gt.0
      end if

      h1dia => xdummy
      b2dia => xdummy
      x2dia => xdummy
      ! this assumption is probably not too bad:
      if (use_h)
     &     ifree = mem_alloc_real(h1dia,2*orb_info%ntoob,'h1dia')
      if (use_b1)
     &     ifree = mem_alloc_real(b1dia,2*orb_info%ntoob,'b1dia')
      if (use_x1)
     &     ifree = mem_alloc_real(x1dia,2*orb_info%ntoob,'x1dia')
      if (use_b) then        
        ifree = mem_alloc_real(b2dia,4*orb_info%ntoob**2,'b2dia')
        b2nblk = ndisblk_mel(me_b)    * 10
        ifree = mem_alloc_int (b2off,b2nblk,'b2off')
      end if
      if (use_x) then
        ifree = mem_alloc_real(x2dia,4*orb_info%ntoob**2,'x2dia')
        x2nblk = ndisblk_mel(me_x)    * 10
        ifree = mem_alloc_int (x2off,x2nblk,'x2off')
      end if

      if (iprlvl.ge.1)
     &     write(luout,*) 'set up diagonal'//
     &     ' from rank 1 part of ',trim(me_ham%op%name)
      if (iprlvl.ge.1.and.use_b.and.use_x)
     &     write(luout,*) 
     &     '     and the diagonals of ',trim(me_b%op%name),
     &                          ' and ',trim(me_x%op%name)
      if (iprlvl.ge.1.and.use_b.and..not.use_x)
     &     write(luout,*) 
     &     '     and the diagonal of ',trim(me_b%op%name)

      if (open_close_prc)
     &     call file_open(me_prc%fhand)
      if (use_h.and.open_close_ham)
     &     call file_open(me_ham%fhand)
      if (use_b.and.open_close_b)
     &     call file_open(me_b%fhand)
      if (use_x.and.open_close_x)
     &     call file_open(me_x%fhand)

      ! extract the fock-matrix diagonal
      if (use_h)
     &     call onedia_from_op(h1dia,me_ham,orb_info)
      ! diagonal of partial trace of B/X (for R12):
      if (use_b1)
     &     call onedia_from_op(b1dia,me_b,orb_info)
      if (use_x1)
     &     call onedia_from_op(x1dia,me_x,orb_info)


      ! Extract the diagonal elements of the B-matrix for R12.
      if (use_b)
     &     call twodia_from_op(b2dia,!b2off,b2nblk,
     &                         me_b,
     &                         orb_info,str_info)
c dbg
c      if (use_b1) stop 'testing'
c dbg
      ! Extract the diagonal elements of the X-matrix for R12.
      if (use_x)
     &     call twodia_from_op(x2dia,!x2off,x2nblk,
     &                         me_x,
     &                         orb_info,str_info)

      ! set up preconditioner
      if (.not.use_b.and..not.use_x) then
        call dia4op(me_prc,h1dia,str_info,orb_info)
      else
        call dia4op_r12(me_prc,h1dia,b2dia,b2off,x2dia,x2off,use_x,
     &       str_info,orb_info)
      end if

      if (ntest.ge.1000) then
        call wrt_mel_file(luout,5,me_prc,1,
     &       me_prc%op%n_occ_cls,str_info,orb_info)
      end if

      if (open_close_ham)
     &     call file_close_keep(me_ham%fhand)
      if (open_close_prc)
     &     call file_close_keep(me_prc%fhand)

      ifree = mem_flushmark()

      call atim_csw(cpu,sys,wall)

      call prtim(luout,'time for diagonal ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
