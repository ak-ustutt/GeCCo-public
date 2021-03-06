*----------------------------------------------------------------------*
      subroutine dia_from_op(label_out,label_inp,mode,op_info,
     &                      str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for dia_from_blk
*     finds and resets output blocks that correspond to the diagonal of
*     the input operator and does not alter any other non-matching
*     operator blocks
*
*     matthias, fall 2009
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'
      include 'ifc_operators.h'

      character(*), intent(in) ::
     &     label_out, label_inp, mode
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     first, ms_fix, fix_success, open_close_inp,
     &     open_close_out, extend, ext_act, zero_dia, set_dia
      integer ::
     &     ifree, idxout, idxinp, iblkinp, iblkout, i_occ_cls,
     &     iocc_dia(ngastp,2), njoined, split_sign,
     &     lenblkinp, lenblkout, ioffinp, ioffout
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      real(8), pointer ::
     &     buffer_inp(:), buffer_out(:)

      type(me_list), pointer ::
     &     meinp, meout
      type(filinf), pointer ::
     &     ffinp, ffout
      type(operator), pointer ::
     &     opinp, opout
      integer, pointer ::
     &     iocc_inp(:,:,:), iocc_out(:,:)
      integer, external ::
     &     idx_mel_list, ndisblk_mel, iblk_occ
      logical, external ::
     &     occ_is_diag_blk

      call atim_csw(cpu0,sys0,wall0)

      extend = .false.
      ext_act = .false.
      set_dia = .false.
      zero_dia = .false.
      select case(trim(mode))
      case('extend')
        extend = .true.
      case('ext_act')
        extend = .true.
        ext_act = .true.
      case('set_dia')
        extend = .true.
        set_dia = .true.
      case('zero_dia')
        extend = .true.
        zero_dia = .true.
        set_dia = .true.
      end select

      idxout = idx_mel_list(label_out,op_info)
      if (idxout.lt.0) call quit(1,'dia_from_op',
     &                           'label not on list: '//trim(label_out))
      meout => op_info%mel_arr(idxout)%mel
      ffout => meout%fhand
      opout => meout%op
      if (.not.associated(ffout))
     &     call quit(1,'dia_from_op','no file handle defined for '//
     &                  trim(meout%label))
      open_close_out = ffout%unit.le.0

      idxinp = idx_mel_list(label_inp,op_info)
      if (idxinp.lt.0) then
        call quit(1,'dia_from_op',
     &       'label not on list: '//trim(label_inp))
      end if
      meinp => op_info%mel_arr(idxinp)%mel
      ffinp => meinp%fhand
      opinp => meinp%op
      if (.not.associated(ffinp))
     &     call quit(1,'dia_from_op','no file handle defined for '//
     &     trim(meinp%label))
      open_close_inp = ffinp%unit.le.0

      if (open_close_out)
     &     call file_open(ffout)
      if (open_close_inp)
     &     call file_open(ffinp)

      if (ntest.ge.100) then
        write(lulog,*) '============================'
        write(lulog,*) ' dia_from_op messing around'
        write(lulog,*) '============================'
        write(lulog,*) ' input list = ',trim(meinp%label)
        write(lulog,*) ' ffinp: ',trim(ffinp%name)
        write(lulog,*) ' opinp: ',opinp%name(1:len_trim(opinp%name))
        write(lulog,*) ' output list = ',trim(meout%label)
        write(lulog,*) ' ffout: ',trim(ffout%name)
        write(lulog,*) ' opout: ',opout%name(1:len_trim(opout%name))
        write(lulog,*) ' mode: ',trim(mode)
      end if

      if (opout%njoined.ne.1.or.opinp%njoined.gt.3)
     &    call quit(1,'dia_from_op','not made for these operator types')
      njoined = opinp%njoined

      ! loop over occupation classes of operator representing the full matrix
      do i_occ_cls = 1, opinp%n_occ_cls
        iblkinp = (i_occ_cls-1)*njoined+1
        iocc_inp => opinp%ihpvca_occ(1:ngastp,1:2,
     &                                iblkinp:iblkinp+njoined-1)
        ! only for diagonal blocks
        if (.not.
     &       occ_is_diag_blk(iocc_inp(1:ngastp,1:2,1:njoined),njoined))
     &       cycle

        ! find occupation of diagonal. possibilities are:
        ! /a 0 0 0\                                         /e b c d\
        ! \0 f g h/                /a 0 0 0\                \e b c d/
        ! /e b c d\ --> /a b c d\  \0 f g h/ --> /a 0 0 0\      |
        ! \e b c d/     \e f g h/; /0 f g h\     \0 f g h/;    \|/
        ! /0 f g h\                \a 0 0 0/                /0 b c d\
        ! \a 0 0 0/                                         \e 0 0 0/
        iocc_dia = 0
        split_sign = 1
        if (njoined.eq.1) then
          split_sign = 1-2*mod(iocc_inp(1,1,1),2)
          iocc_dia(1,2) = iocc_inp(1,2,1)
          iocc_dia(2:ngastp,1) = iocc_inp(2:ngastp,1,1)
        else if (njoined.eq.2) then
          iocc_dia(1,1) = iocc_inp(1,1,1)
          iocc_dia(2:ngastp,2) = iocc_inp(2:ngastp,2,1)
        else if (njoined.eq.3) then
          split_sign = 1-2*mod(iocc_inp(1,1,2),2)
          iocc_dia(1,1) = iocc_inp(1,1,1)
          iocc_dia(1,2) = iocc_inp(1,1,2)
          iocc_dia(2:ngastp,1) = iocc_inp(2:ngastp,1,2)
          iocc_dia(2:ngastp,2) = iocc_inp(2:ngastp,2,1)
        else
          call quit(1,'dia_from_op','not made for these operator types')
        end if
        if (ntest.ge.100) then
          write(lulog,'(a,i4)') 'diagonal of occupation class',i_occ_cls
          call wrt_occ(lulog,iocc_dia)
          write(lulog,*) 'split_sign: ',split_sign
        end if

        ! is there a matching output block?
        do iblkout = 1, opout%n_occ_cls
          iocc_out => opout%ihpvca_occ(1:ngastp,1:2,iblkout)
          if (iocc_equal(iocc_dia,.false.,iocc_out,.false.).or.
     &        extend.and. ! allow additional lines of different types only
     &        iocc_bound('<=',iocc_dia,.false.,iocc_out,.false.).and.
     &        iocc_zero(iocc_overlap(iocc_add(-1,iocc_dia,.false.,
     &            1,iocc_out,.false.),.false.,iocc_dia,.false.)).and.
     &        iocc_zero(iocc_overlap(iocc_add(-1,iocc_dia,.false.,
     &            1,iocc_out,.false.),.false.,iocc_dia,.true.)).and.
c           ad hoc: add scalar contrib. only to purely inactive blks
     &        (ext_act.or..not.(iocc_zero(iocc_dia).and.
     &              any(iocc_out(IVALE,1:2).ne.0)))) then
            if (ntest.ge.100) write(lulog,'(a,i2)')
     &           'found matching output block: # ',iblkout

            ! get the input, extract diagonal and write to output
            ifree = mem_setmark('dia_from_op')
            lenblkinp = meinp%len_op_occ(i_occ_cls)
            lenblkout = meout%len_op_occ(iblkout)
            ifree = mem_alloc_real(buffer_inp,lenblkinp,'buffer_inp')
            ifree = mem_alloc_real(buffer_out,lenblkout,'buffer_out') 
            ioffinp = meinp%off_op_occ(i_occ_cls)
     &       +ffinp%length_of_record*(ffinp%current_record-1)
            ioffout = meout%off_op_occ(iblkout)
     &       +ffout%length_of_record*(ffout%current_record-1)
            call get_vec(ffinp,buffer_inp,ioffinp+1,ioffinp+lenblkinp)
            if (extend.and..not.zero_dia) then ! no reset of output, only add!
              call get_vec(ffout,buffer_out,ioffout+1,ioffout+lenblkout)
            else ! reset output buffer
              buffer_out(1:lenblkout) = 0d0
            end if
            call dia_from_blk(buffer_out,buffer_inp,set_dia,
     &                        meinp,meout,i_occ_cls,iblkout,
     &                        iocc_dia,str_info,orb_info)
            if (set_dia) then
              call put_vec(ffinp,buffer_inp,ioffinp+1,ioffinp+lenblkinp)
            else
              call put_vec(ffout,buffer_out,ioffout+1,ioffout+lenblkout)
            end if
            ifree = mem_flushmark('dia_from_op')

            if (.not.extend) exit ! not more than one matching block
          end if
        end do
      end do

      if (open_close_inp)
     &     call file_close_keep(meinp%fhand)
      if (open_close_out)
     &     call file_close_keep(meout%fhand)

      call atim_csw(cpu,sys,wall)

      call prtim(lulog,'time for extracting diagonal ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
