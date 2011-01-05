*----------------------------------------------------------------------*
      subroutine reo_mel(label_out,label_inp,op_info,
     &                      str_info,strmap_info,orb_info,fromto,dag)
*----------------------------------------------------------------------*
*     wrapper for reo_mel_blk
*     reorder vertices of input op according to instruction on fromto
*     e.g. fromto = 13 --> reorder vertex 1 to 3
*     result operator list must match for every occupation class
*
*     matthias, dec 2009 (adopted from dia_from_op)
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
      include 'def_strmapinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      character(*), intent(in) ::
     &     label_out, label_inp
      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     fromto
      logical, intent(in) ::
     &     dag

      logical ::
     &     open_close_inp,
     &     open_close_out, reorder
      integer ::
     &     ifree, idxout, idxinp, iblkinp, iblkout, i_occ_cls,
     &     njinp, njout, 
     &     lenblkinp, lenblkout, ioffinp, ioffout,
     &     ifrom, ito, ica, ij, j_occ_cls
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0, xdum, tra_sign
      real(8), pointer ::
     &     buffer_inp(:), buffer_out(:)

      type(me_list), pointer ::
     &     meinp, meout
      type(filinf), pointer ::
     &     ffinp, ffout
      type(operator), pointer ::
     &     opinp, opout
      integer, pointer ::
     &     iocc_inp(:,:,:), iocc_out(:,:,:)
      logical, pointer ::
     &     transposed(:)
      integer, external ::
     &     idx_mel_list
      logical, external ::
     &     iocc_zero, iocc_equal_n

      call atim_csw(cpu0,sys0,wall0)

      ifree = mem_setmark('reo_mel')

      idxout = idx_mel_list(label_out,op_info)
      if (idxout.lt.0) call quit(1,'reo_mel',
     &                           'label not on list: '//trim(label_out))
      meout => op_info%mel_arr(idxout)%mel
      ffout => meout%fhand
      opout => meout%op
      if (.not.associated(ffout))
     &     call quit(1,'reo_mel','no file handle defined for '//
     &                  trim(meout%label))
      open_close_out = ffout%unit.le.0

      idxinp = idx_mel_list(label_inp,op_info)
      if (idxinp.lt.0) then
        call quit(1,'reo_mel',
     &       'label not on list: '//trim(label_inp))
      end if
      meinp => op_info%mel_arr(idxinp)%mel
      ffinp => meinp%fhand
      opinp => meinp%op
      if (.not.associated(ffinp))
     &     call quit(1,'reo_mel','no file handle defined for '//
     &     trim(meinp%label))
      open_close_inp = ffinp%unit.le.0

      if (open_close_out)
     &     call file_open(ffout)
      if (open_close_inp)
     &     call file_open(ffinp)

      ifrom = fromto/10
      ito = fromto-ifrom*10

      if (ntest.ge.100) then
        write(luout,*) '============================'
        write(luout,*) ' reo_mel messing around'
        write(luout,*) '============================'
        write(luout,*) ' input list = ',trim(meinp%label)
        write(luout,*) ' ffinp: ',trim(ffinp%name)
        write(luout,*) ' opinp: ',opinp%name(1:len_trim(opinp%name))
        write(luout,*) ' adj. : ',dag
        write(luout,*) ' output list = ',trim(meout%label)
        write(luout,*) ' ffout: ',trim(ffout%name)
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name))
        write(luout,*) ' reorder from ',ifrom,' to ',ito
      end if

      if (dag.and.max(iprlvl,ntest).ge.3) write(luout,*)
     &         'Input list will be overwritten by its adjoint.'

      njinp = opinp%njoined
      njout = opout%njoined

      if (ifrom.le.0.or.ito.le.0.or.ifrom.gt.njinp.or.ito.gt.njinp)
     &      call quit(1,'reo_mel','from/to index wrong')

      ! loop over occupation classes
      ! transpose input operator if requested
      if (dag) then
        allocate(transposed(opinp%n_occ_cls))
        transposed = .false.
        do i_occ_cls = 1, opinp%n_occ_cls
          if (transposed(i_occ_cls)) cycle
          iblkinp = (i_occ_cls-1)*njinp+1
          iocc_inp => opinp%ihpvca_occ(1:ngastp,1:2,
     &                                iblkinp:iblkinp+njinp-1)
          ! search corresponding block
          do j_occ_cls = i_occ_cls, opinp%n_occ_cls
            iblkout = (j_occ_cls-1)*njinp+1
            iocc_out => opinp%ihpvca_occ(1:ngastp,1:2,
     &                                iblkout:iblkout+njinp-1)
            if (iocc_equal_n(iocc_inp,.true.,iocc_out,.false.,njinp))
     &         exit
            if (j_occ_cls.eq.opinp%n_occ_cls)
     &         call quit(1,'reo_mel','no adjungate block found!')
          end do

          transposed(i_occ_cls) = .true.
          if (i_occ_cls.eq.j_occ_cls) then
            ! diagonal block: just transpose
            call add_opblk_transp(xdum,0,1d0,meinp,meinp,dag,.false.,
     &                          i_occ_cls,i_occ_cls,
     &                          op_info,str_info,orb_info,.true.)
          else
c dbg
cmh         CAUTION: AD HOC FIX!!!
cmh         (matters for icMRCCSD with singles and doubles)
cmh         We need to understand why we need these sign changes!!!
            tra_sign = 1d0
            if (i_occ_cls.eq.11.and.j_occ_cls.eq.12) then
              tra_sign = -1d0
              print *,'changing sign for blocks ',i_occ_cls,j_occ_cls
            end if
c dbgend
            ! two non-diagonal blocks: transpose both blocks
            transposed(j_occ_cls) = .true.
            ! a) save first block
            lenblkinp = meinp%len_op_occ(i_occ_cls)
            ioffinp = meinp%off_op_occ(i_occ_cls)
            ifree = mem_alloc_real(buffer_inp,lenblkinp,'buffer_inp')
            call get_vec(ffinp,buffer_inp,ioffinp+1,ioffinp+lenblkinp)
            ! b) transpose second block --> first block
            call add_opblk_transp(xdum,0,tra_sign,meinp,meinp,dag,
     &                          .false.,
     &                          j_occ_cls,i_occ_cls,
     &                          op_info,str_info,orb_info,.true.)
            ! c) save transposed first block
            ifree = mem_alloc_real(buffer_out,lenblkinp,'buffer_out')
            call get_vec(ffinp,buffer_out,ioffinp+1,ioffinp+lenblkinp)
            ! d) restore original first block
            call put_vec(ffinp,buffer_inp,ioffinp+1,ioffinp+lenblkinp)
            ! e) transpose first block --> second block
            call add_opblk_transp(xdum,0,tra_sign,meinp,meinp,dag,
     &                          .false.,
     &                          i_occ_cls,j_occ_cls,
     &                          op_info,str_info,orb_info,.true.)
            ! f) restore transposed first block
            call put_vec(ffinp,buffer_out,ioffinp+1,ioffinp+lenblkinp)
          end if
        end do
        deallocate(transposed)
      end if
      ! we assume that the corresponding blocks are in the same order
      do i_occ_cls = 1, opinp%n_occ_cls
        iblkinp = (i_occ_cls-1)*njinp+1
        iblkout = (i_occ_cls-1)*njout+1
        iocc_inp => opinp%ihpvca_occ(1:ngastp,1:2,
     &                                iblkinp:iblkinp+njinp-1)
        iocc_out => opout%ihpvca_occ(1:ngastp,1:2,
     &                                iblkout:iblkout+njout-1)


        ! does this block have to be reordered?
        ! reorder, if from vertex is not zero
        reorder = .not.iocc_zero(iocc_inp(1:ngastp,1:2,ifrom))

        ! get the input, extract diagonal and write to output
        lenblkinp = meinp%len_op_occ(i_occ_cls)
        lenblkout = meout%len_op_occ(i_occ_cls)
        if (lenblkinp.ne.lenblkout) call quit(1,'reo_mel',
     &           'input and output blocks should have same length')
        ioffinp = meinp%off_op_occ(i_occ_cls)
        ioffout = meout%off_op_occ(i_occ_cls)
        if (reorder) then
          if (ntest.ge.100) write(luout,*) 'reorder block no',i_occ_cls
          call reo_mel_blk(meinp,meout,i_occ_cls,i_occ_cls,
     &                     str_info,strmap_info,orb_info,ifrom,ito,
     &                      idxinp)
        else
          ifree = mem_alloc_real(buffer_inp,lenblkinp,'buffer_inp')
          call get_vec(ffinp,buffer_inp,ioffinp+1,ioffinp+lenblkinp)
          call put_vec(ffout,buffer_inp,ioffout+1,ioffout+lenblkinp)
        end if

      end do

      call touch_file_rec(ffout)

      if (open_close_inp)
     &     call file_close_keep(meinp%fhand)
      if (open_close_out)
     &     call file_close_keep(meout%fhand)

      ifree = mem_flushmark('reo_mel')

      call atim_csw(cpu,sys,wall)

      call prtim(luout,'time for reordering ME-list ',
     &                cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
