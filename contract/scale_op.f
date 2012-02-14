      subroutine scale_op(label_res,mode,idx_blk,fac,label_inp,nblk,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     scale blocks of operator list
*     mode == 1:
*       by factor fac,
*       if nblk==-1, all blocks are scaled with the same factor
*     mode == 2:
*       label_inp(2) contains scalar ME-list with respective scaling
*       factor
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

      integer, intent(in) ::
     &     mode, nblk, idx_blk(*)
      real(8), intent(in) ::
     &     fac(*)
      character(*), intent(in) ::
     &     label_res, label_inp(2)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info

      type(me_list), pointer ::
     &     me_res, me_inp, me_fac

      integer ::
     &     idx_res, idx_inp, idx_fac, idx, idxnd,
     &     njoined, iblk, ipri
      logical ::
     &     open_close_res, open_close_inp, open_close_fac,
     &     same
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0, xnorm2, factor
      integer, pointer  ::
     &     occ(:,:,:)

      integer, external ::
     &     idx_mel_list, vtx_type

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(luout,*) '===================='
        write(luout,*) ' scale operators   '
        write(luout,*) '===================='
        write(luout,*) 'Result: ',trim(label_res)
        if (mode.eq.1) then
          write(luout,*) 'The factors and summands: '
          do idx = 1, nblk
            write(luout,'(3x,f12.6,x,i8)') fac(idx),idx_blk(idx)
          end do
          if (nblk.lt.0) then
            write(luout,'(3x,f12.6,x,a)')
     &           fac(idx),'applied to all blocks'
          end if
        else
          write(luout,*) 'Scaling factor(s) on: ',trim(label_inp(2))
        end if
      endif

      idx_res = idx_mel_list(label_res,op_info)
      idx_inp = idx_mel_list(label_inp(1),op_info)

      if (mode.ge.2) then
        idx_fac = idx_mel_list(label_inp(2),op_info)
      else
        idx_fac = 1
      end if

      if (idx_res.lt.0) then
        write(luout,*) '"',trim(label_res),'"'
        write(luout,*) idx_res
        call quit(1,'scale_op','label not on list (1)')
      end if
      if (idx_inp.lt.0) then
        write(luout,*) '"',trim(label_inp(1)),'"'
        write(luout,*) idx_inp
        call quit(1,'scale_op','label not on list (2)')
      end if
      if (idx_fac.lt.0) then
        write(luout,*) '"',trim(label_inp(2)),'"'
        write(luout,*) idx_fac
        call quit(1,'scale_op','label not on list (3)')
      end if

      same = idx_res.eq.idx_inp

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      if (.not.associated(me_res%fhand))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = me_res%fhand%unit.le.0

      if(open_close_res)then
        call file_open(me_res%fhand)
      endif

      me_inp => op_info%mel_arr(idx_inp)%mel
      if (.not.same) then
        if (.not.associated(me_inp%fhand))
     &     call quit(1,'add_op','no file handle defined for '//
     &                  trim(me_inp%label))
        open_close_inp = me_inp%fhand%unit.le.0
        if (open_close_inp) then
          call file_open(me_inp%fhand)
        endif
      else
        open_close_inp = .false.
      end if

      if (mode.eq.1) then
        idxnd = nblk
        if (nblk.lt.0) idxnd = me_inp%op%n_occ_cls
        open_close_fac = .false.
      else
        idxnd = me_inp%op%n_occ_cls
        
        me_fac => op_info%mel_arr(idx_fac)%mel

        if (vtx_type(me_fac%op).ne.vtxtyp_scalar)
     &       call quit(1,'scale_op',
     &       trim(me_fac%label)//' is not a scalar')

        open_close_fac = me_fac%fhand%unit.le.0
        if (open_close_fac) then
          call file_open(me_fac%fhand)
        endif
        
      end if

      njoined = me_res%op%njoined

      ! needed: outer loop over active records

      ! load factor
      if (mode.ge.2) then
        if (me_fac%fhand%buffered) then
          factor = me_fac%fhand%buffer(1)
        else
          call get_vec(me_fac%fhand,factor,1,1)
        end if
        if (mode.ge.3) factor = 1d0/factor
        if (ntest.ge.10) write(luout,*)
     &       'factor from list: ',factor
      else
        factor = fac(1)
      end if

      do idx = 1, idxnd
        if (mode.eq.1.and.nblk.gt.0) then
          iblk = idx_blk(idx)
          factor = fac(idx)
        else
          iblk = idx
        end if

        call scale_opblk(xnorm2,factor,me_inp,me_res,
     &       iblk,iblk,orb_info)

      end do

      ! needed: close loop over active records

      if (open_close_fac)
     &     call file_close_keep(me_fac%fhand)

      if (open_close_res)
     &     call file_close_keep(me_res%fhand)

      if (open_close_inp)
     &     call file_close_keep(me_inp%fhand)

      if (ntest.ge.10) then
        write(luout,*) 'dump of scaled list:'
        if (ntest.ge.10) ipri = 1
        if (ntest.ge.50) ipri = 2
        if (ntest.ge.100) ipri = 3
        if (ntest.ge.500) ipri = 4
        if (ntest.ge.1000) ipri = 5
        call wrt_mel_file(luout,ipri,me_res,
     &       1,me_res%op%n_occ_cls,
     &       str_info,orb_info)
      end if

      call atim_csw(cpu,sys,wall)
      call prtim(luout,'time for scaling',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
