      subroutine add_1b_op_ham(label_add,fac,nadd,
     &     op_info,orb_info,str_info)
*----------------------------------------------------------------------*
*     add ME of one body operator to hamiltonian 
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'
      include 'def_target.h'

      integer, intent(in) ::
     &     nadd
      real(8), intent(in) ::
     &     fac(nadd)
      character(*), intent(in) ::
     &     label_add(nadd)
      type(operator_info), intent(inout) ::
     &     op_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(strinf), intent(in) ::
     &     str_info

      type(me_list), pointer ::
     &     me_res, me_current

      character(len_command_par) ::
     &     label_res
      integer ::
     &     iadd,idx_res,idx_add(nadd),
     &     njoined, nblk, iblk, iblkoff, iblk_res

      logical ::
     &     open_close_res, open_close_sum
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0,xnorm
      integer, pointer  ::
     &     occ(:,:,:)

      integer, external ::
     &     idx_mel_list, iblk_occ

      call atim_csw(cpu0,sys0,wall0)

      if(ntest.ge.100)then
        write(lulog,*) '===================='
        write(lulog,*) ' Add up operators   '
        write(lulog,*) '===================='
        write(lulog,*) 'Result: ',trim(mel_ham)
        write(lulog,*) 'The factors and summands:', nadd
        do iadd = 1, nadd
          write(lulog,'(3x,f12.6,x,a)') fac(iadd),
     &         trim(label_add(iadd))
        end do
      endif


      idx_res = idx_mel_list(mel_ham,op_info)

      if (idx_res.lt.0) then
        write(lulog,*) '"',trim(label_res),'"'
        write(lulog,*) idx_res
        call quit(1,'add_op','hamiltonian not on list')
      end if

      ! Point to the relevant operators and their associated files.
      me_res => op_info%mel_arr(idx_res)%mel
      if (.not.associated(me_res%fhand))
     &     call quit(1,'add_1b_op_ham','no file handle defined for '//
     &                  trim(me_res%label))
      open_close_res = me_res%fhand%unit.le.0


      do iadd = 1, nadd
        idx_add(iadd) = idx_mel_list(label_add(iadd),op_info)
        if (idx_add(iadd).lt.0) then
          write(lulog,*) '"',trim(label_add(iadd)),'"'
          write(lulog,*) idx_add
          call quit(1,'add_1b_op_ham','label not on list')
        end if
        me_current => op_info%mel_arr(idx_add(iadd))%mel
        if (.not.associated(me_current%fhand))
     &       call quit(1,'inv_op','no file handle defined for '//
     &                  trim(me_current%label)) 
      end do

      if(open_close_res)then
        call file_open(me_res%fhand)
      endif

      do iadd = 1, nadd
        me_current => op_info%mel_arr(idx_add(iadd))%mel 
        nblk    = me_current%op%n_occ_cls
        njoined = me_current%op%njoined
        open_close_sum = me_current%fhand%unit.le.0
        if(open_close_sum)then
          call file_open(me_current%fhand)
        endif
        
        ! loop over result blocks
        do iblk = 1, nblk
          iblkoff = (iblk-1)*njoined
          occ => me_current%op%ihpvca_occ(1:ngastp,1:2,
     &                               iblkoff+1:iblkoff+njoined)

          iblk_res = iblk_occ(occ,.false.,me_res%op,
     &                        me_current%op%blk_version(iblk))

          if (iblk_res.lt.1) then
            call wrt_occ_n(lulog,occ,njoined)
            call quit(1,'add_1b_op_ham',
     &           'block not found: '//trim(me_res%op%name))
          end if

          call add_opblk(xnorm,1,fac(iadd),me_current,me_res,
     &         iblk,iblk_res,orb_info,.false.)

        end do

        if (open_close_sum)
     &     call file_close_keep(me_current%fhand)

      end do

      if (open_close_res)
     &     call file_close_keep(me_res%fhand)

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'time for adding',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
