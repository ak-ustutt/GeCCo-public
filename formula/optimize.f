*----------------------------------------------------------------------*
      subroutine optimize(fl_opt,
     &     op_info,str_info,orb_info,label)
*----------------------------------------------------------------------*
*     automatic optimization of the sequence of binary contractions
*     on fl_opt
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'

      integer, parameter ::
     &     ntest = 00
     
      character(len=8), parameter ::
     &     i_am = 'optimize'
      integer, parameter :: maxpass=10

      type(operator_info), intent(inout) ::
     &     op_info
      type(strinf), intent(inout) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(formula_item), intent(inout), target ::
     &     fl_opt
      character*(*), intent(in) ::
     &     label


      logical ::
     &     opt
      integer ::
     &     idxfl, iprint, ipass, idxopres, lti_cnt
      real(8) ::
     &     cpu, sys, wall, cpu0, sys0, wall0
      type(formula_item), pointer ::
     &     fl_pnt


      ! for each formula target:
      ! loop over formula and find
      !   identical binary contractions
      !     a) both factors are identical
      !        (the same intermediate->calc. only once)
      !     b) one identical factor and same result target
      !        (factor out: sum non-identical operators and multiply only once)
      !     c) one identical factor and result is an intermediate
      !        (potential factor out, if intermediates are summed up due to b))

      if (ntest.gt.0) call write_title(lulog,wst_dbg_subr,i_am)

      if (ntest.ge.100) then
        write(lulog,*) 'formula on entry'
        call print_form_list(lulog,fl_opt,op_info)
      end if

      call atim_csw(cpu0,sys0,wall0)

      iprint = max(ntest,iprlvl)

      fl_pnt => fl_opt
      idxfl = 0

      ! counter for long-term intermediates
      lti_cnt = 0

      main_loop: do

        if (ntest.ge.1000) 
     &      call print_form_item(lulog,idxfl,fl_pnt,op_info)

        if (fl_pnt%command.eq.command_end_of_formula) exit main_loop

        if (fl_pnt%command.ne.command_set_target_init) then
          write(lulog,*) 'I expected find [INIT] at this place, found'
          call print_form_item(lulog,idxfl,fl_pnt,op_info)
          call quit(1,i_am,'[INIT] (or [END]) expected')
        end if

        idxopres = fl_pnt%target        ! op index of result

        if (iprint.ge.3)
     &     write(lulog,*) 'Optimizing for target operator: ',
     &                    trim(op_info%op_arr(idxopres)%op%name)
        
        if (.not.associated(fl_pnt%next))
     &       call quit(1,i_am,'unexpected end of formula list')
        fl_pnt => fl_pnt%next
 
        ipass = 1 !0
        opti_loop: do
          
          if (ntest.ge.100) then
            write(lulog,*) 'calling kernel for ipass = ',ipass
          end if 

          call optimize_kernel(ipass,lti_cnt,opt,fl_pnt,
     &                         op_info,str_info,orb_info)        

          if (ntest.ge.100.and..not.opt) 
     &           write(lulog,*) 'no further optimizations found ...' 

          if (ntest.ge.1000) then
             write(lulog,*) 'formula after pass = ',ipass
             call print_form_list(lulog,fl_opt,op_info)
          end if
          if ((ipass.gt.0.and..not.opt).or.ipass.ge.maxpass) 
     &       exit opti_loop
          ipass = ipass+1
        end do opti_loop

        ! forward to next [INIT] or [END]
        do while(
     &     fl_pnt%command.ne.command_set_target_init.and.
     &     fl_pnt%command.ne.command_end_of_formula)
           if (.not.associated(fl_pnt%next))
     &           call quit(1,i_am,'formula list is broken ... ups!')
           fl_pnt => fl_pnt%next
        end do

      end do main_loop

      if (ntest.ge.150) then
        write(lulog,*) 'formula on exit'
        call print_form_list(lulog,fl_opt,op_info)
      end if
c dbg
       call check_formula_list(fl_opt,op_info)
c dbg

      call atim_csw(cpu,sys,wall)
      call prtim(lulog,'time for optimization',
     &           cpu-cpu0,sys-sys0,wall-wall0)

      return
      end

