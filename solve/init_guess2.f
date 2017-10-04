*----------------------------------------------------------------------*
      subroutine init_guess2(nopt,init,nroots,
     &                me_opt,me_trv,me_dia,me_special,
     &                nspecial,nextra,idxspc,
     &                fl_mvp,depend,fl_spc,nspcfrm,choice,
     &                opti_info,orb_info,op_info,str_info,strmap_info)
*----------------------------------------------------------------------*
*
* get suitable initial guess vectors for the EVP solver
* 
* matthias, nov 2012
*
* improved version, when multi-component vectors are present (iopt>1)
*
* andreas, may 2014
*
* new argument 'choice' added by Pradipta, march 2015
* 'choice' is used to prefer a particular operator over other for the guess
* choice = 0 (default), consider both the operator
* choice = iopt, consider operator iopt, iopt=1,nopt
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'ioparam.h'
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      include 'def_optimize_info.h'
      include 'def_contraction.h'
      include 'def_formula_item.h'
      include 'def_file_array.h'
      include 'mdef_formula_info.h'
      include 'def_dependency_info.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 000
      character(len=12), parameter ::
     &     i_am = 'init_guess2 '

      integer, intent(in) ::
     &     nopt, nroots, nspecial, nextra, idxspc, nspcfrm, choice
      logical, intent(in) ::
     &     init(nopt)
      type(me_list_array), intent(inout) :: 
     &     me_opt(*), me_trv(*), me_dia(*), me_special(nspecial+nextra)
      type(formula_item), intent(inout) ::
     &     fl_mvp, fl_spc(nspcfrm)
      type(dependency_info), intent(in) ::
     &     depend
      type(optimize_info), intent(in) ::
     &     opti_info
      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info

      integer ::
     &     ifree, iopt, jopt, iroot, maxtrials, idx, nset,
     &     nselect, 
     &     jroot, iguess, ntrials_all
      logical ::
     &     trafo, read_in
      real(8) ::
     &     xnrm, xnrm2(nroots), xover, 
     &     fac1, fac2
      type(me_list_array):: 
     &     me_pnt(nopt)

      integer ::
     &     idxset(2)
      real(8) ::
     &     valset(2), trf_nrm

      integer, pointer ::
     &     isign(:), ntrials(:),
     &     idxlist(:,:), idxlist_all(:,:), idxlist_ba(:), idxscr(:)
      real(8), pointer ::
     &     xbuf1(:), xbuf2(:), xlist_all(:), xlist(:,:)
      logical::
     &     list_set
      real(8), external ::
     &     da_ddot

      if (ntest.gt.100) call write_title(luout,wst_dbg_subr,i_am)

      maxtrials = max(1000,4*nroots)
      
      allocate( isign(nopt), ntrials(nopt))
      allocate(idxlist(maxtrials,nopt),idxlist_all(2,maxtrials),
     &     idxlist_ba(maxtrials))
      allocate(xlist(maxtrials,nopt),xlist_all(maxtrials))
      
      do iopt = 1, nopt
         call find_nmin_list(xlist(1:maxtrials,iopt),
     &        idxlist(1:maxtrials,iopt),
     &        maxtrials,me_dia(iopt)%mel)
         ntrials(iopt) = min(maxtrials,me_dia(iopt)%mel%len_op)
      end do
      
      call merge_min_lists(xlist_all,idxlist_all,ntrials_all,
     &     xlist,idxlist,
     &     nopt,maxtrials,ntrials,choice)

      if (has_to_read_from_file_h(init,nopt)) then
         call read_from_file_h(init, nopt, nroots, me_trv, me_opt)
      else

        do iopt =1,nopt
          isign(iopt) = me_trv(iopt)%mel%absym
        end do
        call create_inverted_spin_indexlist_h(
     &       idxlist_ba,idxlist,idxlist_all,ntrials_all,
     &       maxtrials,ntrials,me_dia,isign,nopt,
     &       op_info,str_info,strmap_info,orb_info)


        iroot = 1
        do iguess = 1,ntrials_all
          iopt = idxlist_all(1,iguess)

! transformed preconditioner => transformed initial guess vector
          do jopt = 1,nopt
            me_pnt(jopt)%mel => me_trv(jopt)%mel
          end do
          
          if (opti_info%typ_prc(iopt).eq.optinf_prc_traf) then
            me_pnt(iopt)%mel => me_special(1)%mel
            trafo = .true.
          else if (opti_info%typ_prc(iopt).eq.optinf_prc_traf_spc)then
            me_pnt(iopt)%mel => me_special(4)%mel
            trafo = .true.
          else
            me_pnt(iopt)%mel => me_trv(iopt)%mel
            trafo = .false.
          end if
          
          ! switching all lists to the relevant record
          do jopt = 1,nopt
            call switch_mel_record(me_trv(jopt)%mel,iroot)
            call switch_mel_record(me_pnt(jopt)%mel,iroot)
          end do
          
          if (.not. generate_guess_list_h(
     &         me_pnt,isign,iopt,nopt,
     &         iroot, 
     &         idxlist_all, idxlist_ba, iguess, ntrials_all) ) cycle

          ! if requested, back-transformation of initial guess vector
          if (trafo) then
             call transform_back_wrap(fl_mvp,depend,
     &            me_special,me_pnt(iopt)%mel,me_trv(iopt)%mel,
     &            trf_nrm,
     &            iopt,nspecial,
     &            me_trv(iopt)%mel,
     &            op_info, str_info, strmap_info, orb_info, opti_info)
            
! guess vectors of wrong spin symmetry and twinned guess vectors will be discarded
             if (should_discard_vector_h(trf_nrm, nopt,nroots,
     &            iroot,xnrm2,me_trv, opti_info) )then
               me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
               cycle
             else
               xnrm2(iroot) = trf_nrm**2
             end if 
           end if               ! trafo

          ! project out spin contaminations or other components?
          if (opti_info%typ_prc(iopt).eq.optinf_prc_spinp.or.
     &        opti_info%typ_prc(iopt).eq.optinf_prc_prj.or.
     &         opti_info%typ_prc(iopt).eq.optinf_prc_spinrefp) then
             call  spin_proj_h(opti_info%typ_prc(iopt),
     &            opti_info%nwfpar(iopt),me_trv(iopt)%mel,me_special,
     &            nspecial, nextra, idxspc,
     &            fl_spc, xnrm,
     &            opti_info,orb_info, op_info, str_info, strmap_info)
            if (xnrm**2.lt.1d-12) then
               if (iprlvl.ge.5) write(lulog,*)
     &              'Discarding guess vector due to projection.'
               me_trv(iopt)%mel%fhand%last_mod(iroot) = -1
               cycle
            else if (iroot.gt.1) then
               call orthogonalize_roots_h(iroot, opti_info%nwfpar(iopt),
     &              xnrm2,xnrm**2,me_trv(iopt)%mel)
            else
               xnrm2(iroot) = xnrm**2
            end if
          end if

          ! Assumption: If the controlflow reaches here, a new guess vector was created
          if (iroot.eq.nroots) exit
          iroot = iroot +1
        end do !  iguess
        if (iroot.ne.nroots) call quit(1,i_am,
     &        'Could not find enough guess vectors')
      end if

      deallocate(isign,ntrials)
      deallocate(idxlist,idxlist_all,idxlist_ba,xlist,xlist_all)

      return
      contains
      

      
*----------------------------------------------------------------------*
!>    orthogonalizes the newest root (iroot) to all previous roots
*----------------------------------------------------------------------*
      subroutine orthogonalize_roots_h(iroot,len_op,xnrm2,xnrm,me_root)
*----------------------------------------------------------------------*
      implicit none
      integer, intent(in)::
     &     len_op
      integer, intent(inout)::
     &     iroot
      type(me_list)::
     &     me_root
  
      integer::
     &     ifree,jroot,len_buf
      real(8)::
     &     xnrm2(*),xnrm
      
      real(8),pointer::
     &     xbuf1(:),xbuf2(:)
      real(8)::
     &     xover,fac1,fac2

      ifree = mem_setmark('init_guess.check_guess')
      len_buf = len_op
      ifree = mem_alloc_real(xbuf1,len_buf,
     &     'xbuf1')
      ifree = mem_alloc_real(xbuf2,len_buf,
     &     'xbuf2')

      do jroot = 1,iroot-1
         xover = da_ddot(
     &        me_root%fhand,jroot,
     &        me_root%fhand,iroot,
     &        len_op,
     &        xbuf1, xbuf2,
     &        len_buf)
         if(abs(xover**2/xnrm2(jroot)-xnrm).lt.1d-6)then
            if (iprlvl.ge.5) write(lulog,*)
     &           'Discarding redundant guess vector.'
            me_root%fhand%last_mod(iroot) = -1
            iroot = iroot - 1
            exit
         else if (abs(xover).ge.1d-12) then
! remove this component using norm-conserving factors
            if (iprlvl.ge.5) write(lulog,*)
     &           'Removing overlap to root',jroot,'(',xover,')'
            fac1 = 1d0/sqrt(1d0-xover**2/(xnrm2(jroot)*xnrm))
            fac2 = -sign(1d0/sqrt((xnrm2(jroot)/xover)**2
     &           -xnrm2(jroot)/xnrm),xover)
            call da_vecsum(me_root%fhand,iroot,
     &           me_root%fhand,iroot,fac1,
     &           me_root%fhand,jroot,fac2,
     &           len_op,
     &           xbuf1,xbuf2,
     &           len_buf)
c                  xnrm = xnrm - xover**2/xnrm2(jroot) ! new norm**2
         end if
         if (jroot.eq.iroot-1)xnrm2(iroot) = xnrm
      end do
      ifree = mem_flushmark()
      end subroutine 
*----------------------------------------------------------------------*
!!    determines if the trialvector should be read from a file
*----------------------------------------------------------------------*
      pure function has_to_read_from_file_h(init,nopt) 
*----------------------------------------------------------------------*
      implicit none

      logical ::
     &     has_to_read_from_file_h

      integer,intent(in)::
     &     nopt

      logical,intent(in)::
     &     init(nopt)
      integer::
     &     iopt
      
      has_to_read_from_file_h = .false.
      do iopt = 1, nopt
         has_to_read_from_file_h = has_to_read_from_file_h
     &        .or..not.init(iopt)
      end do
      return
      end function
      
*----------------------------------------------------------------------*
!>    reads a trialvector from the file of the optimized operator
*----------------------------------------------------------------------*
      subroutine read_from_file_h(init,nopt, nroots, me_trv, me_opt) 
*----------------------------------------------------------------------*
      implicit none
      
      
      integer,intent(in)::
     &     nopt, nroots
      logical,intent(in)::
     &     init(nopt)
      type(me_list_array), intent(inout) :: 
     &     me_opt(nopt), me_trv(nopt)

      integer::
     &     iopt, iroot
      
      do iopt = 1, nopt
         if (.not.init(iopt)) then
            do iroot = 1, nroots
               call switch_mel_record(me_trv(iopt)%mel,iroot)
               call switch_mel_record(me_opt(iopt)%mel,iroot)
               call list_copy(me_opt(iopt)%mel,me_trv(iopt)%mel,.false.)
            end do
         else
            do iroot = 1, nroots
               call switch_mel_record(me_trv(iopt)%mel,iroot)
               call zeroop(me_trv(iopt)%mel)
            end do
         end if
      end do
      
      end subroutine
      
*----------------------------------------------------------------------*
!>    reads a trialvector from the file of the optimized operator
*----------------------------------------------------------------------*
      function check_last_guess_h( nopt, iopt,nroots,iroot,guess_norm2,
     &     xnrm2, me_trv,opti_info) 
*----------------------------------------------------------------------*
      implicit none
      logical ::
     &     check_last_guess_h
      
      integer,intent(in)::
     &     nopt, iopt,
     &     nroots,iroot
      real(8),intent(in)::
     &     guess_norm2,xnrm2(nroots)

      type(me_list_array)::
     &     me_trv(nopt)

      type(optimize_info), intent(in)::
     &     opti_info

      
      integer::
     &     ifree
      real(8)::
     &     xover
      real(8),pointer::
     &     xbuf1(:), xbuf2(:)
      real(8),external::
     &     da_ddot

      if (iroot .eq.1) then
         check_last_guess_h = .true.
         return
      end if  

      ifree = mem_setmark('init_guess.check_guess')
      ifree = mem_alloc_real(xbuf1,opti_info%nwfpar(iopt),
     &     'xbuf1')
      ifree = mem_alloc_real(xbuf2,opti_info%nwfpar(iopt),
     &     'xbuf2')
      xover = da_ddot(me_trv(iopt)%mel%fhand,iroot-1,
     &     me_trv(iopt)%mel%fhand,iroot,
     &     opti_info%nwfpar(iopt),
     &     xbuf1, xbuf2,
     &     opti_info%nwfpar(iopt))
      ifree = mem_flushmark()
      if (abs(
     &     xover**2/xnrm2(iroot-1)-guess_norm2).lt.1d-6) then
         check_last_guess_h = .false.
      else
         check_last_guess_h = .true.
      end if  
      return
      end function
      
*----------------------------------------------------------------------*
!>    reads a trialvector from the file of the optimized operator
*----------------------------------------------------------------------*
      subroutine spin_proj_h( typ_prc,len_list,me_root,me_special,
     &     nspecial, nextra, idxspc,
     &     fl_spc, xnrm,
     &     opti_info, orb_info, op_info,str_info, strmap_info) 
*----------------------------------------------------------------------*
      implicit none
      integer,intent(in)::
     &     typ_prc,len_list, nspecial, nextra, idxspc
      real(8),intent(inout)::
     &     xnrm

      type(me_list)::
     &     me_root
      type(me_list_array)::
     &     me_special(nspecial+nextra)
      type(formula_item)::
     &     fl_spc(*)
      type(optimize_info)::
     &     opti_info
      type(orbinf)::
     &     orb_info
      type(operator_info)::
     &     op_info
      type(strinf)::
     &     str_info
      type(strmapinf)::
     &     strmap_info
      
      real(8),pointer::
     &     xbuf1(:),xbuf2(:)
      
      integer::
     &     ifree, len_buf
      len_buf = len_list
      ifree = mem_setmark('init_guess.spin_proj')
      ifree = mem_alloc_real(xbuf1,len_buf,'xbuf1')
      ifree = mem_alloc_real(xbuf2,len_buf,'xbuf2')
      select case(typ_prc)
      case(optinf_prc_spinp)
         call spin_project(me_root,me_special(idxspc)%mel,
     &        fl_spc(1),len_buf,
     &        xbuf1,xbuf2,.true.,xnrm,
     &        opti_info,orb_info,
     &        op_info,str_info,strmap_info)
      case(optinf_prc_spinrefp)
         call spin_project(me_root,me_special(idxspc)%mel,
     &        fl_spc(2),len_buf,
     &        xbuf1,xbuf2,.true.,xnrm,
     &        opti_info,orb_info,
     &        op_info,str_info,strmap_info)
         call evaluate2(fl_spc(1),.false.,.false.,
     &        op_info,str_info,strmap_info,orb_info,
     &        xnrm,.true.)
      case(optinf_prc_prj)
         call evaluate2(fl_spc(1),.false.,.false.,
     &        op_info,str_info,strmap_info,orb_info,
     &        xnrm,.true.)
      end select
      ifree = mem_flushmark()
      return
      end subroutine
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine create_inverted_spin_indexlist_h(
     &     idxlist_ba,idxlist,idxlist_all, ntrials_all,
     &     maxtrials,ntrials,
     &     me_dia,isign,nopt,
     &     op_info, str_info, strmap_info, orb_info)
*----------------------------------------------------------------------*
      implicit none
      
      integer,intent(out)::
     &     idxlist_ba(ntrials_all)

      integer,intent(in)::
     &     maxtrials,nopt,ntrials_all,
     &     ntrials(nopt),
     &     isign(nopt),
     &     idxlist(maxtrials,nopt),
     &     idxlist_all(2,ntrials_all)

      type(me_list_array)::
     &     me_dia(nopt)

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      
      integer::
     &     iopt, idx, iguess
      integer,pointer::
     &     idxscr(:)
      
      allocate(idxscr(maxtrials))
      do iopt = 1, nopt
        if (isign(iopt).eq.0) cycle
        
        call symidx_ab(idxscr,
     &       idxlist(1,iopt),ntrials(iopt),me_dia(iopt)%mel,
     &       op_info,str_info,strmap_info,orb_info)
        
        if (ntest.ge.200) then
          write(lulog,*) 'ntrials(iopt),maxtrials: ',
     &         ntrials(iopt),maxtrials
          write(lulog,*) 'idxscr = '
          write(lulog,'(1x,5i4,x,5i4)') idxscr(1:ntrials(iopt))
        end if
        
        idx = 0
        do iguess = 1, ntrials_all
          if (idxlist_all(1,iguess) .ne.iopt) cycle
          idx = idx+1
          idxlist_ba(iguess) = idxscr(idx)
        end do
      end do
      deallocate(idxscr)
      return
      end subroutine
*----------------------------------------------------------------------*

*----------------------------------------------------------------------*
      function generate_guess_list_h(
     &     me_trv,isign,iopt,nopt,
     &     iroot,
     &     idxlist_all,idxlist_ba,iguess,ntrials_all)
*----------------------------------------------------------------------*
      implicit none

      
      
      integer,intent(in)::
     &     isign(nopt),iguess,nopt,iopt,
     &     iroot,
     &     ntrials_all,
     &     idxlist_all(2,ntrials_all),
     &     idxlist_ba(ntrials_all)

      logical::
     &     generate_guess_list_h

      type(me_list_array)::
     &     me_trv(nopt)
      
      integer::
     &     idxset(2),
     &     nset,
     &     jopt

      real(8)::
     &     valset(2)


      nset = 0
      
      if (isign(iopt) .ne. 0) then
        if (idxlist_all(2,iguess).lt.abs(idxlist_ba(iguess)))then
          nset = 2
          idxset(1) = abs(idxlist_all(2,iguess))
          idxset(2) = abs(idxlist_ba(iguess))
          valset(1) = 1d0/sqrt(2d0)
          valset(2) = dble(isign(iopt))*
     &         dble(sign(1,idxlist_ba(iguess)))
     &         /sqrt(2d0)
        else if (idxlist_all(2,iguess).eq.abs(idxlist_ba(iguess))
     &         .and.isign(iopt).eq.+1) then
          if (idxlist_ba(iguess).lt.0)
     &         call quit(1,i_am,'unexpected case')
! set a single element
          nset = 1
          idxset(1) = abs(idxlist_all(2,iguess))
          valset(1) = 1d0
        end if
      else
! set a single element
        nset = 1
        idxset(1) = abs(idxlist_all(2,iguess))
        valset(1) = 1d0
      end if

      generate_guess_list_h=.false.
      if (nset.eq.0) return
      generate_guess_list_h=.true.
      
      if (ntest.ge.100) then
        write(lulog,*) 'iopt:   ',iopt
        write(lulog,*) 'idxset: ',idxset(1:nset)
        write(lulog,*) 'valset: ',valset(1:nset)
        write(lulog,*) '=> switching to iroot = ',iroot
      end if
      
      do jopt = 1, nopt
        if (jopt.eq.iopt) then
          call set_list(me_trv(jopt)%mel,idxset,valset,nset)
        else
          call zeroop(me_trv(jopt)%mel)
        end if
      end do
      
      return
      end function
      
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      function should_discard_vector_h(trf_nrm, nopt, nroots,iroot,
     &     xnrm2,me_trv,opti_info)
*----------------------------------------------------------------------*
      implicit none

      logical::
     &     should_discard_vector_h
      real(8),intent(in)::
     &     trf_nrm,xnrm2(nroots)

      integer,intent(in)::
     &     nopt,nroots,iroot

      type(me_list_array)::
     &     me_trv(nopt)
      type(optimize_info), intent(in)::
     &     opti_info

      if (trf_nrm.lt.1d-12) then
        if (iprlvl.ge.5) write(lulog,*)
     &       'Discarding guess vector with wrong spin symmetry.'
        should_discard_vector_h =.true.
      else
        if(check_last_guess_h(nopt,iopt,nroots,iroot,
     &       (trf_nrm**2),
     &       xnrm2, me_trv,opti_info) ) then
          should_discard_vector_h =.false.
        else
          if (iprlvl.ge.5) write(lulog,*)
     &         'Discarding twin guess vector.'
          should_discard_vector_h =.true.
        end if
        
      end if 
      return
      end function
      end
