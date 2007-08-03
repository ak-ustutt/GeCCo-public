*----------------------------------------------------------------------*
      subroutine set_cc_operators(op_info,orb_info)
*----------------------------------------------------------------------*
*     hard-wired set-up of the basic operators needed for the
*     coupled-cluster Lagrangian
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'cc_routes.h'
      include 'ifc_input.h'
      include 'par_opnames_gen.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'

      type(operator_info), intent(inout), target ::
     &     op_info
      type(orbinf) ::
     &     orb_info
      type(operator), pointer ::
     &     op_pnt, top_pnt, tbar_pnt, ham_pnt, cclg_pnt
      logical ::
     &     dagger
      character ::
     &     name*(len_opname)
      integer ::
     &     absym, casym, s2, ms, min_rank, max_rank, ncadiff,
     &     gamma, iarr(1), isim, idx, iformal
      integer ::
     &     ihpv_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      integer, external ::
     &     idx_oplist2


      ! new entry: the T operator
      call add_operator(op_top,op_info)
      idx = idx_oplist2(op_top,op_info)
      op_pnt => op_info%op_arr(idx)%op
      top_pnt => op_pnt

      name = op_top
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      ncadiff = 0
      iformal=max_rank+1
      call set_hpvx_and_restr_for_xop()

      call set_genop(op_pnt,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)

      ! new entry: the Tbar operator
      call add_operator(op_tbar,op_info)
      idx = idx_oplist2(op_tbar,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,top_pnt,orb_info)
      ! we define an excitation operator to ensure same
      ! storage sequence as for T
      op_pnt%dagger = .true.  ! but we consider the conjugate

      ! new entry: the CC residual OMG
      call add_operator(op_omg,op_info)
      idx = idx_oplist2(op_omg,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,top_pnt,orb_info)
     
      ! new entry: the DIAgonal
      call add_operator(op_dia1,op_info)
      idx = idx_oplist2(op_dia1,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,top_pnt,orb_info)

      ! new entry: the CC-Lagrangian (scalar)
      ! looks like overkill, but makes life easier
      call add_operator(op_cclg,op_info)
      idx = idx_oplist2(op_cclg,op_info)
      op_pnt => op_info%op_arr(idx)%op
      cclg_pnt => op_pnt

      name = op_cclg
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      min_rank=0
      max_rank=0
      ncadiff = 0
      iformal=1
      call set_hpvx_and_restr_for_xop()

      call set_genop(op_pnt,name,optyp_operator,
     &     dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,iformal,
     &     orb_info)

      ! ... and the CC energy
      call add_operator(op_ccen,op_info)
      idx = idx_oplist2(op_ccen,op_info)
      op_pnt => op_info%op_arr(idx)%op
      call clone_operator(op_pnt,cclg_pnt,orb_info)

      if (solve_tbar) then

        ! new entry: the Tbar.A transform ...
        call add_operator(op_tbar_a,op_info)
        idx = idx_oplist2(op_tbar_a,op_info)
        op_pnt => op_info%op_arr(idx)%op
        ! same as tbar
        idx = idx_oplist2(op_tbar,op_info)
        tbar_pnt => op_info%op_arr(idx)%op
        call clone_operator(op_pnt,tbar_pnt,orb_info)

        ! ... and the RHS for the tbar-equations (eta)
        call add_operator(op_eta,op_info)
        idx = idx_oplist2(op_eta,op_info)
        op_pnt => op_info%op_arr(idx)%op
        call clone_operator(op_pnt,tbar_pnt,orb_info)

      end if

      return
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*
*     some loops extracted for better overview:
*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_xop()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas

      do ica = 1, 2
        do igastp = 1, ngastp
          if (orb_info%nactt_hpv(igastp).gt.0.and.
     &        ((ica.eq.1.and.(igastp.eq.ipart.or.igastp.eq.ivale)).or.
     &         (ica.eq.2.and.(igastp.eq.ihole.or.igastp.eq.ivale)) ))
     &           then
            ihpv_mnmx(1,igastp,ica) = 0
            ihpv_mnmx(2,igastp,ica) = max_rank
          else
            ihpv_mnmx(1,igastp,ica) = 0
            ihpv_mnmx(2,igastp,ica) = 0
          end if
        end do
      end do
      irestr(1:2,1:orb_info%ngas,1:2,1:2) = 0
      do ica = 1, 2
        do igas = 1, orb_info%ngas
          irestr(1,igas,ica,1) = 0
          irestr(2,igas,ica,1) = max_rank
        end do
      end do

      return
      end subroutine set_hpvx_and_restr_for_xop
*----------------------------------------------------------------------*

      end
