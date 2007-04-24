*----------------------------------------------------------------------*
      subroutine set_cc_operators(op_list,nops,orb_info)
*----------------------------------------------------------------------*
*     hard-wired set-up of the basic operators needed for the
*     coupled-cluster Lagrangian
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'ifc_input.h'
      include 'stdunit.h'
      include 'def_orbinf.h'
      include 'def_operator.h'
      include 'def_operator_list.h'

      type(operator_list), intent(inout), target ::
     &     op_list
      integer, intent(inout) ::
     &     nops
      type(orbinf) ::
     &     orb_info
      type(operator_list), pointer ::
     &     list_pnt
      logical ::
     &     dagger
      character ::
     &     name*(len_opname)
      integer ::
     &     absym, casym, s2, ms, min_rank, max_rank, ncadiff,
     &     gamma, iarr(1)
      integer ::
     &     ihpv_mnmx(2,ngastp,2), irestr(2,orb_info%ngas,2,2)

      ! advance to end of operator list:
      list_pnt => op_list
      do while (associated(list_pnt%next))
        list_pnt => list_pnt%next
      end do
      ! is last entry already in use?
      if (associated(list_pnt%op)) then
        allocate(list_pnt%next)
        list_pnt%next%prev => list_pnt
        list_pnt => list_pnt%next
        nullify(list_pnt%next)
c        nullify(list_pnt%op)
      end if
      allocate (list_pnt%op)

      nops = nops+1
      ! new entry: the Hamiltonian
      name = 'H'
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      min_rank = 0
      max_rank = 2
      ncadiff = 0

      call set_hpvx_and_restr_for_h()

      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      ! new entry: the T operator
      nops = nops+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
      nullify(list_pnt%op)
      allocate (list_pnt%op)

      name = 'T'
      dagger = .false.
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      call get_argument_value('method.CC','minexc',ival=min_rank)
      call get_argument_value('method.CC','maxexc',ival=max_rank)
      ncadiff = 0
      call set_hpvx_and_restr_for_xop()

      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      ! new entry: the Tbar operator
      nops = nops+1
      allocate(list_pnt%next)
      list_pnt%next%prev => list_pnt
      list_pnt => list_pnt%next
      nullify(list_pnt%next)
c      nullify(list_pnt%op)
      allocate (list_pnt%op)
     
      name = 'TBAR'
      ! we define an excitation operator to ensure same
      ! storage sequence as for T
      dagger = .true.  ! but we consider the conjugate
      absym = 0
      casym = 0
      gamma = 1
      s2 = 0
      ms = 0
      ! min_rank and max_rank are still set
      ncadiff = 0
      call set_hpvx_and_restr_for_xop()

      call set_genop(list_pnt%op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,ihpv_mnmx,irestr,
     &     orb_info%iad_gas,orb_info%ihpvgas,orb_info%ngas)

      return
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      contains
*----------------------------------------------------------------------*
*     some loops extracted for better overview:
*----------------------------------------------------------------------*
      subroutine set_hpvx_and_restr_for_h()
*----------------------------------------------------------------------*

      implicit none

      integer ::
     &     ica, igastp, igas

      orb_info%nactt_hpv(ngastp)=2
      do ica = 1, 2
        do igastp = 1, ngastp
          if (orb_info%nactt_hpv(igastp).gt.0) then
            ihpv_mnmx(1,igastp,ica) = min_rank
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
          irestr(1,igas,ica,1) = min_rank
          irestr(2,igas,ica,1) = max_rank
        end do
      end do

      return
      end subroutine set_hpvx_and_restr_for_h
*----------------------------------------------------------------------*

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
