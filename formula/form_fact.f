*----------------------------------------------------------------------*
      subroutine form_fact(contr,
     &     op_info,str_info,orb_info)
*----------------------------------------------------------------------*
*     find optimum factorization of a given contraction
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      
      integer, parameter ::
     &     ntest = 0

      type(contraction), intent(inout) ::
     &     contr

      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info

      integer ::
     &     ngas, nsym,
     &     idx, idxmin, narc, nvtx, nfact
      integer, pointer ::
     &     ihpvgas(:)
      integer, allocatable ::
     &     iconn(:,:), ifact(:,:), iperm(:), imult(:)
*----------------------------------------------------------------------*
*     iconn(3,*): describes contraction graph by giving:
*            ( arc-number,  vertex-number 1, vertex-number 2 )
*          only plain (i.e. not fused) arc and vertex numbers occur
*     ifact(4,*): describes factorization by giving
*            ( vertex-number 1, vertex-number 2, 
*                         result-vertex-number, arc-number )
*          all numbers might refer to fused vertices/arcs
*          the number is then an ordered list of all vertices/arcs fused
*          packed using the number of plain arcs (narc) or vertices (nvtx)
*          in the current contraction (as given in contr)
*----------------------------------------------------------------------*
      real(8) ::
     &     cost(3), costmin(3)
      integer ::
     &     iscale(ngastp,3), iscalemin(ngastp,3)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      logical, external ::
     &     next_fact
      integer, external ::
     &     int_pack
     
      if (ntest.ge.100) then
        write(luout,*) '===================='
        write(luout,*) ' form_fact at work!'
        write(luout,*) '===================='
      end if

      call atim(cpu0,sys0,wall0)

      ngas = orb_info%ngas
      nsym = orb_info%nsym
      ihpvgas => orb_info%ihpvgas

      narc = contr%narc
      nvtx = contr%nvtx

      ! if no or only 1 arc is present, we need not bother too much
      if (narc.le.1) then
        ! save factorization info
        if (contr%mxfac.gt.0) deallocate(contr%inffac)
        contr%nfac = narc
        contr%mxfac = narc
        if (narc.eq.1) then
          allocate(contr%inffac(4,narc))
          contr%inffac(1,1) = 1
          contr%inffac(2,1) = 2
          ! naming convention, see function int_pack
          contr%inffac(3,1) = int_pack(contr%inffac,2,3)
          contr%inffac(4,1) = 1
        end if
        return
      end if

      allocate(iconn(3,narc),ifact(4,narc),iperm(narc),imult(narc))

      ! store in more convenient connectivity array
      do idx = 1, narc
        iconn(1,idx) = idx
        iconn(2:3,idx) = contr%arc(idx)%link(1:2)
      end do
      
      ! set up first permutation
      do idx = 1, narc
        iperm(idx) = idx
        imult(idx) = 1
      end do

      idx = 1
      costmin = huge(costmin)
      do 
c        write(luout,'(x,i3,a,20i4)') idx,' next iperm = ',iperm(1:narc)
        ! get next possible factorization
        call mk_fact(ifact,nfact,iperm,iconn,narc,nvtx)

        ! get computational cost of factorization
c        call fact_cost_old(cost,iscale,ifact,nfact,
        call fact_cost(cost,iscale,ifact,nfact,
     &       contr,op_info,str_info,ihpvgas,ngas,nsym)

        ! evaluate cost
        ! preliminary: only cost(1)==time counts
        if (cost(1).lt.costmin(1)) then
          idxmin = idx
          costmin = cost
          iscalemin = iscale
          ! save factorization info
          if (contr%mxfac.gt.0) deallocate(contr%inffac)
          contr%nfac = nfact
          contr%mxfac = nfact
          allocate(contr%inffac(4,nfact))
          contr%inffac(1:4,1:nfact) = ifact(1:4,1:nfact)

        end if

        if (.not.next_fact(iperm,imult,narc,narc,iconn)) exit
        idx = idx+1
      end do

      if (ntest.ge.100) then
        write(luout,*) 'optimal factorization: permutation ',idxmin
        do idx = 1, contr%nfac
          write(luout,*) contr%inffac(1:4,idx)
        end do
        write(luout,'(x,a,g15.10)') 'flops: ',costmin(1)
        write(luout,'(x,a,2g15.10)') 'mem:   ',costmin(2:3)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2)')
     &       'contraction scaling:  ',iscalemin(1:3,1)
        write(luout,'(x,a,"H^",i2," P^",i2,"V^",i2)')
     &       'intermediate scaling: ',iscalemin(1:3,2)
      end if

      deallocate(iconn,ifact,iperm,imult)

      call atim(cpu,sys,wall)
      if (iprlvl.ge.5)
     &    call prtim(luout,'time in form_fact',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return
      end
