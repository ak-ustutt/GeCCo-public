*----------------------------------------------------------------------*
      subroutine get_exc_restr(excrestr,maxh,maxp,nactel,nactorb)
*----------------------------------------------------------------------*
*     get the restrictions for the excitation classes from input
*
*     matthias, march 2011
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'ifc_input.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     maxh, maxp, nactel, nactorb
      integer, intent(out) ::
     &     excrestr(0:maxh,0:maxp,1:2)

      integer ::
     &     icnt, ncnt, ip, ih, excrestr0(6),
     &     minexc, maxexc, minh, minp, maxv, maxvv, gno
      logical ::
     &     pure_vv, l_icci
      character ::
     &     idx_str*1, triples*1

      ! get minimum and maximum numbers of excitations, holes, particles,
      ! valence-valence excitations
      call get_argument_value('method.MR','minh',
     &     ival=minh)
      call get_argument_value('method.MR','minp',
     &     ival=minp)
      call get_argument_value('method.MR','maxv',
     &     ival=maxv)
      call get_argument_value('method.MR','maxvv',
     &     ival=maxvv)
      call get_argument_value('method.MR','minexc',
     &     ival=minexc)
      call get_argument_value('method.MR','maxexc',
     &     ival=maxexc)
      ! which normal ordering is used?
      call get_argument_value('method.MR','GNO',
     &     ival=gno)
      ! icMRCI calculation?
      l_icci = is_keyword_set('method.MRCI').gt.0

      ! minimum excitation: maxexc for GNO=0
      if (gno.eq.0.and.l_icci) minexc = maxexc
      if (maxv.lt.0) maxv = 2*maxexc
      if (maxvv.lt.0) maxvv = maxexc

      call get_argument_value('method.MR','pure_vv',
     &     lval=pure_vv)

      call get_argument_value('method.MR','triples',
     &     str=triples(1:1))

      excrestr(0:maxh,0:maxp,1) = minexc
      excrestr(0:maxh,0:maxp,2) = maxexc
      if (minh.gt.0) then
        excrestr(0:minh-1,0:maxp,1) = 1
        excrestr(0:minh-1,0:maxp,2) = 0
      end if
      if (minp.gt.0) then
        excrestr(0:maxh,0:minp-1,1) = 1
        excrestr(0:maxh,0:minp-1,2) = 0
      end if
      do ip = minp, maxp
        do ih = minh, maxh
          excrestr(ih,ip,1) = max(ih,ip)
          excrestr(ih,ip,2) = min(excrestr(ih,ip,2),(ih+ip+maxv)/2)
          excrestr(ih,ip,2) = min(excrestr(ih,ip,2),max(ih,ip)+maxvv)
          ! cannot annihilate more active electrons than defined
          excrestr(ih,ip,2) = min(excrestr(ih,ip,2),nactel+ih)
          ! cannot in total create more act. el. than empty active orbs
          if (ih-ip.gt.2*nactorb-nactel) excrestr(ih,ip,2) = 0
          if (excrestr(ih,ip,2).lt.excrestr(ih,ip,1)) then
            excrestr(ih,ip,1) = 1
            excrestr(ih,ip,2) = 0
          end if
        end do
      end do
      ncnt = is_keyword_set('method.MR')
      do icnt = 1, ncnt
        call get_argument_value('method.MR','excrestr',
     &       keycount=icnt,iarr=excrestr0(1:6))
        if (excrestr0(1).lt.0) excrestr0(1) = minh
        if (excrestr0(2).lt.0) excrestr0(2) = maxh
        if (excrestr0(3).lt.0) excrestr0(3) = minp
        if (excrestr0(4).lt.0) excrestr0(4) = maxp
        if (excrestr0(5).lt.0) excrestr0(5) = minexc
        if (excrestr0(6).lt.0) excrestr0(6) = maxexc
        do ip = excrestr0(3), excrestr0(4)
          do ih = excrestr0(1), excrestr0(2)
            excrestr(ih,ip,1) = max(excrestr(ih,ip,1),excrestr0(5))
            excrestr(ih,ip,2) = min(excrestr(ih,ip,2),excrestr0(6))
            if (excrestr(ih,ip,2).lt.excrestr(ih,ip,1)) then
              excrestr(ih,ip,1) = 1
              excrestr(ih,ip,2) = 0
            end if
          end do
        end do
      end do
      if (.not.pure_vv) then
        ! CI: allow scalar operator
        excrestr(0,0,1:2) = 0
        ! else: forbid scalar operator
        if (.not.l_icci) excrestr(0,0,1) = 1
      end if

      ! triples models: exclude certain triples blocks if requested
      select case(triples)
      case('F','f') ! do nothing
      case('0','A','a','B','b','C','c','D','d','E','e')
        ! exclude a_uvi^wxy, a_uvw^xya (5 active lines)
        if (maxh.ge.1) then
          if (excrestr(1,0,2).ge.3) excrestr(1,0,2) = 2
        end if
        if (maxp.ge.1) then
          if (excrestr(0,1,2).ge.3) excrestr(0,1,2) = 2
        end if
        select case(triples)
        case('0','A','a','B','b','C','c')
          ! exclude a_uvi^wxa (4 active lines)
          if (maxh.ge.1.and.maxp.ge.1) then
            if (excrestr(1,1,2).ge.3) excrestr(1,1,2) = 2
          end if
        end select
        select case(triples)
        case('0','A','a','B','b','D','d')
          ! exclude a_uij^vwx, a_uvw^xab (4 active lines)
          if (maxh.ge.2) then
            if (excrestr(2,0,2).ge.3) excrestr(2,0,2) = 2
          end if
          if (maxp.ge.2) then
            if (excrestr(0,2,2).ge.3) excrestr(0,2,2) = 2
          end if
        end select
        select case(triples)
        case('0','A','a')
          ! exclude a_uij^vwa, a_uvi^xab (3 active lines)
          if (maxh.ge.2.and.maxp.ge.1) then
            if (excrestr(2,1,2).ge.3) excrestr(2,1,2) = 2
          end if
          if (maxp.ge.2.and.maxh.ge.1) then
            if (excrestr(1,2,2).ge.3) excrestr(1,2,2) = 2
          end if
        end select
        if (triples.eq.'0') then
          ! exclude a_uij^vab (2 active lines)
          if (maxh.ge.2.and.maxp.ge.2) then
            if (excrestr(2,2,2).ge.3) excrestr(2,2,2) = 2
          end if
        end if
      case default
        call quit(1,'get_exc_restr','unknown triples model: '//triples)
      end select

      if (ntest.ge.100) then
        write(luout,'(40("-"))')
        write(luout,'(x,a)') 'Excitation class restrictions:'
        write(idx_str(1:1),'(i1)') maxp+1
        write(luout,'(x,a,'//idx_str//'i8)') 'n_h\n_p',(ip,ip=0,maxp)
        do ih = 0, maxh
          write(luout,'(x,i7,'//idx_str//'(i6,"-",i1))')
     &         ih,((excrestr(ih,ip,icnt),icnt=1,2),ip=0,maxp)
        end do
        write(luout,'(40("-"))')
      end if

      return
      end
