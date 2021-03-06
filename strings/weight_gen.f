*----------------------------------------------------------------------*
      subroutine weight_gen(ipass,
     &     len_y4sg_total,len_w4sg_max,len_info,len_wexit, ndis,
     &     iy4sg,iy_info,iyssg,iwssg,lenstr,
     &     norb,igamorb,mostnd,mnmxspc,nelmax,ngam,nspc,
     &     iw4sg_scr,iwexit,iwex_info)
*----------------------------------------------------------------------*
*     driver routine for weight array generation
*
*     ipass = 1: get dimensions (iyssg and iwssg are needed at this
*                stage: length is iyssg(nelmax*nspc) and
*                                 iwssg((nelmax+1)*nspc) )
*                iy_info has 3*len_info
*                iwex_info needs 1*len_info
*     ipass = 2: set arrays
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer ::
     &     ntest = 00

      integer, intent(in) ::
     &     ipass, 
     &     norb, nelmax, ngam, nspc,
     &     igamorb(norb), mostnd(2,ngam,nspc), mnmxspc(2,nspc)

      integer, intent(out) ::
     &     len_y4sg_total, len_w4sg_max, len_info, len_wexit, ndis,
     &     iy4sg(*), iy_info(3,*),
     &     iyssg(nelmax,nspc), iwssg(0:nelmax,nspc), lenstr(*)

      integer, intent(inout) ::
     &     iw4sg_scr(*), iwexit(*), iwex_info(*)

      integer ::
     &     msmax, idxy, ispc, norb_spc, idx_spc, nel_prev, nel_left,ims,
     &     idx_csg, idx_psg, idxw, idxwp, leny, lenw, idis, nel_psg,
     &     ispcp, iprint, nelmax_spc
      integer ::
     &     ientry(-nelmax:nelmax,ngam), idss(nelmax)
      integer, external ::
     &     idxssd, idxssg, lensubspc
      logical, external ::
     &     next_ssd, allow_sbsp_dis
c dbg
      integer idxww, imss
c dbg

      iprint = max(ntest,iprlvl)

      if (ntest.ge.50) then
        write(lulog,*) '========================='
        write(lulog,*) ' Weight array generation'
        write(lulog,*) '========================='
        write(lulog,*)
        write(lulog,*) ' Number of particles/holes: ', nelmax
        write(lulog,*) ' Number of orbital spaces:  ', nspc
        write(lulog,*) ' Number of orbitals:        ', norb
      end if
      if (ntest.gt.50) write(lulog,*) ' pass = ',ipass

      msmax = nelmax

      ! generate subspace graph vertex and arc weights 
      call weight_ssg(iyssg,iwssg,nelmax,nspc,mnmxspc)
      ! number of distributions = vertex weight(nelmax,nspc)
      ndis = iwssg(nelmax,nspc)

      if (ndis.le.0) then
        write(lulog,*) 'found zero distributions (or less?)'
        write(lulog,*) ' number of electrons: ',nelmax
        write(lulog,*) ' check restriction array: '
        write(lulog,'(x,10(2i3,x))') mnmxspc(1:2,1:nspc)
        call quit(1,'weight_gen','ndis.le.0')
      end if

      if (iprint.ge.50) then
        write(lulog,*) ' Number of subspace distributions: ', ndis
      end if

      ! index of first arc weight array
      idxy = 1
      ! index of first vertex-weight-on-exit array
      idxw = 1

      len_y4sg_total=0
      len_w4sg_max=0
      len_info=0
      len_wexit=0

      ! initialize lenstr
      if (ipass.eq.2) then
        lenstr(1:ndis*ngam*(msmax+1)) = 0
      end if
 
      !-----------------------
      ! loop over subspaces
      !-----------------------
      do ispc = 1, nspc
        norb_spc = mostnd(2,ngam,ispc)-mostnd(1,1,ispc)+1
        idx_spc = mostnd(1,1,ispc)
        nelmax_spc = mnmxspc(2,ispc)

        if (ntest.ge.100) then
          write(lulog,*) '===================================='
          write(lulog,*) 'now in subspace: ',ispc
          write(lulog,*) ' #orbitals: ',norb_spc
          write(lulog,*) ' offset:    ',idx_spc
          write(lulog,*) ' max #elec.:',nelmax_spc
          write(lulog,*) '===================================='
        end if
c dbg
c        if (ipass.eq.2) then
c          print *,'>>>>>> watch at top of ispc-loop'
c        write(lulog,*) 'lengths per dis, IRREP, Ms'
c        idxww=1
c        do imss = 1, msmax+1
c          write(lulog,*) 'ms = ',-msmax+2*(imss-1)
c          write(lulog,*) ' row: distrib, col: IRREP'
c          call iwrtma(lenstr(idxww),ndis,ngam,ndis,ngam)
c          idxww = idxww+ndis*ngam
c        end do
c          print *,'<<<<<<< end of watch'
c        end if
c dbg

        ! no orbitals in this subspace?
        if (norb_spc.le.0) cycle

        !----------------------------------------------------
        ! loop over possible substrings in previous graphs
        !----------------------------------------------------
        ! start with minimal number of electrons in previous subspaces
        if (ispc.eq.1) then
          nel_prev = 0
        else
          nel_prev = mnmxspc(1,ispc-1)
          idss(1:nel_prev) = 1
        end if

        ssd_loop: do
          nel_left = nelmax-nel_prev
          if (ntest.ge.100) 
     &        write(lulog,*) 'ssd_loop: nel_left = ',nel_left
          if (nel_left.le.0) exit ssd_loop

          nel_psg = 0

          ! get address of current and previous graph
          if (ispc.eq.1) then
            ! first space: current index is 1, no previous graph
            idx_csg = 1
            idx_psg = -1
          else
            if (ispc.eq.2) then
              ! second space: previous graph is graph 1, for sure!
              idx_psg = 1
            else
              ! else: get the last subspace from the
              ! subspace distribution
              ispcp = idss(nel_prev)
              ! get the subspace distribution leading to that subspace
              ispcp = ispcp-1
              nel_psg = lensubspc(ispcp,idss,nel_prev)
              ! get the index of the graph of that 
              ! subspace distribution ( idss(1:nel_psg) )
              idx_psg = idxssg(nel_psg,idss,iyssg,iwssg,nelmax,nspc)
            end if
            ! get the index of the graph of the 
            ! subspace distribution idss(1:nel_prev)
            if (nel_left.gt.0) idss(nel_prev+1:nelmax) = ispc
            idx_csg = idxssg(nel_prev,idss,iyssg,iwssg,nelmax,nspc)
          end if

          if (ntest.ge.100) then
            write(lulog,*) '-------------------------------------------'
            write(lulog,*) 'current subspace, ssg: ',ispc, idx_csg
            write(lulog,*) 'previous ssg: ',idx_psg, nel_psg
            write(lulog,*) 'info on ssd: nel_prev = ',nel_prev
            write(lulog,*) '    idss: ',idss(1:nel_prev)
            write(lulog,*) '-------------------------------------------'
          end if

          ! set ientry (info from previous subspace graph)
          if (ipass.eq.2.and.idx_psg.eq.-1) then
            ientry(-msmax:msmax,1:ngam) = 0
            ientry(0,1) = 1
          else if (ipass.eq.2) then
            ! index of exit weights of previous graph
            idxwp = iwex_info(idx_psg)
            ! we have nel_prev electrons set up to now, 
            ! nel_psg have been set in before-previous graphs,
            ! so nel_prev-nel_psg electron were set in prev. graph;
            ! get the appropriate line from iwexit
            call set_ientry(ientry,iwexit(idxwp),
     &           msmax,ngam,nel_prev-nel_psg)
          end if

          if (ntest.ge.100.and.ipass.eq.2) then
            write(lulog,*) 'ientry set to:'
            call iwrtma(ientry,2*msmax+1,ngam,2*msmax+1,ngam)
          end if

          ! generate arc weights
          if (ipass.eq.2) then

            call weight_4sg(iy4sg(idxy),iw4sg_scr,ientry,
     &         norb_spc,igamorb(idx_spc),nel_left,msmax,ngam)

            iy_info(1,idx_csg) = idxy
            iy_info(2,idx_csg) = nel_left
            iy_info(3,idx_csg) = msmax
          end if

          leny = 3*(nel_left+1)*(2*msmax+1)*ngam*norb_spc
          idxy = idxy+leny

          lenw = (nel_left+1)*(2*msmax+1)*ngam*norb_spc
          len_y4sg_total=len_y4sg_total+leny
          len_w4sg_max=max(len_w4sg_max,lenw)
          len_info=len_info+1

          if (ntest.ge.100) then
            write(lulog,*) 'current y length: ',leny
            write(lulog,*) 'current w length: ',lenw
          end if

          ! complete subspace distribution; get index
          if (ipass.eq.2) then
            if (nel_left.gt.0) idss(nel_prev+1:nelmax) = ispc
            idis = idxssd(idss,iyssg,nelmax,nspc)
            if (allow_sbsp_dis(idss,nelmax,nspc,mnmxspc))
     &       ! store length per subspace distribution, IRREP, Ms
     &       call set_lenstr(lenstr,idis,iw4sg_scr,ndis,
     &           nel_left,norb_spc,ngam,msmax)
          end if
          ! not last subspace? 
          if (ispc.lt.nspc) then
              ! keep weights for last orbital level (for ientry in next levels)
            if (ipass.eq.2) then
              call get_iexit(iwexit(idxw),iw4sg_scr,
     &             nel_left,msmax,ngam,norb_spc)
              iwex_info(idx_csg) = idxw
            end if
            idxw = idxw + (nel_left+1)*(2*msmax+1)*ngam
            len_wexit=len_wexit+(nel_left+1)*(2*msmax+1)*ngam

          end if

c dbg
c           if (ipass.eq.2) then
c          print *,'>>>>>> watch at bot of ssd_loop'
c        write(lulog,*) 'lengths per dis, IRREP, Ms'
c        idxww=1
c        do imss = 1, msmax+1
c          write(lulog,*) 'ms = ',-msmax+2*(imss-1)
c          write(lulog,*) ' row: distrib, col: IRREP'
c          call iwrtma(lenstr(idxww),ndis,ngam,ndis,ngam)
c          idxww = idxww+ndis*ngam
c        end do
c          print *,'<<<<<<< end of watch'
c          end if
c dbg

          ! get next subspace distribution
          if (.not.next_ssd(idss,nel_prev,nelmax,ispc-1,mnmxspc))
     &         exit ssd_loop
        end do ssd_loop

      end do

      if (iprint.ge.5.and.ipass.eq.1) then
        write(lulog,'(x,a,i10,a,f6.2,a)')
     &       'Memory for current weight array: ',len_y4sg_total,
     &       ' words (=',dble(len_y4sg_total)*8d0/(1024d0**2),' Mb)'
      end if
      if (ntest.ge.100) then
        write(lulog,*) 'len_y4sg_total = ',len_y4sg_total
        write(lulog,*) 'len_w4sg_max = ',len_w4sg_max
        write(lulog,*) 'len_info = ',len_info
        write(lulog,*) 'len_wexit = ',len_wexit
        write(lulog,*) 'ndis = ',ndis
        write(lulog,*) 'norb = ',norb
      end if

      ! exception: empty space
      if (ipass.eq.2.and.norb.eq.0) then
        lenstr(1:ndis*ngam*(msmax+1)) = 0
      end if

      if (ipass.eq.2.and.(iprint.ge.20.or.ntest.ge.100)) then
        write(lulog,*) 'lengths per dis, IRREP, Ms'
        idxw=1
        do ims = 1, msmax+1
          write(lulog,*) 'ms = ',-msmax+2*(ims-1)
          write(lulog,*) ' row: distrib, col: IRREP'
          call iwrtma(lenstr(idxw),ndis,ngam,ndis,ngam)
          idxw = idxw+ndis*ngam
        end do
      end if

      if (ntest.ge.100) then
        write(lulog,*) '==================='
        write(lulog,*) ' end of weight_gen'
        write(lulog,*) '==================='        
      end if

      end
