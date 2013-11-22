*----------------------------------------------------------------------*
      subroutine set_graph2(ipass,igraph,
     &     str_info,leny,iwscr,lenwscr,
     &     ngas_hpv,norb_hpv,
     &     igamorb,mostnd,idx_gas,ngas,ngam)
*----------------------------------------------------------------------*
*
*     set up weight arrays using the present info on str_info
*
*     andreas, summer 2006
*
*     patched version for treating specific graph #igraph
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     ipass, igraph, ngas, ngam,
     &     ngas_hpv(ngastp), norb_hpv(ngastp), igamorb(*),
     &     mostnd(2,ngam,ngas), idx_gas(ngastp)
      type(strinf), intent(inout) ::
     &     str_info
      integer, intent(inout) ::
     &     leny(4), lenwscr(3), iwscr(*)

      integer ::
     &     iprint, ihpv, nexc, most, 
     &     len_w4sg, leninf, len_wexit, ndis,
     &     idx_dgm, ims, igam, idis, isum, isum_tot,
     &     ndis_a, idum
      integer, allocatable ::
     &     idss(:), idmin(:), idmax(:)

      logical, external ::
     &     allow_sbsp_dis, next_rvlex

      if (ntest.ge.100) then
        write(lulog,*) '====================='
        write(lulog,*) 'here speaks set_graph'
        write(lulog,*) '====================='
        write(lulog,*) ' ipass = ',ipass
        write(lulog,*) ' setting: '
        write(lulog,*) '   ihpv =', str_info%ispc_typ(igraph)
        write(lulog,*) '   occ  =', str_info%ispc_occ(igraph)
        ihpv = str_info%ispc_typ(igraph)
        write(lulog,*) '   restr = ', 
     &              str_info%igas_restr(1:2,1:ngas_hpv(ihpv),1,1,igraph)
        write(lulog,*) '           ', 
     &              str_info%igas_restr(1:2,1:ngas_hpv(ihpv),2,1,igraph)
      end if

      iprint = max(ntest,iprlvl)
      
      if (iprint.ge.5.and.ipass.eq.2) then
        write(lulog,*) 'Number of strings'
        write(lulog,*) '-----------------'
        write(lulog,*) '  Graph space electr.  MS      length   '
      end if

      if (ipass.eq.1) then
        lenwscr(1:3) = 0
      end if


        ihpv = str_info%ispc_typ(igraph)
        nexc = str_info%ispc_occ(igraph)
        most = mostnd(1,1,idx_gas(ihpv))

        ! Generate weights. First pass is used for dimensioning, the
        ! second for setting the arrays.

        if (ipass.eq.1) then
          call weight_gen(ipass,
     &         leny(1),len_w4sg,leninf,len_wexit, ndis,
     &         idum,idum,
     &         str_info%g(igraph)%yssg,
     &         str_info%g(igraph)%wssg,
     &         idum,
     &         norb_hpv(ihpv),igamorb,
     &         mostnd(1,1,idx_gas(ihpv)),
     &         str_info%igas_restr(1,1,1,1,igraph),
c              !   ADAPT FOR OPEN-SHELL ^^^  (ispin set to 1 here)
     &         nexc,ngam,ngas_hpv(ihpv),
     &         iwscr,iwscr(lenwscr(1)+1),iwscr(lenwscr(1)+lenwscr(2)+1))
        else
          call weight_gen(ipass,
     &         leny(1),len_w4sg,leninf,len_wexit, ndis,
     &         str_info%g(igraph)%y4sg,
     &         str_info%g(igraph)%yinf,
     &         str_info%g(igraph)%yssg,
     &         str_info%g(igraph)%wssg,
     &         str_info%g(igraph)%lenstr_dgm,
     &         norb_hpv(ihpv),igamorb,
     &         mostnd(1,1,idx_gas(ihpv)),
     &         str_info%igas_restr(1,1,1,1,igraph),
c              !   ADAPT FOR OPEN-SHELL ^^^
     &         nexc,ngam,ngas_hpv(ihpv),
     &         iwscr,iwscr(lenwscr(1)+1),iwscr(lenwscr(1)+lenwscr(2)+1))
        end if
        if (ipass.eq.1) then
          ! Set array sizes.
          leny(2) = 3*leninf
          leny(3) = ndis*ngam*(nexc+1)
          leny(4) = (nexc+1)
          lenwscr(1) = max(lenwscr(1),len_w4sg)
          lenwscr(2) = max(lenwscr(2),len_wexit)
          lenwscr(3) = max(lenwscr(3),leninf)
          ! set this one as well
          str_info%g(igraph)%ndis = ndis
        else if (ipass.eq.2) then
          ! set offsets and other info
          str_info%g(igraph)%ndis = ndis

          allocate(idss(nexc),idmin(nexc),idmax(nexc))

          idx_dgm = 1
          do ims = 1, nexc+1
            isum_tot = 0
            do igam = 1, ngam
              isum = 0
              ndis_a = 0
              
              idss(1:nexc) = 1
              idmin(1:nexc) = 1
              idmax(1:nexc) = ngas_hpv(ihpv)
C              ! find lowest allowed distribution
C              lowest_dss: do
C                if (allow_sbsp_dis(idss,nexc,ngas_hpv(ihpv),
C     &               str_info%igas_restr(1,1,1,1,igraph)))
Cc                                   !         ^^^
C     &                             exit lowest_dss
C                if (.not.next_rvlex(nexc,idss,idmin,idmax))
C     &               call quit(1,'set_graph','unexpected case (a)')
C              end do lowest_dss
              
              ! loop over distributions, identify and mark the allowed ones
              ! and sum up lengths of allowed distr. to get offset arrays
              do idis = 1, ndis
c dbg
c               print '(x,a,i3,a,5i3)','idis=',idis,' idss=',idss(1:nexc)
c               print *,'allowed? ',
c     &                       allow_sbsp_dis(idss,nexc,ngas_hpv(ihpv),
c     &                       str_info%igas_restr(1,1,1,1,igraph))
c               print *,'masked? ',
c     &                       allow_sbsp_dis(idss,nexc,ngas_hpv(ihpv),
c     &                       str_info%igas_restr(1,1,2,1,igraph))
c dbg
                if (allow_sbsp_dis(idss,nexc,ngas_hpv(ihpv),
     &                       str_info%igas_restr(1,1,2,1,igraph))) then
                              ! ADAPT FOR OPEN SHELL  ^^^
                  ! this is a masked distribution
                  str_info%g(igraph)%idis_m(idis) = 0
                else
                  str_info%g(igraph)%idis_m(idis) = 1
                  str_info%g(igraph)%ioffstr_dgm(idx_dgm) = isum
                  isum = isum + str_info%g(igraph)%lenstr_dgm(idx_dgm)
c                  idx_dgm = idx_dgm+1
                  ndis_a=ndis_a+1
c dbg
c                  print *,' # str: ',
c    &                      str_info%g(igraph)%lenstr_dgm(idx_dgm)
c dbg
                end if

                idx_dgm = idx_dgm+1

                if (.not.next_rvlex(nexc,idss,idmin,idmax)
     &               .and.idis.ne.ndis)
     &               call quit(1,'set_graph','unexpected case (b)')
              end do
              str_info%g(igraph)%ndis_a = ndis_a
              str_info%g(igraph)%lenstr_gm(igam,ims) = isum
              isum_tot = isum_tot + isum
            end do

            if (iprint.gt.5) then
              if (ims.eq.1) then
                write(lulog,'(4x,4(i3,3x),i12)') igraph, ihpv, nexc,
     &               -nexc+(ims-1)*2, isum_tot
              else
                write(lulog,'(10x,3(i3,3x),i12)')        ihpv, nexc,
     &               -nexc+(ims-1)*2, isum_tot
              end if
            end if
            
          end do
          
          deallocate(idss,idmin,idmax)
          
        else

          call quit(1,'set_graph','illegal value for ipass')

        end if
      
      return
      end

