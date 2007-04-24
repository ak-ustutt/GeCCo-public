*----------------------------------------------------------------------*
      subroutine set_genop(op,name,dagger,absym,casym,gamma,s2,ms,
     &     min_rank,max_rank,ncadiff,hpvx_mnmx,irestr,
     &     iad_gas,hpvxgas,ngas)
*----------------------------------------------------------------------*
*     set up occupations for a general operator described by
*     min_rank, 
*     max_rank:  minimum and maximum rank of operator defined as 
*                max(C ops,A ops)
*     ncadiff:   number of C minus number of A (= net. number of
*                created particles)
*     hpvx_mnmx: minimum and maximum number of operators per H/P/V/X
*                space per C/A
*     irestr:    restriction on subspaces 
*                min, max. number of operators after completion of
*                subspace within H/P/V/X, for C/A
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 100

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      logical, intent(in) ::
     &     dagger
      integer, intent(in) ::
     &     ngas
      integer, intent(in) ::
     &     absym, casym, gamma, s2, ms,
     &     min_rank, max_rank, ncadiff
      integer, intent(in) ::
     &     hpvx_mnmx(2,ngastp,2), irestr(2,ngas,2,2),
     &     iad_gas(ngas), hpvxgas(ngas)

      logical, parameter ::
     &     inv_hole = .true.

      logical ::
     &     init_c, init_a, ok
      integer ::
     &     ifree, ipass, irank, na, nc, ica, igas, igasl, idiff, imaxr,
     &     iocc, igastp, iprint
      integer ::
     &     a_distr(ngastp), c_distr(ngastp), 
     &     a_distr_rv(ngastp), c_distr_rv(ngastp),
     &     hpvxprint(ngastp)


      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        write(luout,*) '==================='
        write(luout,*) ' set_genop at work'
        write(luout,*) '==================='
        write(luout,*) ' min_rank, max_rank: ',min_rank, max_rank
        write(luout,*) ' ncadiff: ',ncadiff
        write(luout,'(x,a,4(i2,x,i2,2x))')
     &       ' hpvx_mnmx: ',hpvx_mnmx(1:2,1:ngastp,1)
        write(luout,'(x,a,4(i2,x,i2,2x))')
     &       '            ',hpvx_mnmx(1:2,1:ngastp,2)
        call wrt_rstr(luout,irestr,ngas)
      end if

      if (len_trim(name).gt.len_opname)
     &     call quit(1,'set_genop','name too long: "'//trim(name)//'"')

      op%name = '        '
      op%name = name

      op%dagger = dagger
      op%casym = casym
      op%absym = absym

      if (absym.ne.0) stop 'adapt for absym.ne.0'
      if (casym.ne.0) stop 'adapt for casym.ne.0'

      ! set info, some consistency checks would be appropriate as
      ! soon as we seriously use that info
      op%gamt = gamma
      op%s2 =   s2
      op%mst =  ms 

      do ipass = 1, 2

        if (ipass.eq.2) then
          allocate(op%ihpvca_occ(ngastp,2,op%n_occ_cls),
     &             op%ica_occ(2,op%n_occ_cls),
     &             op%igasca_restr(2,ngas,2,2,op%n_occ_cls))
          ifree = mem_register((ngastp*2+2+8*ngas)*op%n_occ_cls,
     &         trim(name)//' occ')
        end if

        op%n_occ_cls = 0
        rank: do irank = min_rank, max_rank
          if (ncadiff.ge.0) then
            nc = irank
            na = nc-ncadiff
            if (na.lt.0) cycle rank
          else
            na = irank
            nc = na+ncadiff
            if (nc.lt.0) cycle rank
          end if
          init_c = .true.
          c_part: do while(next_part_number(init_c,.false.,c_distr_rv,
     &           nc,ngastp,0,nc))

            init_c = .false.

            ! invert sequence delivered by next_part_number:
            do igastp = 1, ngastp
              c_distr(igastp) = c_distr_rv(ngastp+1-igastp)
            end do

            ! check whether distribution is allowed
            ok = .true.
            do igastp = 1, ngastp
              ok = ok.and.hpvx_mnmx(1,igastp,1).le.c_distr(igastp)
     &               .and.hpvx_mnmx(2,igastp,1).ge.c_distr(igastp)
            end do
            if (.not.ok) cycle c_part

            init_a = .true.
            a_part: do while(next_part_number(init_a,.false.,a_distr_rv,
     &           na,ngastp,0,na))
              init_a = .false.

              ! invert sequence delivered by next_part_number:
              do igastp = 1, ngastp
                a_distr(igastp) = a_distr_rv(ngastp+1-igastp)
              end do

              ! check whether distribution is allowed
              ok = .true.
              do igastp = 1, ngastp
                ok = ok.and.hpvx_mnmx(1,igastp,2).le.a_distr(igastp)
     &               .and.hpvx_mnmx(2,igastp,2).ge.a_distr(igastp)
              end do
              if (.not.ok) cycle a_part

              op%n_occ_cls = op%n_occ_cls + 1
              if (ipass.eq.2) then
                ! set occupation of current class
                op%ihpvca_occ(1:ngastp,1,op%n_occ_cls) = c_distr
                op%ihpvca_occ(1:ngastp,2,op%n_occ_cls) = a_distr
                op%ica_occ(1,op%n_occ_cls) = sum(c_distr(1:ngastp))
                op%ica_occ(2,op%n_occ_cls) = sum(a_distr(1:ngastp))

                ! set restrictions
                do ica = 1, 2
                  do igas = 1, ngas
                    ! set a/c rank as upper bound
                    idiff = - irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
                    imaxr = min(irestr(2,igas,ica,1),
     &                       op%ica_occ(ica,op%n_occ_cls))
                    op%igasca_restr(1,igas,ica,1,op%n_occ_cls) =
     &                   max(0,imaxr - idiff)
                    op%igasca_restr(2,igas,ica,1,op%n_occ_cls) =
     &                   imaxr                    
                  end do
                end do
                ! post-processing for frozen shells
                do ica = 1, 2
                  do igas = 1, ngas
                    if (iad_gas(igas).ne.2) then
                      if (igas.eq.1) then                        
                        op%igasca_restr(1:2,igas,ica,1,op%n_occ_cls) = 0
                      else
                        op%igasca_restr(1,igas,ica,1,op%n_occ_cls) =
     &                      op%igasca_restr(2,igas,ica,1,op%n_occ_cls) 
                        op%igasca_restr(1,igas-1,ica,1,op%n_occ_cls) =
     &                      op%igasca_restr(2,igas-1,ica,1,op%n_occ_cls) 
                      end if
                    end if
                  end do
                end do
                ! mask restriction currently unused
                op%igasca_restr(1:2,1:ngas,1:2,2,op%n_occ_cls) = 0
              end if 

            end do a_part
          end do c_part

        end do rank

        if (ipass.eq.1) then
          if (iprlvl.ge.2)
     &         write(luout,'(x,3a,i4)')
     &         'Number of occupation classes for ',
     &         trim(name),': ',op%n_occ_cls
        else if (iprlvl.ge.5) then

          do igas = 1, ngas
            hpvxprint(igas) = igas
          end do
          if (inv_hole) then
            igasl = 0
            do igas = ngas, 1, -1
              if (hpvxgas(igas).eq.1) then
                igasl = igasl+1
                hpvxprint(igas) = igasl
              end if
            end do
          end if
          do iocc = 1, op%n_occ_cls
            write(luout,'(/x,a,i4)') 'Occupation Nr. ',iocc
            call wrt_occ(luout,op%ihpvca_occ(1,1,iocc))
            write(luout,'(/4x,6(2x,i2,x))') hpvxprint(1:ngas)
            call wrt_rstr(luout,op%igasca_restr(1,1,1,1,iocc),ngas)
          end do
        end if
  
      end do

      return
      end
*----------------------------------------------------------------------*
