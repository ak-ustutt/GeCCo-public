*----------------------------------------------------------------------*
      subroutine set_genop2(op,name,type,
     &     min_rank,max_rank,ncadiff,
     &     min_xrank,max_xrank,
     &     hpvx_mnmx,hpvxca_mnmx,
     &     irestr,iformal,freeze,
     &     orb_info)
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
*     iformal:   blocks with this number or more external indices
*                are considered to be to be purely formal.
*     freeze(1): use settings on iadgas to enforce frozen C occupations
*     freeze(2): use settings on iadgas to enforce frozen A occupations
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'stdunit.h'
      include 'ifc_baserout.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(operator), intent(inout) ::
     &     op
      character, intent(in) ::
     &     name*(*)
      integer, intent(in) ::
     &     type, 
     &     min_rank, max_rank, min_xrank, max_xrank,
     &     ncadiff, iformal
      type(orbinf), intent(in), target ::
     &     orb_info
      integer, intent(in) ::
     &     hpvx_mnmx(2,ngastp), hpvxca_mnmx(2,ngastp,2),
     &     irestr(2,orb_info%ngas,2,2)
      logical, intent(in) ::
     &     freeze(2)

      logical, parameter ::
     &     inv_hole = .true.

      logical ::
     &     init_c, init_a, ok
      integer ::
     &     ifree, ipass, irank, na, nc, ica, igas, igasl, idiff, imaxr,
     &     iocc, igastp, iprint, nx, idx, xrank, gasst, gasnd
      integer ::
     &     a_distr(ngastp), c_distr(ngastp), 
     &     a_distr_rv(ngastp), c_distr_rv(ngastp),
     &     hpvxprint(orb_info%ngas)
      integer ::
     &     n_occls_x012(0:2), idx_occls_x012(0:2)
      integer, pointer ::
     &     nspin, ngas, iad_gas(:), hpvxgas(:,:)

      integer, external ::
     &     rank_occ

      iprint = max(iprlvl,ntest)

      if (iprint.ge.100) then
        write(luout,*) '==================='
        write(luout,*) ' set_genop at work'
        write(luout,*) '==================='
        write(luout,*) ' min_rank, max_rank: ',min_rank, max_rank
        write(luout,*) ' min_xrank, max_xrank: ',min_xrank, max_xrank
        write(luout,*) ' ncadiff: ',ncadiff
        write(luout,'(x,a,4(i2,x,i2,2x))')
     &       ' hpvx_mnmx: ',hpvx_mnmx(1:2,1:ngastp)
        write(luout,'(x,a,4(i2,x,i2,2x))')
     &       ' hpvxca_mnmx: ',hpvxca_mnmx(1:2,1:ngastp,1)
        write(luout,'(x,a,4(i2,x,i2,2x))')
     &       '              ',hpvxca_mnmx(1:2,1:ngastp,2)
        call wrt_rstr(luout,irestr,orb_info%ngas)
      end if

      if (len_trim(name).gt.len_opname)
     &     call quit(1,'set_genop','name too long: "'//trim(name)//'"')

      if (type.ne.optyp_operator.and.type.ne.optyp_density) then
        if (type.eq.optyp_intermediate) then
          write(luout,*)
     &         'use set_gen_intermediate to define intermediates'
        else
          write(luout,*) 'type: ',type,' ?'
        end if
        call quit(1,'set_genop','illegal type specification')
      end if

      nspin => orb_info%nspin
      ngas => orb_info%ngas
      iad_gas => orb_info%iad_gas
      hpvxgas => orb_info%ihpvgas

      ! basic settings:
      op%name = '        '
      op%name = name

      op%type = type
      op%njoined = 1  ! always for operators and densities

      op%dagger = .false.

      op%formal=.true.

      ! pass 1: count classes
      ! pass 2: set up occupation information
      do ipass = 1, 2
        
        ! second round: allocate and setup offsets
        if (ipass.eq.2) then
          call init_operator(op,orb_info)
          ! counters according to number of external indices (R12)
          ! to sort operators in the way:
          ! 1st: all operators with no X index (conventional)
          ! 2nd: all operators with  1 X index (auxbasis)
          ! 2nd: all operators with >1 X index (formal operators)
          idx_occls_x012(0) = 0
          idx_occls_x012(1) = n_occls_x012(0)
          idx_occls_x012(2) = n_occls_x012(0)+n_occls_x012(1)
        end if

        op%n_occ_cls = 0
        ! counter for: number of external indices (for R12 stuff)
        n_occls_x012(0:2) = 0
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
              ok = ok.and.hpvxca_mnmx(1,igastp,1).le.c_distr(igastp)
     &               .and.hpvxca_mnmx(2,igastp,1).ge.c_distr(igastp)
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
                ok = ok.and.hpvxca_mnmx(1,igastp,2).le.a_distr(igastp)
     &               .and.hpvxca_mnmx(2,igastp,2).ge.a_distr(igastp)
              end do
              if (.not.ok) cycle a_part

              ! check total number of HPVX
              do igastp = 1, ngastp
                ok = ok.and.
     &                 hpvx_mnmx(1,igastp).le.
     &                        a_distr(igastp)+c_distr(igastp)
     &            .and.hpvx_mnmx(2,igastp).ge.
     &                        a_distr(igastp)+c_distr(igastp)
              end do
              if (.not.ok) cycle a_part

              ! check excitation
              xrank = rank_occ('X',(/c_distr,a_distr/),1)
c dbg
c              print *,'C:',c_distr
c              print *,'A:',a_distr
c              print *,'xrank = ',xrank
c dbg
              ok = xrank.ge.min_xrank.and.xrank.le.max_xrank
              if (.not.ok) cycle a_part

              op%n_occ_cls = op%n_occ_cls + 1
              ! how many X indices?
              if (iextr.gt.0) then ! iextr is set in opdim.h
                nx = min(2,c_distr(iextr)+a_distr(iextr))
              else
                nx = 0
              end if
              n_occls_x012(nx) = n_occls_x012(nx)+1 

              if (ipass.eq.2) then
                ! set occupation of current class
                idx_occls_x012(nx) = idx_occls_x012(nx)+1
                idx = idx_occls_x012(nx)
                ! declare whether this block is formal or not.
                op%formal_blk(idx)=
     &               ((c_distr(iextr)+a_distr(iextr)).ge.iformal)
                op%formal=op%formal.and.op%formal_blk(idx)

                op%ihpvca_occ(1:ngastp,1,idx) = c_distr
                op%ihpvca_occ(1:ngastp,2,idx) = a_distr
                op%ica_occ(1,idx) = sum(c_distr(1:ngastp))
                op%ica_occ(2,idx) = sum(a_distr(1:ngastp))

                ! set restrictions
                do ica = 1, 2
                  do igastp = 1, ngastp
                    gasst = orb_info%idx_gas(igastp)
                    gasnd = gasst-1+orb_info%ngas_hpv(igastp)
                    do igas = gasst, gasnd
                      ! set a/c rank as upper bound
                      idiff = -irestr(1,igas,ica,1)+irestr(2,igas,ica,1)
                      imaxr = min(irestr(2,igas,ica,1),
     &                    op%ihpvca_occ(hpvxgas(igas,1),ica,idx))
c dbg
c                      print *,'igas,idiff,maxr = ',igas,idiff,imaxr
c dbg
                      ! not sure whether this will work in all cases:
                      if (igas.lt.gasnd) then
                        op%igasca_restr(1,igas,ica,1,1:nspin,idx) =
     &                       max(0,imaxr - idiff)
                      else
                        op%igasca_restr(1,igas,ica,1,1:nspin,idx) =
     &                       imaxr
                      end if
                      op%igasca_restr(2,igas,ica,1,1:nspin,idx) =
     &                     imaxr                    
c very quick fix:                    
c                    op%igasca_restr(1:2,igas,ica,1,1:nspin,idx) =
c     &                   imaxr                    
                    end do
                  end do
                end do
                ! post-processing for frozen shells
                do ica = 1, 2
                  if (.not.freeze(ica)) cycle
                  do igas = 1, ngas
                    if (iad_gas(igas).ne.2) then
                      if (igas.eq.1) then                        
                        op%igasca_restr(1:2,igas,ica,1,1:nspin,idx) = 0
                      else
                        op%igasca_restr(1,igas,ica,1,1:nspin,idx) =
     &                      op%igasca_restr(2,igas,ica,1,1:nspin,idx) 
                        op%igasca_restr(1,igas-1,ica,1,1:nspin,idx) =
     &                      op%igasca_restr(2,igas-1,ica,1,1:nspin,idx) 
                      end if
                    end if
                  end do
                end do
                ! mask restriction currently unused
                op%igasca_restr(1:2,1:ngas,1:2,2,1:nspin,idx) = 0
              end if 

            end do a_part
          end do c_part

        end do rank

        if (ipass.eq.1) then
          if (iprlvl.ge.10.or.ntest.ge.10)
     &         write(luout,'(x,3a,i4)')
     &         'Number of occupation classes for ',
     &         trim(name),': ',op%n_occ_cls
        else if (iprlvl.ge.10.or.ntest.ge.100) then

          do igas = 1, ngas
            hpvxprint(igas) = igas
          end do
          if (inv_hole) then
            igasl = 0
            do igas = ngas, 1, -1
              if (hpvxgas(igas,1).eq.1) then
                igasl = igasl+1
                hpvxprint(igas) = igasl
              end if
            end do
          end if
c dbg
c          print *,'formal:',op%formal_blk(1:op%n_occ_cls)
c dbg          
          do iocc = 1, op%n_occ_cls
            call wrt_occ_rstr(luout,iocc,
     &           op%ihpvca_occ(1,1,iocc),
     &           op%igasca_restr(1,1,1,1,1,iocc),ngas,nspin)
c            write(luout,'(/x,a,i4)') 'Occupation Nr. ',iocc
c            call wrt_occ(luout,op%ihpvca_occ(1,1,iocc))
c            write(luout,'(/4x,6(2x,i2,x))') hpvxprint(1:ngas)
c            call wrt_rstr(luout,op%igasca_restr(1,1,1,1,iocc),ngas)
          end do
        end if
  
      end do

      return
      end
*----------------------------------------------------------------------*
