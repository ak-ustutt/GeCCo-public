*----------------------------------------------------------------------*
      subroutine so_expand(idxval,idxspin,idxperm,n_so,max_so,
     &                     idxprqs,nidx,
     &                     set_aa,set_bb,set_ab,set_ba,set_abx,
     &                     set_ccaa,set_ca)
*----------------------------------------------------------------------*
*     expand to spin-orbitals:
*       for each idxprqs set up requested spin cases (if allowed)
*       and permutations of indices (if CCAA or CA symmetry present)
*
*     set_aa: set alpha/alpha
*     set_bb: set beta/beta
*     set_ab: set alpha/beta
*     set_ab: set beta/alpha
*     set_abx: set exchange type contributions for alpha/beta
*     
*     set_ccaa: set Hermitian conjugate 
*                    pr <-> qs
*     set_ca:   set particle-individual Hermitian conjugate
*                    p <-> q  and r <-> s
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 00

      integer, intent(in) ::
     &     nidx, max_so
      integer, intent(out) ::
     &     n_so, 
     &     idxval(max_so), idxspin(max_so), idxperm(max_so)
      integer(2), intent(in) ::
     &     idxprqs(4,nidx)
      logical, intent(in) ::
     &     set_aa, set_bb, set_ab, set_ba, set_abx, set_ccaa, set_ca

      integer ::
     &     iiso, ii, ip, iq, ir, is,
     &     ncase, pattern, icase, ispcase, idxset
      integer ::
     &     idxcase(6)
      logical ::
     &     col_eq, colcol_eq,
     &     row_eq, rowrow_eq,
     &     dia_eq, diadia_eq

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'so_expand')
        write(luout,*) 'nidx = ',nidx
        write(luout,*) 'idxprqs:'
        write(luout,'(5x,4i5)') idxprqs(1:4,1:nidx)
      end if

      ncase = 0
      idxcase(1:6) = 0
      if (set_aa) ncase = ncase+1
      if (set_aa) idxcase(ncase) = 1
      if (set_bb) ncase = ncase+1
      if (set_bb) idxcase(ncase) = 2
      if (set_ab) ncase = ncase+1
      if (set_ab) idxcase(ncase) = 3
      if (set_ba) ncase = ncase+1
      if (set_ba) idxcase(ncase) = 4
      if (set_abx) ncase = ncase+2
      if (set_abx) idxcase(ncase-1:ncase) = (/5,6/)

      if (ntest.ge.100) write(luout,*) 'idxcase: ',idxcase(1:ncase)

      iiso = 0
      do ii = 1, nidx
        ! get present index quadruple
        ip = idxprqs(1,ii)
        ir = idxprqs(2,ii)
        iq = idxprqs(3,ii)
        is = idxprqs(4,ii)

        if (ip.le.0) cycle ! this entry has been flagged to be ignored

        ! set flags for index patterns
        col_eq =    (ip.eq.iq).or. (ir.eq.is)
        colcol_eq = (ip.eq.iq).and.(ir.eq.is)
        row_eq =    (ip.eq.ir).or. (iq.eq.is)
        rowrow_eq = (ip.eq.ir).and.(iq.eq.is)
        dia_eq =    (ip.eq.is).or. (ir.eq.iq)
        diadia_eq = (ip.eq.is).and.(ir.eq.iq)
        
        ! translate to index pattern
        if (.not.row_eq.and..not.col_eq.and..not.dia_eq) then
          pattern = 9 ! prqs 
        else if (rowrow_eq.and.colcol_eq) then
          pattern = 1 ! pppp
        else if (row_eq.and.col_eq) then
          pattern = 2 ! pqqq and sim.
        else if (rowrow_eq) then
          pattern = 6 ! ppqq
        else if (colcol_eq) then
          pattern = 7 ! pqpq
        else if (diadia_eq) then
          pattern = 8 ! pqqp
        else if (row_eq) then
          pattern = 3 ! ppqs and sim.
        else if (col_eq) then
          pattern = 4 ! prps and sim.
        else
          pattern = 5 ! prqp and sim.
        end if

        ! loop over requested spin cases and set output arrays
        do icase = 1, ncase
          ispcase = idxcase(icase)

          ! select permutations to set according to pattern and spin-case
          select case(pattern)
          case(1)
            if (ispcase.le.2) cycle
            idxset = 2
          case(2)
            if (ispcase.le.2) cycle
            idxset = 3
            if (.not.set_ccaa) idxset = 1
          case(3)
            idxset = 6
            if (.not.set_ca)   idxset = 3
            if (.not.set_ccaa) idxset = 1
            if (ispcase.le.2)  idxset = 5
            if (ispcase.le.2.and..not.set_ca) cycle
          case(4)
            idxset = 3
            if (.not.set_ccaa) idxset = 1
          case(5)
            idxset = 6
            if (.not.set_ca)   idxset = 3
            if (ispcase.le.2)  idxset = 3
            if (.not.set_ccaa) idxset = 1
          case(6)
            idxset = 7
            if (.not.set_ca)   idxset = 8
            if (.not.set_ccaa) idxset = 2
            if (ispcase.le.2) idxset = 4
            if (ispcase.le.2.and..not.set_ca) cycle
          case(7)
            idxset = 1
          case(8)
            idxset = 7
            if (.not.set_ca)   idxset = 8 !??
            if (.not.set_ccaa) idxset = 1 !2
            if (ispcase.le.2)  idxset = 1
          case(9)
            idxset = 6
            if (.not.set_ca)   idxset = 3
            if (.not.set_ccaa) idxset = 1
          end select

          ! process the different cases
          select case(idxset)
          case(1)
            iiso = iiso+1
            idxval(iiso)  = ii
            idxspin(iiso) = ispcase
            idxperm(iiso) = 1
          case(2)
            if (icase.gt.1.and.
     &           (ispcase.eq.4.and.idxcase(icase-1).eq.3).or.
     &           (ispcase.eq.6.and.idxcase(icase-1).eq.5))                   
     &           cycle
            iiso = iiso+1
            idxval(iiso)  = ii
            idxspin(iiso) = ispcase
            idxperm(iiso) = 1
          case(3)
            idxval(iiso+1:iiso+2) = ii
            idxspin(iiso+1:iiso+2) = ispcase
            idxperm(iiso+1:iiso+2) = (/1,2/)
            iiso = iiso+2
          case(4)
            iiso = iiso+1
            idxval(iiso)  = ii
            idxspin(iiso) = ispcase
            idxperm(iiso) = 3
          case(5)
            idxval(iiso+1:iiso+2) = ii
            idxspin(iiso+1:iiso+2) = ispcase
            idxperm(iiso+1:iiso+2) = (/3,4/)
            iiso = iiso+2
          case(6)
            idxval(iiso+1:iiso+4) = ii
            idxspin(iiso+1:iiso+4) = ispcase
            idxperm(iiso+1:iiso+4) = (/1,2,3,4/)
            iiso = iiso+4
          case(7)
            if (icase.gt.1.and.
     &           (ispcase.eq.4.and.idxcase(icase-1).eq.3).or.
     &           (ispcase.eq.6.and.idxcase(icase-1).eq.5))                   
     &           cycle
            idxval(iiso+1:iiso+4) = ii
            idxspin(iiso+1:iiso+4) = ispcase
            idxperm(iiso+1:iiso+4) = (/1,2,3,4/)
            iiso = iiso+4
          case(8)
            if (icase.gt.1.and.
     &           (ispcase.eq.4.and.idxcase(icase-1).eq.3).or.
     &           (ispcase.eq.6.and.idxcase(icase-1).eq.5))                   
     &           cycle
            idxval(iiso+1:iiso+2) = ii
            idxspin(iiso+1:iiso+2) = ispcase
            idxperm(iiso+1:iiso+2) = (/1,2/)
            iiso = iiso+2
          end select

        end do

      end do
      n_so = iiso

      if (ntest.ge.100) then
        write(luout,*) 'on exit: '
        write(luout,*) 'n_so: ',n_so
        write(luout,*) 'idxval idxspin idxperm'
        do ii = 1, n_so
          write(luout,'(x,3i10)') idxval(ii),idxspin(ii),idxperm(ii) 
        end do
      end if

      return
      end
