*----------------------------------------------------------------------*
      subroutine set_cc_lagrangian(form_cclag,
     &     nops,ops,idxham,idxlag,idxtop,idxecc)
*----------------------------------------------------------------------*
*
*     set up sequence of operators, integrals and contractions that
*     defines a CC-Lagrangian within the chosen operator space 
*
*     note that definitions of intermediates and optimal sequences
*     is done after the actual target (vector-function, Jacobi-transf.,
*     etc.) has been defined (see routines .... )
*
*     written by andreas, sept. 2006
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_filinf.h'
      include 'ifc_operators.h'
      include 'def_formula.h'
      include 'par_formnames_gen.h'

      type(formula), intent(inout) ::
     &     form_cclag

      integer, intent(in) ::
     &     nops,
     &     idxham,idxlag,idxtop,idxecc

      type(operator), intent(in) ::
     &     ops(nops)

      ! local variables
      logical ::
     &     init_pn, init_pc, len

      type(contraction) ::
     &     contr

      integer ::
     &     lucclag,
     &     maxvtx, maxarc, maxexc, maxcnt, ncommin, ncommax,
     &     iloccls, ihoccls, icomm,
     &     idx, idxh, nterms, ncterm(5),
     &     npart, iexc, ihdh, ihdp, nlop, nlhc, iop, idxarc, ivtxoff
      ! occupations:
      integer ::
     &     iocc_l(ngastp,2), iocc_h(ngastp,2),
     &     iocc_hd(ngastp,2), iocc_hx(ngastp,2),
     &     iocc_hxovl(ngastp,2), iocc_ttot(ngastp,2),
     &     iocc_ltc(ngastp,2), iocc_scr(ngastp,2)
      
      ! some small arrays
      integer, parameter ::
     &     maxpart = 4   ! maximum 4-fold commutators
      integer ::
     &     iexc_part(maxpart), ihd_part(2,maxpart),
     &     ihd_part_idx(maxpart), ihd(2)
      character ::
     &     name*(form_maxlen_label*2)

      ! for timings:
      real(8) ::
     &     cpu, wall, sys, cpu0, wall0, sys0

      ! external functions
      logical, external ::
c     &     iocc_equal,
     &     next_part_number, next_part_pair
      integer, external ::
     &     iopen_nus,
     &     ielsqsum, ielsum, iblk_occ, maxxlvl_op, ifndmin, ifndmax,
     &     ieqfac


      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_cc_lagrangian'
        write(luout,*) '==============================='
      end if

      call atim_csw(cpu0,sys0,wall0)

      contr%idx_res = idxecc
      contr%iblk_res = 1

      ! L + H + 4times T gives maximum 6 vertices
      maxvtx = 6
      ! we can have at most 8 contraction arcs between operators
      maxarc = 8
      ! for convenience, we allocate the maximum number here
      allocate(contr%vertex(maxvtx),contr%arc(maxarc))

      ! get maximum excitation level of T-operators
      maxexc = maxxlvl_op(ops(idxtop))
      if (ntest.ge.100) write(luout,*) 'max. exc.level of T: ', maxexc

      ! assign canonical name and comment
      form_cclag%label = label_cclg0
      form_cclag%comment = title_cclg0

      ! init file
      write(name,'(a,".fml")') label_cclg0
      call file_init(form_cclag%fhand,name,ftyp_sq_unf,0)

      ! open file
      call file_open(form_cclag%fhand)
      lucclag = form_cclag%fhand%unit
      rewind lucclag

      ! first record: a name
      len = len_trim(title_cclg0)
      write(lucclag) len,title_cclg0
      ! second record: define target 
      write(lucclag) 0,idxecc

      ! next records -- how to obtain the result, written by wrt_contr

      nterms = 0
      ! loop over occupation classes of projection space 
      ! (= Lagrange multipliers)
      write(luout,'(2x,42("-"))')
      write(luout,'(3x,a)') '  L    number of   n-fold commutators'
      write(luout,'(3x,a)') 'class    terms    0    1    2    3    4'
      write(luout,'(2x,42("-"))')
      l_loop: do iloccls = 0, ops(idxlag)%n_occ_cls

        ncterm(1:5) = 0

        contr%nvtx=0
        contr%narc=0

        ! we start with zero --
        ! here, we dump how to obtain the CC-energy
        if (iloccls.eq.0) then
          iocc_l(1:ngastp,1:2)=0
          nlop = 0
        else
          if (ops(idxlag)%dagger) then
            ! daggered operator: interchange c<->a
            iocc_l(1:ngastp,2) =
     &           ops(idxlag)%ihpvca_occ(1:ngastp,1,iloccls)
            iocc_l(1:ngastp,1) =
     &           ops(idxlag)%ihpvca_occ(1:ngastp,2,iloccls)
          else
            iocc_l(1:ngastp,1:2) =
     &           ops(idxlag)%ihpvca_occ(1:ngastp,1:2,iloccls)
          end if
          contr%vertex(1)%idx_op=idxlag
          contr%vertex(1)%iblk_op=iloccls
          contr%nvtx = 1
          nlop = 1
        end if

        ! current projection is zeroth order space? find out here ....

        ! loop over blocks of Hamiltonian
        h_loop: do ihoccls = 1, ops(idxham)%n_occ_cls
          
          contr%narc = 0
          contr%nvtx = nlop+1
          iocc_h(1:ngastp,1:2) =
     &         ops(idxham)%ihpvca_occ(1:ngastp,1:2,ihoccls)
          contr%vertex(contr%nvtx)%idx_op=idxham
          contr%vertex(contr%nvtx)%iblk_op=ihoccls
          idxh = contr%nvtx ! remember H vertex number

          ! get excitation part of H
          iocc_hx = iocc_xdn(1,iocc_h)
          ! get overlap with L+ ...
          iocc_hxovl = iocc_overlap(iocc_h,.false.,iocc_l,.true.)
          ! and test whether it is identical with excitation part of H
          ! (as we know that only T operators follow, so excitation part
          ! of H must fully contract with L)
          if (.not.iocc_equal(iocc_hx,.false.,iocc_hxovl,.false.))
     &       cycle h_loop

          ! so, excitation part of H defines the contraction with L
          ! (unless, it is the 0-contraction)
          if (ielsqsum(iocc_hx,ngastp*2).gt.0) then
            nlhc = 1
            contr%narc = 1
            contr%arc(1)%link(1)=1
            contr%arc(1)%link(2)=2
            contr%arc(1)%occ_cnt=iocc_dagger(iocc_hx)
          else
            nlhc = 0
          end if

          ! rest is de-excitation part
          iocc_hd = iocc_add(1,iocc_h,.false.,-1,iocc_hx,.false.)

          ! get number of T that can contract with H (=rank of commutator)
          ncommax = ielsum(iocc_hd,ngastp*2)
          if (ncommax.eq.0) then
            ! excitation part of H must be identical with L+
            if (.not.iocc_equal(iocc_hx,.false.,iocc_l,.true.))
     &           cycle h_loop

            ! set up contraction info (only factor is missing)
            contr%fac = 1d0

            if (ntest.ge.100) then
              call prt_contr(luout,contr,ops)
            end if
            call wrt_contr(lucclag,contr)
            nterms = nterms+1
            ncterm(1) = ncterm(1)+1

          else

            ! total T occupation:
            ! [Ttotal] = [L^\dag] - [Hx] + [Hd^\dag]
            iocc_ttot = iocc_add(1,iocc_l,.true.,-1,iocc_hx,.false.)
            iocc_ttot = iocc_add(1,iocc_ttot,.false.,1,iocc_hd,.true.)
            !  check whether [Ttotal] looks the way we expect:
            iexc = iocc_ttot(ipart,1) ! keep excitation level in mind
            if (.not.(iocc_ttot(ihole,2).eq.iexc.and.
     &           iocc_ttot(ihole,1).eq.0.and.
     &           iocc_ttot(ipart,2).eq.0) ) then
              write(luout,*) 'fishy [Ttot]:'
              call wrt_occ(luout,iocc_ttot)
              call quit(1,'set_cc_lagrangian','fishy occupation')
            end if

            ncommin = iocc_ttot(ipart,1)/(maxexc+1) + 1

            if (ntest.ge.100) write(luout,*) 'ncommin, ncommax: ',
     &           ncommin, ncommax

            do icomm = ncommin, ncommax

              if (icomm.eq.1) then
                contr%fac = 1d0
                contr%nvtx = nlop+2
                contr%vertex(contr%nvtx)%idx_op = idxtop
                idx = iblk_occ(iocc_ttot,.false.,ops(idxtop))
                if (idx.le.0) then
                  write(luout,*) 'occupation not found in list:'
                  write(luout,*) iocc_ttot(1:ngastp,1)
                  write(luout,*) iocc_ttot(1:ngastp,2)
                  stop 'occupation not found in list (1)'
                end if
                contr%vertex(contr%nvtx)%iblk_op = idx

                ! [Hd] defines contraction between H and T
                contr%narc = nlhc+1
                contr%arc(contr%narc)%link(1)=idxh
                contr%arc(contr%narc)%link(2)=contr%nvtx
                contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)
     &               =iocc_hd(1:ngastp,1:2)

                ! [T^\dag]-[Hd] defines contraction between L and T
                iocc_ltc = iocc_add(-1,iocc_hd,.false.,
     &                              +1,iocc_ttot,.true.)
                if (ielsum(iocc_ltc,ngastp*2).gt.0) then
                  contr%narc = contr%narc+1
                  contr%arc(contr%narc)%link(1)=1
                  contr%arc(contr%narc)%link(2)=contr%nvtx
                  contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)=
     &                 iocc_ltc(1:ngastp,1:2)
                end if

                if (ntest.ge.100) then
                  call prt_contr(luout,contr,ops)
                end if
                call wrt_contr(lucclag,contr)
                nterms = nterms+1
                ncterm(2) = ncterm(2)+1

              else
                ! double and higher commutator term:

                ! loop over all n-fold partitionings of [Ttotal] 
                init_pn = .true.
                do while (next_part_number(init_pn,.true.,iexc_part,
     &               iexc,icomm,1,maxexc))
                  init_pn = .false.
                  ! max possible contraction length:
                  maxcnt = min(2,ifndmax(iexc_part,1,icomm,1))

                  ! holes to contract:
                  ihd(1) = iocc_hd(ihole,1)
                  ! particles to contract:
                  ihd(2) = iocc_hd(ipart,2)
                  
                  ! loop over all n-fold partitionings of [Hd]
                  !  (incl. non-equivalent) permutations
                  init_pc = .true.
                  pc_loop: do while (
     &                 next_part_pair(init_pc,.false.,ihd_part,
     &                 ihd,icomm,1,maxcnt))
                    init_pc = .false.

                    ! check, whether all contractions are possible
                    do idx = 1, icomm
                      if (ihd_part(1,idx).gt.iexc_part(idx).or.
     &                    ihd_part(2,idx).gt.iexc_part(idx))
     &                     cycle pc_loop
                    end do
                    
                    ! represent each p/h pair for contractions
                    ! by one integer (for ordering purposes)
                    ihd_part_idx(1:icomm) = ihd_part(1,1:icomm)
     &                                 + ihd_part(2,1:icomm)*(maxcnt+1)

                    ! for equivalent T operators, only one of
                    ! the possible permutations of contractions
                    ! needs be generated: allow only the contr.
                    ! with increasing order of ihd_part_idx
                    do idx = 2, icomm
                      if (iexc_part(idx-1).eq.iexc_part(idx).and.
     &                     ihd_part_idx(idx-1).gt.ihd_part_idx(idx))
     &                     cycle pc_loop
                    end do

                    ! set up contraction information
                    ! get factor from equivalent contractions
                    ! cf. my (andreas) notes on origin of factors
                    contr%fac =
     &                   1d0/dble(ieqfac(iexc_part,ihd_part_idx,icomm))
                    
                    contr%nvtx = nlop+1+icomm
                    ivtxoff = nlop+1
                    idxarc = nlhc
                      
                    do iop = 1, icomm
                      iocc_scr(1:ngastp,1:2) = 0
                      iocc_scr(ihole,2) = iexc_part(iop)
                      iocc_scr(ipart,1) = iexc_part(iop)
                      idx = iblk_occ(iocc_scr,.false.,ops(idxtop))
                      if (idx.le.0) then
                        write(luout,*) 'occupation not found in list:'
                        write(luout,*) iocc_scr(1:ngastp,1)
                        write(luout,*) iocc_scr(1:ngastp,2)
                        stop 'occupation not found in list (2)'
                      end if
                      contr%vertex(ivtxoff+iop)%idx_op = idxtop
                      contr%vertex(ivtxoff+iop)%iblk_op = idx
                      idxarc = idxarc+1
                      contr%arc(idxarc)%link(1)=idxh
                      contr%arc(idxarc)%link(2)=ivtxoff+iop
                      contr%arc(idxarc)%occ_cnt(1:ngastp,1:2) = 0
                      contr%arc(idxarc)%occ_cnt(1,1)
     &                     = ihd_part(1,iop)
                      contr%arc(idxarc)%occ_cnt(2,2)
     &                     = ihd_part(2,iop)
                      ! contraction between L and T
                      iocc_ltc = iocc_add(-1,
     &                   contr%arc(idxarc)%occ_cnt,.false.,
     &                   +1,iocc_scr,.true.)
                      if (ielsum(iocc_ltc,ngastp*2).gt.0) then
                        idxarc = idxarc+1
                        contr%arc(idxarc)%link(1)=1
                        contr%arc(idxarc)%link(2)=nlop+1+iop
                        contr%arc(idxarc)%occ_cnt(1:ngastp,1:2)=
     &                       iocc_ltc(1:ngastp,1:2)
                      end if
                    end do
                    contr%narc = idxarc

                    if (ntest.ge.100) then
                      call prt_contr(luout,contr,ops)
                    end if
                    call wrt_contr(lucclag,contr)
                    nterms = nterms+1
                    ncterm(1+icomm) = ncterm(1+icomm)+1

                  end do pc_loop ! [Hd] partitioning

                end do ! [Ttotal] partitioning

              end if ! number of commutators switch

            end do ! loop over icomm

          end if ! no-commutator-at-all exception

        end do h_loop

c        write(luout,'(4x,i4,3x,iluout,x,5(x,i4))')
c dbg -- add "??" mark for grepping
        write(luout,'(2x,"??",i4,3x,i6,x,5(x,i4))')
     &       iloccls,sum(ncterm(1:5)),ncterm(1:5)

      end do l_loop
      write(luout,'(2x,42("-"))')

      call file_close_keep(form_cclag%fhand)
      deallocate(contr%vertex,contr%arc)

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      end
