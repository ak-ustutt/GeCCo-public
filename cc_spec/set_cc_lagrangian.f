*----------------------------------------------------------------------*
      subroutine set_cc_lagrangian(form_cclag,
     &     nops,ops,idxecc,idxham,idxlag,idxtop,
     &     idxr12,idxc12,idxrba,idxcba,explicit)
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
*     Modified for R12 capability, GWR April 2007.
*
*----------------------------------------------------------------------*
      implicit none

      integer, parameter ::
     &     ntest = 100

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_operator.h'
      include 'def_contraction.h'
      include 'def_filinf.h'
      include 'ifc_operators.h'
      include 'ifc_input.h'
      include 'def_formula.h'
      include 'par_formnames_gen.h'

      type(formula), intent(inout) ::
     &     form_cclag

      integer, intent(in) ::
     &     nops,
     &     idxham,idxlag,idxtop,idxr12,idxc12,idxrba,idxcba,idxecc

      type(operator), intent(in) ::
     &     ops(nops)

      logical, intent(in) ::
     &     explicit

      ! local variables
      logical ::
     &     init_pn, init_pc, len, init

      type(contraction) ::
     &     contr

      integer ::
     &     lucclag,
     &     maxvtx, maxarc, maxexc, maxcnt, ncommin, ncommax,
     &     iloccls, ihoccls, icomm,
     &     idx, idxh, nterms, ncterm(5),
     &     npart, iexc, ihdh, ihdp, nlop, nlhc, iop, idxarc, ivtxoff,
     &     maxl, lagocc, rbaocc, narc, ansatze, ians, nextern, modcom
      ! occupations:
      integer ::
     &     iocc_l(ngastp,2), iocc_h(ngastp,2),
     &     iocc_hd(ngastp,2), iocc_hx(ngastp,2),
     &     iocc_hxovl(ngastp,2), iocc_ttot(ngastp,2),
     &     iocc_ltc(ngastp,2), iocc_scr(ngastp,2),
     &     i,j,k,l
      
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

      ! Allocatable arrays.
      integer, allocatable ::
     &     iocc_cba(:,:),iocc_rbcbov(:,:),iocc_hxovc(:,:),
     &     iocc_temp(:,:),iocc_rcco(:,:),iocc_r12(:,:),iocc_c12(:,:),
     &     iocc_r12c12(:,:), iocc_rtemp(:,:), itrip(:), itrip_part(:,:)

      ! external functions
      logical, external ::
     &     iocc_equal,
     &     next_part_number, next_part_pair, next_part_triple
      integer, external ::
     &     iopen_nus,
     &     ielsqsum, ielsum, iblk_occ, maxxlvl_op, ifndmin, ifndmax,
     &     ieqfac


      if (ntest.eq.100) then
        write(luout,*) '==============================='
        write(luout,*) ' output from set_cc_lagrangian'
        write(luout,*) '==============================='
        write(luout,*) 'idxecc,idxham,idxlag,idxtop:'
        write(luout,*) idxecc,idxham,idxlag,idxtop
        write(luout,*) 'idxr12,idxc12,idxrba,idxcba,explicit:'
        write(luout,*) idxr12,idxc12,idxrba,idxcba,explicit
      end if

      call atim_csw(cpu0,sys0,wall0)

      contr%idx_res = idxecc
      contr%iblk_res = 1

      ! L + H + 4times T gives maximum 6 vertices
      maxvtx = 6
      ! A contraction arc is the contraction matrix between each 
      ! operator.We can have at most 8 contraction arcs in CC, but
      ! up to 11 in CC-R12. (Check this limit, it is at least 11).
      maxarc=8
      if(explicit)then
        maxarc=11
      endif
      ! For convenience, we allocate the maximum number here.
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

c dbg
      write(luout,'("unit")')
      print *,form_cclag%fhand%unit
c dbg


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

      maxl=ops(idxlag)%n_occ_cls
      lagocc=maxl
      if(explicit)then
        rbaocc=ops(idxrba)%n_occ_cls
        maxl=maxl+rbaocc
        allocate(iocc_r12(ngastp,2),iocc_c12(ngastp,2),
     &       iocc_r12c12(ngastp,2),iocc_hxovc(ngastp,2),
     &       iocc_temp(ngastp,2),iocc_rcco(ngastp,2),
     &       iocc_rtemp(ngastp,2),itrip(3),itrip_part(3,maxpart))
      endif  
      l_loop: do iloccls = 0, maxl

c dbg
            write(luout,'("unit top l")')
            print *,form_cclag%fhand%unit
c dbg


        ncterm(1:5) = 0
        contr%nvtx=0
        narc=0

        ! we start with zero --
        ! here, we dump how to obtain the CC-energy
        if (iloccls.eq.0) then
          iocc_l(1:ngastp,1:2)=0
          nlop = 0
        elseif(iloccls.gt.0.and.iloccls.le.lagocc)then
        ! Copy L terms from operator list to local arrays.  
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
        elseif(explicit.and.iloccls.gt.lagocc.and.
     &         iloccls.le.(lagocc+rbaocc))then
        ! If R12 is requested then must deal with the terms in the 
        ! Lagrangian arising from Rbar. Must first contract the Rbar 
        ! term with the coefficient operator, Cbar.
          if(iloccls.eq.(lagocc+1))then
            allocate(iocc_cba(ngastp,2))
            allocate(iocc_rbcbov(ngastp,2))
          endif
        !  Copy Rbar and Cbar to local arrays. Check to ensure that 
        !  both are adjointed or not, they cannot be a mix.
          if(ops(idxrba)%dagger.and.ops(idxcba)%dagger)then
            iocc_l(1:ngastp,2)=
     &           ops(idxrba)%ihpvca_occ(1:ngastp,1,iloccls-lagocc)
            iocc_l(1:ngastp,1)=
     &           ops(idxrba)%ihpvca_occ(1:ngastp,2,iloccls-lagocc)
            iocc_cba(1:ngastp,2)=
     &           ops(idxcba)%ihpvca_occ(1:ngastp,1,1)
            iocc_cba(1:ngastp,1)=
     &           ops(idxcba)%ihpvca_occ(1:ngastp,2,1)
          elseif((.not.ops(idxrba)%dagger).and.
     &           (.not.ops(idxcba)%dagger))then
            iocc_l(1:ngastp,1:2)=
     &           ops(idxrba)%ihpvca_occ(1:ngastp,1:2,iloccls-lagocc)
            iocc_cba(1:ngastp,1:2)=
     &           ops(idxcba)%ihpvca_occ(1:ngastp,1:2,1)
          else
            call quit(1,'set_cc_lagrangian','Only 1 of C+ and R+ dagg.')
          endif
          contr%vertex(1)%idx_op=idxrba
          contr%vertex(1)%iblk_op=iloccls-lagocc
          contr%vertex(2)%idx_op=idxcba
          contr%vertex(2)%iblk_op=1
          contr%nvtx=2
          nlop=2

        ! Contract Cbar and Rbar, then form their resultants.
          iocc_rbcbov=iocc_overlap(iocc_l,.false.,iocc_cba,.false.)
          iocc_l=iocc_add(1,iocc_l,.false.,-1,iocc_rbcbov,.false.)
          iocc_cba=iocc_add(1,iocc_cba,.false.,-1,iocc_rbcbov,.true.)
  
          narc=1
          contr%narc=1
          contr%arc(1)%link(1)=1
          contr%arc(1)%link(2)=2
          contr%arc(1)%occ_cnt=iocc_rbcbov
        end if

        ! current projection is zeroth order space? find out here ....

        ! loop over blocks of Hamiltonian

        h_loop: do ihoccls = 1, ops(idxham)%n_occ_cls

c dbg
            write(luout,'("unit top h")')
            print *,form_cclag%fhand%unit
c dbg

          
          contr%narc = narc
          contr%nvtx = nlop+1
          ! Extract Hamiltonian block and place it into the contraction.
          iocc_h(1:ngastp,1:2) =
     &         ops(idxham)%ihpvca_occ(1:ngastp,1:2,ihoccls)

c dbg
            write(luout,'("unit top arcs 2")')
            print *,form_cclag%fhand%unit
c dbg

            write(luout,*) 'idxham,ihoccls,contr%nvtx',
     &           idxham,ihoccls,contr%nvtx

          contr%vertex(contr%nvtx)%idx_op=idxham
          contr%vertex(contr%nvtx)%iblk_op=ihoccls
          idxh = contr%nvtx ! remember H vertex number

c dbg
            write(luout,'("unit top arcs 3")')
            print *,form_cclag%fhand%unit
c dbg


          ! get excitation part of H
          iocc_hx = iocc_xdn(1,iocc_h)

          if(iloccls.le.lagocc)then

c dbg
            write(luout,'("unit top arcs")')
            print *,form_cclag%fhand%unit
c dbg

          ! get overlap with L+ ...
            iocc_hxovl = iocc_overlap(iocc_h,.false.,iocc_l,.true.)

          ! and test whether it is identical with excitation part of H
          ! (as we know that only T operators follow, so excitation part
          ! of H must fully contract with L)
            if (.not.iocc_equal(iocc_hx,.false.,iocc_hxovl,.false.))
     &           cycle h_loop

          ! so, excitation part of H defines the contraction with L+
          ! (unless, it is the 0-contraction)
            if (ielsqsum(iocc_hx,ngastp*2).gt.0) then
              nlhc = 1
              contr%narc = contr%narc+1
              contr%arc(1)%link(1)=1
              contr%arc(1)%link(2)=2
              contr%arc(1)%occ_cnt=iocc_dagger(iocc_hx)
            else
              nlhc = 0
            end if

c dbg
            write(luout,'("unit bottom 2")')
            print *,form_cclag%fhand%unit
c dbg


          elseif(explicit.and.iloccls.gt.lagocc)then
            ! Overlap Hx with R+ and C+.
            iocc_hxovl=iocc_overlap(iocc_h,.false.,iocc_l,.true.)
            iocc_hxovc=iocc_overlap(iocc_h,.false.,iocc_cba,.true.)

            ! Test whether Hx has fully contracted with L+ and C+. If 
            ! not, cycle over the H indices.
            iocc_temp=iocc_add(1,iocc_hxovl,.true.,
     &           1,iocc_hxovc,.true.)

            if(.not.iocc_equal(iocc_hx,.false.,iocc_temp,.true.))
     &           cycle h_loop

            ! Place the (up to) two arcs into the total contraction.
            nlhc=0
            if(ielsqsum(iocc_hx,ngastp*2).gt.0)then
              if(ielsqsum(iocc_hxovl,ngastp*2).gt.0)then
                nlhc=nlhc+1
                contr%narc=contr%narc+1
                contr%arc(narc+nlhc)%link(1)=1
                contr%arc(narc+nlhc)%link(2)=3
                contr%arc(narc+nlhc)%occ_cnt=iocc_dagger(iocc_hxovl)
              endif  
              if(ielsqsum(iocc_hxovc,ngastp*2).gt.0)then
                nlhc=nlhc+1
                contr%narc=contr%narc+1
                contr%arc(narc+nlhc)%link(1)=2
                contr%arc(narc+nlhc)%link(2)=3
                contr%arc(narc+nlhc)%occ_cnt=iocc_dagger(iocc_hxovc)
              endif  
              if(nlhc.eq.0)then
                call quit(1,'set_cc_lagrangian','R12+ H contraction.')
              endif
c              cycle h_loop
            endif  
c            cycle h_loop
          endif  

c dbg
            write(luout,'("unit bottom arcs")')
            print *,form_cclag%fhand%unit
c dbg


          ! Rest of H is the de-excitation part.
          iocc_hd = iocc_add(1,iocc_h,.false.,-1,iocc_hx,.false.)

          ! Get number of T that can contract with H (=rank 
          ! of commutator). Simply the number of free lines in Hd.
          ncommax = ielsum(iocc_hd,ngastp*2)
          if (ncommax.eq.0) then

            ! Excitation part of H must be identical with L+ (or with
            ! the contracted R+/C+ in an R12 calculation.) if there
            ! are no following T/R operators.
            if(iloccls.le.lagocc)then
              if (.not.iocc_equal(iocc_hx,.false.,iocc_l,.true.))
     &             cycle h_loop
            elseif(explicit.and.iloccls.gt.lagocc)then
              iocc_rcco=iocc_add(1,iocc_l,.false.,1,iocc_cba,.false.)
              if(.not.iocc_equal(iocc_hx,.false.,iocc_rcco,.true.))
     &             cycle h_loop
            endif  

            ! set up contraction info (only factor is missing)
            contr%fac = 1d0

c dbg
            write(luout,'("unit top print")')
            print *,form_cclag%fhand%unit
c dbg


            if (ntest.ge.100) then
              call prt_contr(luout,contr,ops)
            end if
            call wrt_contr(lucclag,contr)
            nterms = nterms+1
            ncterm(1) = ncterm(1)+1

c dbg
            write(luout,'("unit bottom print")')
            print *,form_cclag%fhand%unit
c dbg


          else

            ! total T occupation:
            ! [Ttotal] = [L^\dag] - [Hx] + [Hd^\dag]
            iocc_ttot = iocc_add(1,iocc_l,.true.,-1,iocc_hx,.false.)
            iocc_ttot = iocc_add(1,iocc_ttot,.false.,1,iocc_hd,.true.)
            ! Add in C+ if doing an R12 amplitude (L+=R+).
            if(explicit.and.iloccls.gt.lagocc)then
              iocc_ttot=iocc_add(1,iocc_ttot,.false.,1,iocc_cba,.true.)
            endif  

            !  check whether [Ttotal] looks the way we expect:
            iexc = iocc_ttot(ipart,1) ! keep excitation level in mind
            if(explicit)then
              iexc=iexc+iocc_ttot(iextr,1)
            endif  

            if (.not.(iocc_ttot(ihole,2).eq.iexc.and.
     &           iocc_ttot(ihole,1).eq.0.and.
     &           iocc_ttot(ipart,2).eq.0.and.
     &           iocc_ttot(iextr,2).eq.0)) then
              write(luout,*) 'fishy [Ttot]:'
              call wrt_occ(luout,iocc_ttot)
              call quit(1,'set_cc_lagrangian','fishy occupation')
            end if
            ! Some extra conditions required for R12 calculations.   
            nextern=iocc_ttot(iextr,1)
            if(explicit)then
              call get_argument_value('method.R12','ansatz',
     &             ival=ansatze)
              if(nextern.ne.0)then
                if(iocc_ttot(ihole,2).lt.2)
     &               cycle h_loop
                if(iocc_ttot(ipart,1).eq.0.and.
     &               mod(nextern,2).ne.0) cycle h_loop
                if(ansatze.eq.1)then
                  if(mod(nextern,2).ne.0)cycle h_loop
                endif
              endif  
            endif

            if(nextern.gt.0)then
              write(luout,'("T tot")')
              write(luout,'(4i4)') iocc_ttot(1:ngastp,1)
              write(luout,'(4i4)') iocc_ttot(1:ngastp,2)
            endif  

            ncommin = iocc_ttot(ipart,1)/(maxexc+1) + 1
            ! Account for external orbitals if R12 requested.
            if(explicit.and.iocc_ttot(iextr,1).gt.0)then
              if(iocc_ttot(ipart,1).eq.0)then
                ncommin=nextern/2
              else
                if(iocc_ttot(ipart,1).eq.1)then
                  if(nextern.eq.1)then
                    ncommin=1
                  elseif(nextern.gt.1)then
                    ncommin=1+iocc_ttot(iextr,1)/2
                  endif
                elseif(iocc_ttot(ipart,1).gt.1)then
                  if(nextern.eq.1)then
                    ncommin=2+(iocc_ttot(ipart,1)-1)/(maxexc+1)
                  elseif(nextern.gt.1)then
                    ncommin=ncommin+nextern/2+mod(nextern,2)
                  endif  
                endif  
              endif 
            endif  

            if (ntest.ge.100) write(luout,*) 'ncommin, ncommax: ',
     &           ncommin, ncommax

            ! Loop over possible commutators of T (and R if R12). 
            ! Initially treat the single commutator terms, including
            ! the R12 terms, then the multimple commutators of T, 
            ! and finally the multiple commutators with R terms.
            com_loop: do icomm = ncommin, ncommax

c dbg
            write(luout,'("unit top")')
            print *,form_cclag%fhand%unit
c dbg


              if (icomm.eq.1) then
                if(nextern.eq.0)then
                  contr%fac = 1d0
                  contr%nvtx = nlop+2
                  contr%vertex(contr%nvtx)%idx_op = idxtop
                  idx = iblk_occ(iocc_ttot,.false.,ops(idxtop))
                  if (idx.le.0) then
                    write(luout,*) 'T occupation not found in list:'
                    write(luout,*) iocc_ttot(1:ngastp,1)
                    write(luout,*) iocc_ttot(1:ngastp,2)
                    stop 'occupation not found in list (1)'
                  end if
                  contr%vertex(contr%nvtx)%iblk_op = idx

                  ! [Hd] defines contraction between H and T
                  contr%narc = contr%narc+1
                  contr%arc(contr%narc)%link(1)=idxh
                  contr%arc(contr%narc)%link(2)=contr%nvtx
                  contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)
     &                 =iocc_hd(1:ngastp,1:2)

                  ! [T^\dag]-[Hd] defines contraction between L and T
                  iocc_ltc = iocc_add(-1,iocc_hd,.false.,
     &                 +1,iocc_ttot,.true.)
                  if (ielsum(iocc_ltc,ngastp*2).gt.0) then
                    if(iloccls.le.lagocc)then
                      contr%narc=contr%narc+1
                      contr%arc(contr%narc)%link(1)=1
                      contr%arc(contr%narc)%link(2)=contr%nvtx
                      contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)=
     &                     iocc_ltc(1:ngastp,1:2)
                    elseif(explicit.and.iloccls.gt.lagocc)then
                      if(iocc_ltc(ipart,2).ne.0.or.
     &                     iocc_ltc(iextr,2).ne.0)then
                        contr%narc=contr%narc+1
                        contr%arc(contr%narc)%link(1)=1
                        contr%arc(contr%narc)%link(2)=contr%nvtx
                        contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)=0
                        contr%arc(contr%narc)%occ_cnt(ipart,2)=
     &                       iocc_ltc(ipart,2)
                        contr%arc(contr%narc)%occ_cnt(iextr,2)=
     &                       iocc_ltc(iextr,2)
                      endif
                      if(iocc_ltc(ihole,1).ne.0)then
                        contr%narc=contr%narc+1
                        contr%arc(contr%narc)%link(1)=2
                        contr%arc(contr%narc)%link(2)=contr%nvtx
                        contr%arc(contr%narc)%occ_cnt(1:ngastp,1:2)=0
                        contr%arc(contr%narc)%occ_cnt(ihole,1)=
     &                       iocc_ltc(ihole,1)
                      endif  
                    endif
                  end if

                ! Add in contractions with R12 operators if needed.
                else
                  contr%fac=1d0
                  contr%nvtx=nlop+3
                  contr%vertex(contr%nvtx)%idx_op=idxr12
                  idx=iblk_occ(iocc_ttot,.false.,ops(idxr12))
                  contr%vertex(contr%nvtx+1)%idx_op=idxc12
                  if(idx.le.0)then
                    write(luout,*) 'R occupation not found in list:'
                    write(luout,*) iocc_ttot(1:ngastp,1)
                    write(luout,*) iocc_ttot(1:ngastp,2)
                    stop 'occupation not found in list (1)'
                  end if                    
                  ! Must first contract R and C.
                  iocc_r12(1:ngastp,1:2)=
     &                 ops(idxr12)%ihpvca_occ(1:ngastp,1:2,idx)
                  iocc_c12(1:ngastp,1:2)=
     &                 ops(idxc12)%ihpvca_occ(1:ngastp,1:2,1)
                  
                  contr%vertex(nlop+2)%idx_op=idxr12
                  contr%vertex(nlop+2)%iblk_op=idx
                  contr%vertex(nlop+3)%idx_op=idxc12
                  contr%vertex(nlop+3)%iblk_op=1

                  iocc_r12c12=
     &               iocc_overlap(iocc_r12,.false.,iocc_c12,.false.)
                  iocc_r12=
     &               iocc_add(1,iocc_r12,.false.,-1,iocc_r12c12,.false.)
                  iocc_c12=
     &               iocc_add(1,iocc_c12,.false.,-1,iocc_r12c12,.false.)
                  contr%narc=contr%narc+1
                  contr%arc(contr%narc)%link(1)=nlop+2
                  contr%arc(contr%narc)%link(2)=nlop+3
                  contr%arc(contr%narc)%occ_cnt=iocc_r12c12

                  ! Now contract H with both new R and new C. Hole 
                  ! lines of Hd define the contraction with C and 
                  ! particle/external lines define that with R.
                  if(iocc_hd(ipart,2).gt.0.or.iocc_hd(iextr,2).gt.0)then
                    iocc_temp(1:ngastp,1:2)=0
                    iocc_temp(ipart,2)=iocc_hd(ipart,2)
                    iocc_temp(iextr,2)=iocc_hd(iextr,2)
                    contr%narc=contr%narc+1
                    contr%arc(contr%narc)%link(1)=nlop+1
                    contr%arc(contr%narc)%link(2)=nlop+2
                    contr%arc(contr%narc)%occ_cnt=iocc_temp
                  endif
                  if(iocc_hd(ihole,1).gt.0)then
                    iocc_temp(1:ngastp,1:2)=0
                    iocc_temp(ihole,1)=iocc_hd(ihole,1)                    
                    contr%narc=contr%narc+1
                    contr%arc(contr%narc)%link(1)=nlop+1
                    contr%arc(contr%narc)%link(2)=nlop+3
                    contr%arc(contr%narc)%occ_cnt=iocc_temp
                  endif
                  ! Finally, contract remaining lines with L+ (or R+/C+)
                  ! [T^\dag]-[Hd] defines contraction between L and T
                  iocc_ltc = iocc_add(-1,iocc_hd,.false.,
     &                 +1,iocc_ttot,.true.)
                  ! Hole part of [T+]-Hd defines contraction with Cbar
                  ! and particle/external part of the same defines 
                  ! the contraction with Rbar.
                  if(iocc_ltc(ipart,2).gt.0.or.
     &                 iocc_ltc(iextr,2).gt.0)then
                    iocc_temp(1:ngastp,1:2)=0
                    iocc_temp(ipart,2)=iocc_ltc(ipart,2)
                    iocc_temp(iextr,2)=iocc_ltc(iextr,2)
                    contr%narc=contr%narc+1
                    contr%arc(contr%narc)%link(1)=1
                    contr%arc(contr%narc)%link(2)=nlop+2
                    contr%arc(contr%narc)%occ_cnt=iocc_temp
                  endif  
                  if(iocc_ltc(ihole,1).gt.0)then
                    iocc_temp(1:ngastp,1:2)=0
                    iocc_temp(ihole,1)=iocc_ltc(ihole,1)
                    contr%narc=contr%narc+1
                    contr%arc(contr%narc)%link(1)=nlop
                    contr%arc(contr%narc)%link(2)=nlop+3
                    contr%arc(contr%narc)%occ_cnt=iocc_temp
                  endif   
                endif  
                if (ntest.ge.100) then
                  call prt_contr(luout,contr,ops)
                end if
                call wrt_contr(lucclag,contr)
                nterms = nterms+1
                ncterm(2) = ncterm(2)+1

              elseif(icomm.gt.1.)then
                if(nextern.eq.0)then
                  ! double and higher commutator term:
                  ! loop over all n-fold partitionings of [Ttotal] 
                  init_pn = .true.
                  do while (next_part_number(init_pn,.true.,iexc_part,
     &                 iexc,icomm,1,maxexc))
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
     &                   next_part_pair(init_pc,.false.,ihd_part,
     &                   ihd,icomm,1,maxcnt))
                      init_pc = .false.

                      write(luout,'("Hd Partition")')
                      write(luout,'(4i4)')ihd_part(1,1:icomm)
                      write(luout,'(4i4)')ihd_part(2,1:icomm)

                      ! check, whether all contractions are possible
                      do idx = 1, icomm
                        if (ihd_part(1,idx).gt.iexc_part(idx).or.
     &                       ihd_part(2,idx).gt.iexc_part(idx))
     &                       cycle pc_loop
                      end do
                    
                      ! represent each p/h pair for contractions
                      ! by one integer (for ordering purposes)
                      ihd_part_idx(1:icomm) = ihd_part(1,1:icomm)
     &                     + ihd_part(2,1:icomm)*(maxcnt+1)

                      ! for equivalent T operators, only one of
                      ! the possible permutations of contractions
                      ! needs be generated: allow only the contr.
                      ! with increasing order of ihd_part_idx
                      do idx = 2, icomm
                        if (iexc_part(idx-1).eq.iexc_part(idx).and.
     &                       ihd_part_idx(idx-1).gt.ihd_part_idx(idx))
     &                       cycle pc_loop
                      end do

                      ! set up contraction information
                      ! get factor from equivalent contractions
                      ! cf. my (andreas) notes on origin of factors
                      contr%fac =
     &                    1d0/dble(ieqfac(iexc_part,ihd_part_idx,icomm))
                    
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
                        contr%arc(idxarc)%occ_cnt(ihole,1)
     &                       = ihd_part(1,iop)
                        contr%arc(idxarc)%occ_cnt(ipart,2)
     &                       = ihd_part(2,iop)
                        
                        ! contraction between L and T
                        iocc_ltc = iocc_add(-1,
     &                       contr%arc(idxarc)%occ_cnt,.false.,
     &                       +1,iocc_scr,.true.)
                        
                        if (ielsum(iocc_ltc,ngastp*2).gt.0) then
                          if(iloccls.le.lagocc)then
                            idxarc = idxarc+1
                            contr%arc(idxarc)%link(1)=1
                            contr%arc(idxarc)%link(2)=nlop+1+iop
                            contr%arc(idxarc)%occ_cnt(1:ngastp,1:2)=
     &                           iocc_ltc(1:ngastp,1:2)
                          elseif(explicit.and.iloccls.gt.lagocc)then
                            if(iocc_ltc(ipart,2).gt.0)then
                              idxarc=idxarc+1
                              contr%arc(idxarc)%link(1)=1
                              contr%arc(idxarc)%link(2)=nlop+1+iop
                              contr%arc(idxarc)%occ_cnt(1:ngastp,1:2)=0
                              contr%arc(idxarc)%occ_cnt(ipart,2)=
     &                           iocc_ltc(ipart,2)
                            endif
                            if(iocc_ltc(ihole,1).gt.0)then
                              idxarc=idxarc+1
                              contr%arc(idxarc)%link(1)=2
                              contr%arc(idxarc)%link(2)=nlop+1+iop
                              contr%arc(idxarc)%occ_cnt(1:ngastp,1:2)=0
                              contr%arc(idxarc)%occ_cnt(ihole,1)=
     &                             iocc_ltc(ihole,1)
                            endif  
                          endif  
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

                  end do        ! [Ttotal] partitioning

                ! Multiple commutators involving R12 operators.  
                else
                  if(.not.explicit) call quit(1,'set_cc_lagrangian',
     &                 'External lines but R12 not requested.')
                  ans_loop: do ians=1,ansatze
                
                  ! Ascertain the form of the commutator by checking the
                  ! number of external lines present and the chosen 
                  ! ansatz.
                    iocc_rtemp(1:ngastp,1:2)=0
                    if(nextern.eq.1)then
                      iocc_rtemp(1:ngastp,1:2)=
     &                     ops(idxr12)%ihpvca_occ(1:ngastp,1:2,1)
                      iocc_temp=iocc_add(1,iocc_ttot,.false.,
     &                     -1,iocc_rtemp,.false.)
                      modcom=icomm-1
                    elseif(nextern.eq.2)then  
                      iocc_rtemp(1:ngastp,1:2)=
     &                     ops(idxr12)%ihpvca_occ(1:ngastp,1:2,ians)
                      iocc_temp=iocc_add(1,iocc_ttot,.false.,
     &                     -1*(3-ians),iocc_rtemp,.false.)
                      modcom=icomm-(3-iocc_rtemp(iextr,1))
                      if(iocc_temp(ihole,2).lt.0)cycle ans_loop
                    elseif(nextern.eq.3)then
                      iocc_rtemp(1:ngastp,1:2)=
     &                     ops(idxr12)%ihpvca_occ(1:ngastp,1:2,1)
                      iocc_temp=iocc_add(1,iocc_ttot,.false.,
     &                     -1,iocc_rtemp,.false.)
                      modcom=0
                    else
                      iocc_rtemp(1:ngastp,1:2)=
     &                     ops(idxr12)%ihpvca_occ(1:ngastp,1:2,2)
                      iocc_temp=iocc_add(1,iocc_ttot,.false.,
     &                     -2,iocc_rtemp,.false.)  
                      modcom=0
                    endif  

c                    write(luout,'("T-R")')
c                    write(luout,'(4i4)')iocc_temp(1:ngastp,1)
c                    write(luout,'(4i4)')iocc_temp(1:ngastp,2)
c                    cycle h_loop

                    ! Contract the first R operator with C.
                    iocc_c12(1:ngastp,1:2)=
     &                   ops(idxc12)%ihpvca_occ(1:ngastp,1:2,1)
                    contr%nvtx=nlop+3
                    idxarc=contr%narc
                    ivtxoff=nlop+1

                    idx=iblk_occ(iocc_rtemp,.false.,ops(idxr12))
                    iocc_r12c12=
     &                   iocc_overlap(iocc_rtemp,.false.,
     &                   iocc_c12,.false.)
                    idxarc=idxarc+1
                    contr%vertex(ivtxoff+1)%idx_op=idxr12
                    contr%vertex(ivtxoff+1)%iblk_op=idx
                    contr%vertex(ivtxoff+2)%idx_op=idxc12
                    contr%vertex(ivtxoff+2)%iblk_op=1
                    contr%arc(idxarc)%link(1)=ivtxoff+1
                    contr%arc(idxarc)%link(2)=ivtxoff+2
                    contr%arc(idxarc)%occ_cnt=iocc_r12c12
                    ivtxoff=ivtxoff+2

                    ! If more than one R operator required, contract it
                    ! with C here.
                    if(nextern.eq.3)then
                      iocc_rtemp(1:ngastp,1:2)=iocc_temp(1:ngastp,1:2)
                    endif  
                    if(modcom.eq.(icomm-2))then
                      contr%nvtx=contr%nvtx+2
                      idxarc=idxarc+1
                      idx=iblk_occ(iocc_rtemp,.false.,ops(idxr12))
                      iocc_r12c12=
     &                     iocc_overlap(iocc_rtemp,.false.,
     &                     iocc_c12,.false.)
                      contr%vertex(ivtxoff+1)%idx_op=idxr12
                      contr%vertex(ivtxoff+1)%iblk_op=idx
                      contr%vertex(ivtxoff+2)%idx_op=idxc12
                      contr%vertex(ivtxoff+2)%iblk_op=1                      
                      contr%arc(idxarc)%link(1)=ivtxoff+1
                      contr%arc(idxarc)%link(2)=ivtxoff+2
                      contr%arc(idxarc)%occ_cnt=iocc_r12c12
                      ivtxoff=ivtxoff+2
                    endif  

                    ! Partition the de-excitation part of H into up to
                    ! 4 triplets (4 comms. of hole/part/extr).
                    itrip(1)=iocc_hd(ihole,1)
                    itrip(2)=iocc_hd(ipart,2)
                    itrip(3)=iocc_hd(iextr,2)
                    maxcnt = max(2,ifndmax(iexc_part,1,icomm,1))
c                    write(luout,'(3i4)')itrip(1:3)
                    
                    init_pc=.true.
                    do while (next_part_triple(init_pc,.false.,
     &                   itrip_part,itrip,icomm,1,maxcnt))
                      init_pc=.false.
                     
                    enddo  
                    cycle ans_loop
  
                  enddo  ans_loop
                endif  
              end if          ! number of commutators switch
                
c dbg
              write(luout,'("unit bottom")')
              print *,form_cclag%fhand%unit
c dbg


            end do com_loop   ! loop over icomm

          end if ! no-commutator-at-all exception

        end do h_loop

c        write(luout,'(4x,i4,3x,iluout,x,5(x,i4))')
c dbg -- add "??" mark for grepping
        write(luout,'(2x,"??",i4,3x,i6,x,5(x,i4))')
     &       iloccls,sum(ncterm(1:5)),ncterm(1:5)

      end do l_loop
      write(luout,'(2x,42("-"))')

c dbg
      print *,form_cclag%fhand%unit
c dbg
      call file_close_keep(form_cclag%fhand)
      deallocate(contr%vertex,contr%arc)
      if(explicit)then
        deallocate(iocc_cba,iocc_rbcbov,iocc_hxovc,iocc_temp,iocc_rcco,
     &       iocc_r12,iocc_c12,iocc_r12c12,iocc_rtemp,itrip,itrip_part)
      endif  

      call atim_csw(cpu,sys,wall)
      write(luout,*) 'Number of generated terms: ',nterms
      call prtim(luout,'CC Lagrangian',cpu-cpu0,sys-sys0,wall-wall0)

      end
