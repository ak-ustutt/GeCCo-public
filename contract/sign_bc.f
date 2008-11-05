*----------------------------------------------------------------------*
      subroutine sign_bc(bc_sign,
     &     idxsuper1,idxsuper2,svertex,nvtx,
     &     occ_vtx,iocc_op1,iocc_op2,iocc_cnt,
     &     njoined_op1,njoined_op2,njoined_op1op2,njoined_cnt,
     &     merge_ex1cnt,merge_ex2cnt,merge_ex1ex2,
     &     ld_m1c,ld_m2c,ld_m12)
*----------------------------------------------------------------------*
*     calculate sign for binary multi-component contraction
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'

      integer, parameter ::
     &     ntest = 00

      logical ::
     &     first
      integer, intent(in) ::
     &     nvtx,
     &     idxsuper1, idxsuper2, svertex(nvtx),
     &     njoined_op1,njoined_op2,njoined_op1op2,njoined_cnt,
     &     occ_vtx(ngastp,2,nvtx),
     &     ld_m1c, ld_m2c, ld_m12, 
     &     iocc_op1(ngastp,2,njoined_op1),
     &     iocc_op2(ngastp,2,njoined_op2),
     &     iocc_cnt(ngastp,2,njoined_cnt),
     &     merge_ex1cnt(ld_m1c,2,*),
     &     merge_ex2cnt(ld_m2c,2,*),
     &     merge_ex1ex2(ld_m12,2,*)
      real(8), intent(out) ::
     &     bc_sign

      integer ::
     &     iop1op2, ivtx1, ivtx2, idx, ivtx, ica, hpvx, hpvx2,
     &     ivtx1m, ivtx2m, icamod,
     &     ivtx1raw, ivtx2raw, idx12,
     &     ncntc, ncnta, nex1c, nex1a, nex2c, nex2a, nencl, nrem,
     &     icasign, ihpvxsign, icnt
      integer ::
     &     iocc_prim(ngastp,2,nvtx),
     &     idx_ex1ex2(2,nvtx)

      integer, external ::
     &     idxlist

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,'this is sign_bc')
        write(luout,*) 'maps:'
        call iwrtma(merge_ex1cnt,ld_m1c,2*njoined_op1,
     &                           ld_m1c,2*njoined_op1)
        call iwrtma(merge_ex2cnt,ld_m2c,2*njoined_op2,
     &                           ld_m2c,2*njoined_op2)
        call iwrtma(merge_ex1ex2,ld_m12,2*njoined_op1op2,
     &                           ld_m12,2*njoined_op1op2)
      end if

      iocc_prim = occ_vtx

      ! loop over svertex and assemble occupations of
      ! primitive vertices in original order:
      iop1op2 = 0
      ivtx1 = 0
      ivtx2 = 0
      do idx = 1, nvtx
        if (svertex(idx).eq.idxsuper1) then
          iop1op2 = iop1op2+1
          ivtx1 = ivtx1+1
c          iocc_prim(1:ngastp,1:2,iop1op2) = iocc_op1(1:ngastp,1:2,ivtx1)
          idx_ex1ex2(1,ivtx1) = idx !iop1op2
c          idx_ex1ex2(1,iop1op2) = 1
c          idx_ex1ex2(1,iop1op2) = ivtx1
        else if (svertex(idx).eq.idxsuper2) then
          iop1op2 = iop1op2+1
          ivtx2 = ivtx2+1
c          iocc_prim(1:ngastp,1:2,iop1op2) = iocc_op2(1:ngastp,1:2,ivtx2)
          idx_ex1ex2(2,ivtx2) = idx !iop1op2
c          idx_ex1ex2(1,iop1op2) = 2
c          idx_ex1ex2(1,iop1op2) = ivtx2
        end if
      end do
c      nex_pvtx = iop1op2
c dbg
      if (iop1op2.ne.njoined_op1+njoined_op2) stop 'what is this?'
c dbg
      if (ntest.ge.100) then
        write(luout,*) 'Contractions:'
        call wrt_occ_n(luout,iocc_cnt,njoined_cnt)
        write(luout,*) 'Vertices in original order'
        call wrt_occ_n(luout,iocc_prim,nvtx)
      end if

      ! 
      ! loop over primitive contractions

      icasign = 0
      ihpvxsign = 0
      do icnt = 1, njoined_cnt
        ! get indices of involved primitive vertices
        ivtx1 = 0
        do ivtx = 1, njoined_op1
          if (idxlist(icnt,merge_ex1cnt(1,1,ivtx),ld_m1c,1).gt.0) then
            ivtx1 = idx_ex1ex2(1,ivtx)
            exit
          end if
        end do
        ! dto for vertex 2
        ivtx2 = 0
        do ivtx = 1, njoined_op2
          if (idxlist(icnt,merge_ex2cnt(1,1,ivtx),ld_m2c,1).gt.0) then
            ivtx2 = idx_ex1ex2(2,ivtx)
            exit
          end if
        end do
        if (ntest.ge.100) write(luout,*) 'CNT #',icnt,': ',ivtx1,ivtx2
c dbg
        if (ivtx1*ivtx2.eq.0) stop 'sbc: oha!'
c dbg
        ! CNT and sequence of vertices refers to sequence of super-vertices
        ! the sign rules, however, refer to the original sequence of
        ! vertices; so, if ivtx1>ivtx2, we have to switch things around:
        if (ivtx1.gt.ivtx2) then
          ncntc = sum(iocc_cnt(1:ngastp,2,icnt))
          ncnta = sum(iocc_cnt(1:ngastp,1,icnt))
        else
          ! usual behaviour:
          ncntc = sum(iocc_cnt(1:ngastp,1,icnt))
          ncnta = sum(iocc_cnt(1:ngastp,2,icnt))
        end if

        ! remove CNT from vertices
        iocc_prim(1:ngastp,1:2,ivtx1) =
     &       iocc_prim(1:ngastp,1:2,ivtx1) - iocc_cnt(1:ngastp,1:2,icnt)
        iocc_prim(1:ngastp,1,ivtx2) =
     &       iocc_prim(1:ngastp,1,ivtx2) - iocc_cnt(1:ngastp,2,icnt)
        iocc_prim(1:ngastp,2,ivtx2) =
     &       iocc_prim(1:ngastp,2,ivtx2) - iocc_cnt(1:ngastp,1,icnt)

        ! ensure increasing sequence: see above
        ivtx1m = min(ivtx1,ivtx2)
        ivtx2m = max(ivtx1,ivtx2)

        ! number of indices on enclosed operators
        nencl = 0
        do ivtx = ivtx1m+1, ivtx2m-1
          nencl = nencl + sum(iocc_prim(1:ngastp,1:2,ivtx))
        end do

        ! n(CNT(A)) * n(enclosed)
        icasign = mod(icasign + ncnta*nencl,2)

        ! enclosed + remaining on contracted vertices
        nrem = sum(iocc_prim(1:ngastp,1:2,ivtx1m))+
     &         sum(iocc_prim(1:ngastp,1:2,ivtx2m))
        ! n(CNT(C)) * [ n(enclosed) + n(vertices remaining)]
        icasign = mod(icasign + ncntc*(nencl+nrem),2)

        if (ntest.ge.100) then
          write(luout,*) 'ncntc,ncnta,nencl,nrem: ',
     &         ncntc,ncnta,nencl,nrem
          write(luout,*) 'updated CA sign:   ',icasign
        end if

        do ica = 1, 2
          icamod = ica
          ! see above
          if (ivtx1.gt.ivtx2) icamod = 3-ica
          do hpvx = 2, ngastp
            if (iocc_cnt(hpvx,icamod,icnt).gt.0) then
              do hpvx2 = 1, hpvx-1
                ihpvxsign = mod(ihpvxsign 
     &         +iocc_cnt(hpvx,icamod,icnt)*iocc_prim(hpvx2,ica,ivtx1m)
     &         +iocc_cnt(hpvx,icamod,icnt)*iocc_prim(hpvx2,3-ica,ivtx2m)
     &                          , 2)
              end do
            end if
          end do
        end do

        if (ntest.ge.100) then
          write(luout,*) 'updated HPVX sign: ',ihpvxsign
          write(luout,*) 'updated vertices:'
          call wrt_occ_n(luout,iocc_prim,nvtx)
        end if

      end do

      ! loop over primitive vertices of final array
      do iop1op2 = 1, njoined_op1op2

        first = .true.
c        ! we merge on the first vertex
c        ivtx1raw = merge_ex1ex2(1,1,iop1op2)
c
c        ivtx1 = idx_ex1ex2(1,ivtx1raw)

        ! get the parties to be merged:
        do idx12 = 1, 2

          do idx = 1, ld_m12

            ivtx2raw = merge_ex1ex2(idx,idx12,iop1op2)

            if (ivtx2raw.eq.0) exit

            if (first) then
              ivtx1 = idx_ex1ex2(idx12,ivtx2raw)
              first = .false.
              cycle
            else
              ivtx2 = idx_ex1ex2(idx12,ivtx2raw)
            end if

c            if (ivtx1.eq.ivtx2) cycle

            if (ntest.ge.100) write(luout,*) 'merging: ',ivtx1,ivtx2
            
            ! same story as above: in case of supervertices,
            ! ivtx1>ivtx2 may happen so we do the following
            ! fix: we want to merge on ivtx1, but for the sign
            !      we need ivtx1 < ivtx2, so:
            ivtx1m = min(ivtx1,ivtx2)
            ivtx2m = max(ivtx1,ivtx2)

            ! get n(EX1(A)), n(EX2(C)), n(EX2(A))
            nex1c = sum(iocc_prim(1:ngastp,1,ivtx1m))
            nex1a = sum(iocc_prim(1:ngastp,2,ivtx1m))
            nex2c = sum(iocc_prim(1:ngastp,1,ivtx2m))
            nex2a = sum(iocc_prim(1:ngastp,2,ivtx2m))
            ! count number of enclosed indices
            ! but only those which belong the presently contracting
            ! super-vertices (others are taken care of by sh_sign)
            ! I promise to come up with a better treatment of signs as soon
            ! as possible ...
            nencl = 0
            do ivtx = ivtx1m+1, ivtx2m-1
c new
              if (svertex(ivtx).ne.idxsuper1.and.
     &            svertex(ivtx).ne.idxsuper2) cycle
c new
              nencl = nencl + sum(iocc_prim(1:ngastp,1:2,ivtx))
            end do
            if (ivtx1.lt.ivtx2) then
              ! 1 2 sequence
              ! n(EX2(C))*(n(EX1(A))+nenclosed) + n(EX2(A))*nenclosed
              icasign = mod(icasign+nex2c*(nex1a+nencl)+nex2a*nencl,2)
            else
              ! 2 1 sequence
              ! n(EX1(C))*(n(EX2(C))+nenclosed) 
              ! + n(EX1(A))*(n(EX2(C))+n(EX2(A))+nenclosed)
              icasign = mod(icasign+nex1c*(nex2c+nencl)+
     &                             +nex1a*(nex2c+nex2a+nencl), 2)
            end if
            if (ntest.ge.100) then
              write(luout,'(1x,a,8i3)')
     &             'nencl,nex1a,nex1c,nex2a,nex2c: ',
     &              nencl,nex1a,nex1c,nex2a,nex2c
              write(luout,*) 'updated CA sign:   ',icasign
            end if

            ! hpvx transpositions:
            ! BUGFIX:
            ! after the CA carried out above, we have in all cases
            ! the ordering vtx1, vtx2, so use these (and not vtx1m/vtx2m)
            ! does not affect usual CC and linear R12 (as it seems)
            do hpvx = 2, ngastp
             ! ordering: ex1c,ex2c, but ex2a,ex1a
cold              if (iocc_prim(hpvx,1,ivtx1m).gt.0) then
              if (iocc_prim(hpvx,1,ivtx1).gt.0) then
c dbg
c                print *,'hpvx,ivtx1m,ivtx2m: ',hpvx,ivtx1m,ivtx2m
c                print *,'iocc1:',iocc_prim(hpvx,1,ivtx1m)
c                print *,'iocc2:',iocc_prim(1:ngastp,1,ivtx2m)
c dbg
                do hpvx2 = 1, hpvx-1
                  ihpvxsign = mod(ihpvxsign 
cold     &            +iocc_prim(hpvx,1,ivtx1m)*iocc_prim(hpvx2,1,ivtx2m),2)
     &            +iocc_prim(hpvx,1,ivtx1)*iocc_prim(hpvx2,1,ivtx2),2)
                end do
              end if
cold              if (iocc_prim(hpvx,2,ivtx2m).gt.0) then
              if (iocc_prim(hpvx,2,ivtx2).gt.0) then
                do hpvx2 = 1, hpvx-1
                  ihpvxsign = mod(ihpvxsign 
c     &            +iocc_prim(hpvx,2,ivtx2m)*iocc_prim(hpvx2,2,ivtx1m),2)
     &            +iocc_prim(hpvx,2,ivtx2)*iocc_prim(hpvx2,2,ivtx1),2)
                end do
              end if
            end do

            ! finally we merge the two:
            ! add to ivtx1 and remove from ivtx2
            ! here we need the ivtx1, ivtx2 in their original sequence:
            iocc_prim(1:ngastp,1:2,ivtx1) =
     &           iocc_prim(1:ngastp,1:2,ivtx1) +
     &           iocc_prim(1:ngastp,1:2,ivtx2)
            iocc_prim(1:ngastp,1:2,ivtx2) = 0

            if (ntest.ge.100) then
              write(luout,*) 'updated HPVX sign: ',ihpvxsign
              write(luout,*) 'updated vertices:'
              call wrt_occ_n(luout,iocc_prim,nvtx)
            end if
          end do
        end do
      end do

      if (mod(icasign+ihpvxsign,2).eq.0) bc_sign = 1d0
      if (mod(icasign+ihpvxsign,2).eq.1) bc_sign = -1d0

      if (ntest.ge.100) then
        write(luout,*) 'Final CA-sign:   ',(-1d0)**(dble(icasign))
        write(luout,*) 'Final HPVX-sign: ',(-1d0)**(dble(ihpvxsign))
        write(luout,*) 'Final total BC:  ',bc_sign
      end if

      end
