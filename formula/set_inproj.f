*----------------------------------------------------------------------*
      subroutine set_inproj(proto,occ_vtx,ok,
     &                      inproj,ninproj)
*----------------------------------------------------------------------*
*     set arcs on proto-contraction that result from inner projections;
*     the inproj array is explained in expand_op_product2
*     ok is set to true, if we think that gen_contr will be successful
*
*     initial version for unique cases
*
*     andreas, march 2008
*
*----------------------------------------------------------------------*

      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'
      include 'par_projectors.h'
      include 'ifc_operators.h'

      integer, parameter ::
     &     ntest =  00

      type(contraction), intent(inout) ::
     &     proto
      logical, intent(out) ::
     &     ok
      integer, intent(in) ::
     &     occ_vtx(ngastp,2,proto%nvtx), 
     &     ninproj, inproj(4,ninproj)

      integer ::
     &     narc_old, iproj, ivtx1, ivtx2, rank, type, type2, nidx,
     &     iblk, ihpvx, idx, ivtx, ica1, ica2
      integer ::
     &     occ1(ngastp,2), occ2(ngastp,2),
     &     ovl(ngastp,2), occ_cnt(ngastp,2)

      integer, external ::
     &     ielsum

      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'here comes set_inproj')
        write(lulog,*) 'inproj: '
        do idx = 1, ninproj
          write(lulog,'(2x,i3," - ",i3," r=",i3," t=",i3)')
     &         inproj(1:4,idx)
        end do
        call wrt_occ_n(lulog,occ_vtx,proto%nvtx)
      end if

      ok = .true.

      narc_old = proto%narc

      ! loop over requests
      do iproj = 1, ninproj
        ivtx1 = min(inproj(1,iproj),inproj(2,iproj))
        ivtx2 = max(inproj(1,iproj),inproj(2,iproj))
        rank  = inproj(3,iproj)
        type  = abs(inproj(4,iproj))
        type2 = 0
        if (inproj(4,iproj).ge.0) then
          ! the usual case: OCC1(A) connects to OCC2(C)
          ica1 = 2
          ica2 = 1
        else
          ! reversed case
          ica1 = 1
          ica2 = 2
        end if
        if (rank.le.0.or.type.lt.0.or.(rank.eq.2.and.type.gt.nproj))
     &       call quit(1,'set_inproj','watch your input')
        if (rank.eq.1) then
          ! NOTE: type2 does not work properly ...
          type2 = type/(ngastp+1)
          type = mod(type,ngastp+1)
        end if
c dbg
c        print *,'v1,v2,r,t: ',ivtx1,ivtx2,rank,type,type2
c dbg

        ! get A of first operator
        occ1 = 0
        occ1(1:ngastp,ica1) = occ_vtx(1:ngastp,ica1,ivtx1)
        ! get C of second operator
        occ2 = 0
        occ2(1:ngastp,ica2) = occ_vtx(1:ngastp,ica2,ivtx2)
        ovl = iocc_overlap(occ1,.false.,occ2,.true.)

c dbg
c        write(lulog,*) 'occ1, occ2'
c        call wrt_occ(lulog,occ1)
c        call wrt_occ(lulog,occ2)
c        write(lulog,*) 'ovl'
c        call wrt_occ(lulog,ovl)
c dbg

        ! not enough indices?
        nidx = ielsum(ovl(1,ica1),ngastp)
c dbg
c        print *,'nidx, rank: ',nidx,rank
c        print *,'type,type2',type,type2
c        print *,'ovl',ovl(type,ica1)
c dbg
        if (nidx.lt.rank) then
          ok = .false.
          exit
        end if

        occ_cnt = 0
        if (type.eq.0.and.nidx.eq.rank) then
          occ_cnt = ovl
        else if(rank.eq.1.and.ovl(type,ica1).ge.1) then
          occ_cnt = 0
          occ_cnt(type,ica1) = 1
        else if(rank.eq.1.and.type2.gt.0.and.ovl(type2,ica1).eq.1) then
          occ_cnt = 0
          occ_cnt(type2,ica1) = 1
        else if (rank.eq.2) then
c dbg
c          print *,'nblk: ',nblk(type)
c dbg
          do iblk = 1, nblk(type)
            ok = .true.
            do ihpvx = 1, ngastp
              if (occ_prj(ihpvx,ioff(type)+iblk).eq.0) cycle
              ok = ok.and.ovl(ihpvx,ica1).eq.
     &                    occ_prj(ihpvx,ioff(type)+iblk)
            end do
c dbg
c            print *,'block: ',iblk
c            print *,'projector: ',occ_prj(1:ngastp,ioff(type)+iblk)
c            print *,'ovl      : ',ovl(1:ngastp,ica1)
c            print *,'ok = ',ok
c dbg
            occ_cnt = 0
            occ_cnt(1:ngastp,ica1) = occ_prj(1:ngastp,ioff(type)+iblk)
            if (ok) exit
          end do
          if (.not.ok) exit
        else if (rank.gt.1) then
          ok = .false.
c presently:
          call quit(1,'set_inproj','I am spoilt for choice ...')
        else
          ok = .false.
          exit
        end if

        if (ntest.ge.100) then
          write(lulog,*) 'present ok: ',ok
          if (ok) then
            write(lulog,*)
     &       'I will establish the following connection of ',ivtx1,ivtx2
            call wrt_occ(lulog,occ_cnt)
          end if
        end if
        ! it is assumed that proto has enough space for us
        proto%narc = proto%narc+1
        idx = proto%narc
        proto%arc(idx)%link(1) = ivtx1
        proto%arc(idx)%link(2) = ivtx2
        proto%arc(idx)%occ_cnt = occ_cnt
      end do

      if (ok.and.proto%narc.ne.narc_old+ninproj)
     &     call quit(1,'set_inproj','consistency check failed')

      if (ok) then
        ! test whether the contraction might work
        occ2 = 0
        do ivtx = 1, proto%nvtx
          occ1 = occ_vtx(1:ngastp,1:2,ivtx)
          ! remove special contractions
          do idx = narc_old+1, narc_old+ninproj
            ivtx1 = proto%arc(idx)%link(1)
            ivtx2 = proto%arc(idx)%link(2)
            if (ivtx1.eq.ivtx) occ1 = occ1 - proto%arc(idx)%occ_cnt
            if (ivtx2.eq.ivtx)
     &             occ1 = occ1 - iocc_dagger(proto%arc(idx)%occ_cnt)
          end do
          occ2 = occ2
     &         + iocc_xdn(1,occ1)
     &         - iocc_dagger(iocc_xdn(2,occ1))
        end do
        ok = iocc_zero(occ2)
      end if

      if (ntest.ge.100)
     &     write(lulog,*) 'final OK: ',ok

      return
      end
