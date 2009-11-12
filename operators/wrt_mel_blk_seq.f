*----------------------------------------------------------------------*      
      subroutine wrt_mel_blk_seq(luwrt,buffer,
     &     mel,iblk,igam,idxms,idxdis,nel,nindex,maxlen,
     &     str_info,orb_info)
*----------------------------------------------------------------------*      
*     write ME-list block to unit luwrt with complete index for
*     each element
*----------------------------------------------------------------------*      

      implicit none

      include 'opdim.h'
      include 'hpvxseq.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'multd2h.h'

      integer, intent(in) ::
     &     luwrt, iblk, igam, idxms, idxdis, nel, nindex, maxlen
      type(me_list), intent(in) ::
     &     mel
      real(8), intent(in) ::
     &     buffer(*)
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      logical ::
     &     first
      integer ::
     &     did, idxstr, idxbuf, iel, ioff, msa, ihpv,
     &     idx_occ, njoined, idxst, idxnd, iadd, ijoin, idum
      integer ::
     &     msdst(ngastp,2,mel%op%njoined),
     &     igamdst(ngastp,2,mel%op%njoined),
     &     lexlscr(nel,3),idorb(nel), idspn(nel), idspc(nel),
     &     nelc(mel%op%njoined), nela(mel%op%njoined)

      integer(8) ::
     &     lenbuf, nel8
      integer(2) ::
     &     spins(nel,maxlen), indices(nel,maxlen)
      real(8) ::
     &     val(maxlen)

      type(operator), pointer ::
     &     op

      logical, external ::
     &     next_tupel_ca

      njoined = mel%op%njoined
      idx_occ = (iblk-1)*njoined+1
      op => mel%op

c          if(nindex.ne.nel) call quit(1,'wrt_mel_blk_seq',
c     &                                'wrong no of spins')

      do ijoin = 1, njoined
        nelc(ijoin) = sum(op%ihpvca_occ(1:ngastp,1,idx_occ-1+ijoin))
        nela(ijoin) = sum(op%ihpvca_occ(1:ngastp,2,idx_occ-1+ijoin))
      end do
      ioff = mel%off_op_gmox(iblk)%d_gam_ms(idxdis,igam,idxms)
      ! get Ms and IRREP distribution info from
      ! distribution ID
      did = mel%off_op_gmox(iblk)%did(idxdis,igam,idxms)
      call did2msgm(msdst,igamdst,did,
     &     op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
 
      ! loop over all possible index tuples:
      first = .true.
      idxstr = 0
      idxbuf = 0
      do while(next_tupel_ca(idorb,idspn,idspc,
     &     nel,njoined,op%ihpvca_occ(1,1,idx_occ),
     &     mel%idx_graph(1,1,idx_occ),
     &     msdst,igamdst,first,
     &     str_info%igas_restr,
     &     orb_info%mostnd,orb_info%igamorb,
     &     orb_info%nsym,orb_info%ngas,
     &     orb_info%ngas_hpv,orb_info%idx_gas,
     &     hpvxseq,lexlscr))
        first = .false.
        idxstr = idxstr+1
        idxbuf = idxbuf+1

        idxnd = 0
        iadd = 0
        do ijoin = 1, njoined
          idxst = idxnd+1
          idxnd = idxst+nelc(ijoin)+nela(ijoin)-1

          spins(idxst:idxnd,idxbuf) = idspn(idxst:idxnd)

c          do iel = idxst, idxnd
c            if (idspn(iel).eq.1)
c     &           write(spnstr(iel+iadd:iel+iadd),'("+")')
c            if (idspn(iel).eq.-1)
c     &           write(spnstr(iel+iadd:iel+iadd),'("-")')
c            if (idspn(iel).eq.2)
c     &           write(spnstr(iel+iadd:iel+iadd),'("2")')
c          end do
c          iadd = iadd+1
c          idxst = idxnd+1
c          idxnd = idxst+nela(ijoin)-1
c          do iel = idxst, idxnd
c            if (idspn(iel).eq.1)
c     &           write(spnstr(iel+iadd:iel+iadd),'("+")')
c            if (idspn(iel).eq.-1)
c     &           write(spnstr(iel+iadd:iel+iadd),'("-")')
c            if (idspn(iel).eq.2)
c     &           write(spnstr(iel+iadd:iel+iadd),'("2")')
c          end do
c          iadd = iadd+1

        end do

        indices(1:nel,idxbuf) = idorb(1:nel)
        val(idxbuf) = buffer(idxstr)

        if(idxbuf.lt.maxlen)cycle

        lenbuf = idxbuf
c        nel8   = nel
        nel8   = nindex
        write(luwrt) lenbuf,
c        write(luwrt) nel8,lenbuf,
     &       indices(1:nel8,1:lenbuf), spins(1:nel8,1:lenbuf),
     &       val(1:lenbuf)
c dbg
        print *,'nel8 = ',nel8
c dbg
        idxbuf = 0

      end do

      if(idxbuf.gt.0)then
        lenbuf = idxbuf
        nel8 = nindex
c        nel8 = nel
c dbg
c        print *,'nel8* = ',nel8
c dbg
        write(luwrt) lenbuf,
c        write(luwrt) nel8,lenbuf,
     &       indices(1:nel8,1:lenbuf), spins(1:nel8,1:lenbuf),
     &       val(1:lenbuf)
      endif

c dbg
c      do idum = 1,idxstr
c        write(6,*)indices(1:nindex,idum),spins(1:nindex,idum),
c     &       val(idum)
c      enddo
c dbg

      return
      end
