*----------------------------------------------------------------------*
      subroutine set_flipmap_kernel(strmap,
     &     iocc,irestr,idxms,igam,
     &     g_y4sg,g_yinf,
     &     g_yssg,g_wssg,
     &     g_off_dgm,g_ndis,
     &     mostnd_cur,nsym,ngas_cur,igamorb)
*----------------------------------------------------------------------*
*     core routine for setting up the (str)->(str_spin_flipped) mapping
*
*     andreas, july 2008
*      
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 000

      integer, intent(out) ::
     &     strmap(*)
      integer, intent(in) ::
     &     iocc, irestr(2,ngas_cur,2), idxms, igam,
     &     nsym, ngas_cur,
     &     mostnd_cur(2,nsym),igamorb(*),
     &     g_y4sg(*),g_yinf(*),
     &     g_yssg(*),g_wssg(*),
     &     g_off_dgm(*),g_ndis

      logical ::
     &     first
      integer ::
     &     idorb(iocc),idspn(iocc),idgam(iocc),idss(iocc),
     &     idspn_flipped(iocc),
     &     istr, nstr, idx, ms, idxmap

      integer, external ::
     &     idx4sg, std_spsign
      logical, external ::
     &     next_string

      ms = iocc-(idxms-1)*2

      idxmap = 0

      first = .true.
      ! loop over original strings
      do while(next_string(idorb,idspn,idss,
     &     iocc,ms,igam,first,
     &     irestr,
     &     mostnd_cur,igamorb,
     &     nsym,ngas_cur))

        first = .false.
        
        ! flip spin indices
        idspn_flipped(1:iocc) = idspn(1:iocc)
        do idx = 1, iocc
          if (idspn_flipped(idx).ne.2)
     &         idspn_flipped(idx) = -idspn_flipped(idx)
        end do

        do idx = 1, iocc
          idgam(idx) = igamorb(idorb(idx))
        end do

        idxmap = idxmap+1

        strmap(idxmap) = std_spsign(idspn        ,iocc)*
     &                   std_spsign(idspn_flipped,iocc)
c dbg
c        print *,'idorb:   ',idorb(1:iocc)
c        print *,'idspn:   ',idspn(1:iocc)
c        print *,'idss:    ',idss (1:iocc)
c        print *,'idspn_fl:',idspn_flipped(1:iocc)
c        print *,'sign: ',strmap(idxmap)
c dbg
        strmap(idxmap) = strmap(idxmap)*
     &        (idx4sg(iocc,idss,idorb,idspn_flipped,idgam,
     &                g_y4sg,g_yinf,
     &                g_yssg,g_wssg,
     &                g_off_dgm,g_ndis,
     &                mostnd_cur,
     &                iocc,nsym,ngas_cur)+1)
c dbg
c        print *,'idx = ',abs(strmap(idxmap))
c dbg

      end do

      if (ntest.ge.1000) then
        write(luout,*) 'set_flipmap_kernel: the result of my work is'
        nstr = idxmap
        do istr = 1, nstr
          write(luout,'(x,i5,"->",i5)') istr,strmap(istr)
        end do
      end if

      return
      end

      integer function std_spsign(idspn,nel)
*
*     set up standard sign for spin-string
*     i.e. the sign for bringing spin string into order aaaabbb
*
      integer, intent(in) ::
     &     nel, idspn(nel)

      integer ::
     &     idx, trp, nb, sig

      trp = 0  ! sum of transpositions
      nb = 0   ! number of betas encountered up to current index
      idx = 1  ! index counter
      do while(idx.le.nel)
        sig = idspn(idx) ! get present sigma
        ! alpha or paired index (interpreted as ab):
        ! increment transposition count with number of prev. betas
        if (sig.eq. 1.or.sig.eq.2) trp = trp + nb
        if (sig.eq. 1.or.sig.eq.2) idx = idx+1 ! inc counter
        ! beta or second part of paired index:
        ! increment number of betas
        if (sig.eq.-1.or.sig.eq.2) nb = nb+1
        if (sig.eq.-1.or.sig.eq.2) idx = idx+1 ! inc counter
      end do

      std_spsign = +1
      if (mod(trp,2).ne.0) std_spsign = -1

      return
      end
