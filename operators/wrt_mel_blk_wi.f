*----------------------------------------------------------------------*      
      subroutine wrt_mel_blk_wi(luout,buffer,
     &     mel,iblk,igam,idxms,idxdis,nel,str_info,orb_info)
*----------------------------------------------------------------------*      
*     write ME-list block to unit luout with complete index for
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

      ! for debugging it is sometimes convenient to have at most:
      integer, parameter ::
     &     maxlines = 50 !5 !75
      ! set to -1 if you want the full output

      integer, intent(in) ::
     &     luout, iblk, igam, idxms, idxdis, nel
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
     &     did, idxstr, iel, ioff, msa, ihpv,
     &     idx_occ, njoined, idxst, idxnd, iadd, ijoin
      integer ::
     &     msdst(ngastp,2,mel%op%njoined),
     &     igamdst(ngastp,2,mel%op%njoined),
     &     lexlscr(nel,3),idorb(nel), idspn(nel), idspc(nel),
     &     nelc(mel%op%njoined), nela(mel%op%njoined)
      character ::
     &     spnstr*(nel+1+mel%op%njoined),fmtstr*256

      type(operator), pointer ::
     &     op

      logical, external ::
     &     next_tupel_ca

      njoined = mel%op%njoined
      idx_occ = (iblk-1)*njoined+1
      op => mel%op

      do ijoin = 1, njoined
        nelc(ijoin) = sum(op%ihpvca_occ(1:ngastp,1,idx_occ-1+ijoin))
        nela(ijoin) = sum(op%ihpvca_occ(1:ngastp,2,idx_occ-1+ijoin))
      end do
c      if (op%off_op_gmox(iblk)%maxd.eq.1) then
c        ! we cannot be sure that the extended information on
c        ! gmox is set for this operator, so ...
c        ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
c        ! ... and we set the Ms and IRREP distributions
c        ! ourselves .... 
c        msdst(1:ngastp,1:2) = 0
c        igamdst(1:ngastp,1:2) = 0
c        msa = nela - (idxms-1)*2
c        do ihpv = 1, ngastp
c          if (op%ihpvca_occ(ihpv,1,iblk).gt.0) then
c            msdst(ihpv,1) = msa + op%mst
c            igamdst(ihpv,1) = multd2h(igam,op%gamt)
c          end if
c          if (op%ihpvca_occ(ihpv,2,iblk).gt.0) then
c            msdst(ihpv,2) = msa
c            igamdst(ihpv,2) = igam
c          end if
c        end do
c      else
        ioff = mel%off_op_gmox(iblk)%d_gam_ms(idxdis,igam,idxms)
        ! get Ms and IRREP distribution info from
        ! distribution ID
        did = mel%off_op_gmox(iblk)%did(idxdis,igam,idxms)
        call did2msgm(msdst,igamdst,did,
     &               op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
c      end if
        
      fmtstr(1:) = '(x,i5,2x,'
      idxst = 10
      do ijoin = 1, njoined
        if (nelc(ijoin).gt.0) then
          write(fmtstr(idxst:),'(i2,"i3,x,")')
     &         nelc(ijoin)
          idxst = len_trim(fmtstr)+1
        end if
        if (nela(ijoin).gt.0) then
          write(fmtstr(idxst:),'(i2,"i3,x,")')
     &         nela(ijoin)
          idxst = len_trim(fmtstr)+1
        end if
      end do
      fmtstr(idxst:) = 'x,a,2x,'
      idxst = len_trim(fmtstr)+1
      do ijoin = 1, njoined
        if (nelc(ijoin).gt.0) then
          write(fmtstr(idxst:),'(i2,"i1,x,")')
     &         nelc(ijoin)
          idxst = len_trim(fmtstr)+1
        end if
        if (nela(ijoin).gt.0) then
          write(fmtstr(idxst:),'(i2,"i1,x,")')
     &         nela(ijoin)
          idxst = len_trim(fmtstr)+1
        end if
      end do
      fmtstr(idxst:) = 'g20.10)'
c dbg
c      print *,'fmtstr:"',trim(fmtstr),'"'
c dbg


c      write(fmtstr,'("(x,i5,2x,",i2,"i3,x,",i2,"i3,2x,a,2x,",'//
c     &                       'i2,"i1,x,",i2,"i1,x,g20.10)")')
c     &     nelc,nela,nelc,nela
      
      spnstr(1:nel+1+njoined) = ' '  

      ! loop over all possible index tuples:
      first = .true.
      idxstr = 0
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
        if (maxlines.gt.0.and.idxstr.gt.maxlines) exit

        idxnd = 0
        iadd = 0
        do ijoin = 1, njoined
          idxst = idxnd+1
          idxnd = idxst+nelc(ijoin)-1
          do iel = idxst, idxnd
            if (idspn(iel).eq.1)
     &           write(spnstr(iel+iadd:iel+iadd),'("+")')
            if (idspn(iel).eq.-1)
     &           write(spnstr(iel+iadd:iel+iadd),'("-")')
            if (idspn(iel).eq.2)
     &           write(spnstr(iel+iadd:iel+iadd),'("2")')
          end do
          iadd = iadd+1
          idxst = idxnd+1
          idxnd = idxst+nela(ijoin)-1
          do iel = idxst, idxnd
            if (idspn(iel).eq.1)
     &           write(spnstr(iel+iadd:iel+iadd),'("+")')
            if (idspn(iel).eq.-1)
     &           write(spnstr(iel+iadd:iel+iadd),'("-")')
            if (idspn(iel).eq.2)
     &           write(spnstr(iel+iadd:iel+iadd),'("2")')
          end do
          iadd = iadd+1
        end do

        write(luout,fmtstr) idxstr+ioff,
     &       idorb(1:nel), spnstr, idspc(1:nel),
     &       buffer(idxstr)

      end do

      return
      end
