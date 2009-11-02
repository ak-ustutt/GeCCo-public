      subroutine add_unity(fac,mel_out,iblkout,orb_info,str_info)
*----------------------------------------------------------------------*
*
*     Add (antisymmetrised) unit operator
*
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_orbinf.h'
      include 'def_me_list.h'
      include 'ifc_memman.h'
      include 'hpvxseq.h'
      include 'multd2h.h'

      integer, parameter ::
     &     ntest = 00

      real(8), intent(in)::
     &     fac
      type(me_list), intent(inout) ::
     &     mel_out
      type(strinf), intent(in), target ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      integer, intent(in) ::
     &     iblkout

      logical ::
     &     bufout,closeit, first
      integer ::
     &     len_str, idum, ifree, lblk, nblkmax, nblk, nbuff,
     &     ioffin, ioffout, idxst, idxnd, njoined, join_off,
     &     idoffin, idoffout, idxmsa, msmax, msa, msc, igama, igamc,
     &     idx, ngam, len_gam_ms, ioff, len_blk, ioff_blk,
     &     igrph, idxms, igamstr, len(ngastp), idx1, idx2, idx3, idx4,
     &     ihpv, ica, iblk, len_dst, ioff_fac(ngastp), ihpvdx,
     &     ioff1, ioff2, ioff3, ioff4, len_c
      integer ::
     &     opout_temp(ngastp,2), msdst(ngastp,2), igamdst(ngastp,2),
     &     iocc(ngastp,2), idx_graph(ngastp,2)
      real(8) ::
     &     xsign
      
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:)

      type(filinf), pointer ::
     &     ffout
      type(operator), pointer ::
     &     opout

      logical, external ::
     &     iocc_equal, irestr_equal, occ_is_diag_blk,
     &     next_msgamdist

      ffout => mel_out%fhand
      opout => mel_out%op

      if (ntest.ge.100) then
        write(luout,*) '=========================='
        write(luout,*) ' add_unity messing around'
        write(luout,*) '=========================='
        write(luout,*) ' fac = ',fac
        write(luout,*) ' list = ',trim(mel_out%label)
        write(luout,*) ' ffout: ',trim(ffout%name),
     &                   ' rec: ',ffout%current_record
        write(luout,*) ' opout: ',opout%name(1:len_trim(opout%name)),
     &       '  block: ',iblkout
      end if

      closeit = .false.
      njoined = opout%njoined

      join_off=(iblkout-1)*njoined

      if (ngastp.ne.4) call quit(1,'add_unity','only for ngastp=4')
      iocc = 0
      idx_graph = 0
      do iblk = 1, njoined
        do ihpv = 1, ngastp
          do ica = 1, 2
            if (opout%ihpvca_occ(ihpv,ica,join_off+iblk).ne.0) then
              if (iocc(ihpv,ica).ne.0)
     &            call quit(1,'add_unity','cannot handle this operator')
              iocc(ihpv,ica) = opout%ihpvca_occ(ihpv,ica,join_off+iblk)
              idx_graph(ihpv,ica) = mel_out%idx_graph(ihpv,ica,
     &                               join_off+iblk)
            end if
          end do
        end do
      end do

      ! check whether the out operator is a diagonal block:
      if (.not.
     &     occ_is_diag_blk(iocc(1:ngastp,1:2),1))
     &     return ! else we just tacitly return

      bufout = .false.
      if(ffout%buffered) bufout = ffout%incore(iblkout).gt.0

      ifree = mem_setmark('add_unity')

      ngam = orb_info%nsym

      ! Find total block length.
      ioff_blk = mel_out%off_op_occ(iblkout)
      len_blk  = mel_out%len_op_occ(iblkout)

      if(len_blk.gt.ifree)
     &     call quit(1,'add_unity','insufficient space')
      
      ! Allocations made to maximum block length to save time.
      if(.not.bufout)then
        nbuff=len_blk
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
      else
        if(ntest.ge.100)
     &       write(luout,*)'Add_unity: Not incore'
        buffer_out => ffout%buffer(ioff_blk+1:)
      endif

      if(.not.bufout)then
        if(ffout%unit.le.0)then
          call file_open(ffout)
          closeit = .true.
        endif
        call get_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  

      ! due to normal ordering, hole lines introduce sign flips:
      xsign = dble(1-2*mod(iocc(1,1),2))

      ! Loop over Ms of annihilator string.
      idxmsa = 0
      ioff = 0
c      msmax = 2
      msmax = opout%ica_occ(2,iblkout)
      msa_loop : do msa = msmax, -msmax, -2

        idxmsa = idxmsa+1
        msc = msa + mel_out%mst
        ! Usually have mst=0 operators => Ms(c)=Ms(a)
      
        ! Loop over Irrep of annihilator string.
        igama_loop: do igama =1, ngam
          igamc = multd2h(igama,mel_out%gamt)


          first = .true.
          distr_loop0: do
            if (.not.next_msgamdist(first,
     &         msa,msc,igama,igamc,iocc,ngam,
     &         msdst,igamdst)) exit distr_loop0
            first = .false.

            len(1:ngastp) = 1
            ioff_fac(1:ngastp) = 1
            len_dst = 1
            do ihpvdx = 1, ngastp
              ihpv = hpvxseq(ihpvdx)
              if (iocc(ihpv,1).eq.0) cycle
              igrph = idx_graph(ihpv,1)
              igamstr = igamdst(ihpv,1)
              idxms = (iocc(ihpv,1)-msdst(ihpv,1))/2+1
              len(ihpv) = str_info%g(igrph)%lenstr_gm(igamstr,idxms)
              igrph = idx_graph(ihpv,2)
              igamstr = igamdst(ihpv,2)
              idxms = (iocc(ihpv,2)-msdst(ihpv,2))/2+1
              len_c = str_info%g(igrph)%lenstr_gm(igamstr,idxms)
              ioff_fac(ihpv) = len_dst
              len_dst = len_dst*len(ihpv)*len_c
            end do
            if (len_dst.eq.0) cycle distr_loop0
            idxst = ioff + 1
            ioff = ioff + len_dst

            ! ms-dst and gamma-dst should also be diagonal:
            if (.not.occ_is_diag_blk(msdst(1:ngastp,1:2),1).or.
     &          .not.occ_is_diag_blk(igamdst(1:ngastp,1:2),1))
     &           cycle distr_loop0

            if (ntest.ge.100) then
              write(luout,*) 'msc,msa,igamc,igama: ',
     &             msc,msa,igamc,igama
              write(luout,*) 'current distribution (MS,GAMMA):'
              call wrt_occ(luout,msdst)
              call wrt_occ(luout,igamdst)
            end if

            ! add constant to diagonal elements only
            do idx1 = 1, len(hpvxseq(4))
              ioff1 = (idx1-1)*(len(hpvxseq(4))+1)
     &                        *ioff_fac(hpvxseq(4))
              do idx2 = 1, len(hpvxseq(3))
                ioff2 = (idx2-1)*(len(hpvxseq(3))+1)
     &                          *ioff_fac(hpvxseq(3))
                do idx3 = 1, len(hpvxseq(2))
                  ioff3 = (idx3-1)*(len(hpvxseq(2))+1)
     &                            *ioff_fac(hpvxseq(2))
                  do idx4 = 1, len(hpvxseq(1))
                    ioff4 = (idx4-1)*(len(hpvxseq(1))+1)
     &                              *ioff_fac(hpvxseq(1))
                    buffer_out(ioff1+ioff2+ioff3+ioff4+idxst)=xsign*fac+
     &                buffer_out(ioff1+ioff2+ioff3+ioff4+idxst)
                  end do
                end do
              end do
            end do

          end do distr_loop0

c          len_gam_ms=int(sqrt(dble(mel_out%
c     &         len_op_gmo(iblkout)%gam_ms(igama,idxmsa))))
c
c          ioff=mel_out%off_op_gmo(iblkout)%gam_ms(igama,idxmsa)-
c     &         ioff_blk
c          idx_loop: do idx =1,len_gam_ms
c
c            buffer_out((idx-1)*len_gam_ms+idx+ioff) = 1d0*fac+
c     &         buffer_out((idx-1)*len_gam_ms+idx+ioff)
c
c          enddo idx_loop
          
        enddo igama_loop
          
      enddo msa_loop

      if(.not.bufout)then
        call put_vec(ffout,buffer_out,ioff_blk+1,ioff_blk+len_blk)
      endif  
      if(closeit)
     &     call file_close_keep(ffout)

      ifree = mem_flushmark('add_unity')

      return
      end
