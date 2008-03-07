*----------------------------------------------------------------------*
      subroutine optc_prc_special(me_grd,me_special,nspecial,
     &                            nincore,xbuf1,xbuf2,xbuf3,lenbuf,
     &                            orb_info,str_info)
*----------------------------------------------------------------------*
*     experimental routine for testing special preconditioners
*     if nincore>1: xbuf1 contains the gradient vector on entry
*                   and should contain the preconditioned gradient
*                   on exit
*----------------------------------------------------------------------*
      implicit none

      include 'opdim.h'
      include 'stdunit.h'
      include 'mdef_me_list.h'
      include 'def_orbinf.h'
      include 'hpvxseq.h'
      include 'explicit.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'multd2h.h'
      include 'cc_routes.h'
      include 'ifc_input.h'


      integer, intent(in) ::
     &     nincore, lenbuf, nspecial
      type(me_list) ::
     &     me_grd
      type(me_list_array) ::
     &     me_special(nspecial)
      real(8), intent(inout) ::
     &     xbuf1(lenbuf), xbuf2(lenbuf), xbuf3(lenbuf)
      type(orbinf), intent(in), target ::
     &     orb_info
      type(strinf),intent(in), target ::
     &     str_info

      integer ::
     &     ngrd, nblk_grd, nbmat, nxmat, nfmat, nblk, idx, iblk, len,
     &     idxmsa, msa, msc, igama, igamc,  ngam, mscmax, idis,
     &     iocc_cls, nloop, istr, ica, ihpvdx, ihpv, nidx, ms_str,
     &     igamstr, ncidx, naidx, msamax, maxidx, igrph, idim,
     &     maxstrbuf, nstrbuf, ioffbuf0, nouter, n_inner2, n_inner1,
     &     lenblk, iouter, ioffbuf, ioff_inner1, ioff_inner2, idxbuf,
     &     idx_inner1, len_buf, bnblk, xnblk, ioff, ilen, idxms, jdx,
     &     kdx, nrot,
     &     max_rank, len_trip, off_trip
      logical ::
     &     first, first_str
      real(8) ::
     &     fac, xsum_outer, xsum_i1
      type(me_list), pointer ::
     &     me_bmat, me_xmat, me_fmat
      type(operator), pointer ::
     &     f_op

      integer ::
     &     f_hpvx(ngastp,2), msdst(ngastp,2), igamdst(ngastp,2),
     &     ioff_xsum(2*ngastp), nstr(2*ngastp), nincr(2*ngastp)

      integer, pointer ::
     &     len_blk(:), ld_blk(:), idxorb(:), idxspn(:), idxdss(:),
     &     igamorb(:), igas_restr(:,:,:,:), mostnd(:,:,:),
     &     idx_graph(:,:), idx_gas(:), ngas_hpv(:), b_off(:), x_off(:)
      real(8), pointer ::
     &     f_dia(:), xsum(:), buffer(:), b_elem(:), x_elem(:), bsum(:),
     &     xmet(:), temp(:), eigen_vec(:,:), eigen_val(:),
     &     xbuf_temp(:)

      integer, external ::
     &     ndisblk_mel, ndisblk_mel2, iblk_occ
      logical,external ::
     &     next_string, next_msgamdist

c dbg
      print *,'in optc_prc_special'
c      do idx =1,nspecial
c        print *,'labels ', me_special(idx)%mel%label
c      enddo
c dbg
      if (nincore.ne.3)
     &     call quit(1,'optc_prc_special',
     &     'currently: special preconditioning for nincore==3, only')

      if (nspecial.lt.1)
     &     call quit(1,'optc_prc_special',
     &     'nspecial must be >= 1')

      ngrd     = me_grd%len_op
      nblk_grd = ndisblk_mel(me_grd)

c dbg
      print *,'ngrd',ngrd
      print *,'nblk_grd',nblk_grd
c dbg

      if(trim(r12_apprx).eq.'A')then

        ! R12 code: here B^-1 was passed as me_special(1)%mel
        me_bmat => me_special(1)%mel
        nbmat = me_bmat%len_op
        ! should be open
        call vec_from_da(me_bmat%fhand,1,xbuf2,nbmat)

        ! Here I use two routines to get the subblock structure
        ! so that we can loop over these without caring too much
        ! about the details
        ! one might instead loop over MS(A), GAMMA(A) [, DISTRIBUTIONS]

        ! What level of R12 are we doing?
        call get_argument_value('method.R12','maxexc',ival=max_rank)

        ! total number of sub blocks in block 1
        nblk = ndisblk_mel2(me_bmat,1,1)

        ! here: the operator lengths must match
        !  will be different for R12-triples etc.
        if (ngrd.ne.nbmat.and.nblk.ne.nblk_grd.and.max_rank.le.2)
     &       call quit(1,'optc_prc_special',
     &       'block structures do not match')

        ! lengths and leading dimensions of sub-blocks
        allocate(len_blk(nblk),ld_blk(nblk))
        call set_disblkdim_mel(len_blk,ld_blk,me_bmat)

        idx = 1
        do iblk = 1, nblk
          len = ld_blk(iblk)
          call dgemm('n','n',len,len,len,
     &         1d0,xbuf2(idx),len,
     &         xbuf1(idx),len,
     &         0d0,xbuf3(idx),len)
          idx = idx + len_blk(iblk)
        end do

        deallocate(len_blk,ld_blk)

        if(do_cc.and.max_rank.gt.2)then
          ! Precondition the triples part of the R12 equations.
          len_trip = me_grd%len_op_occ(2)
          off_trip = me_grd%off_op_occ(2)

c dbg
          print *,'len,off',len_trip,off_trip
c dbg

          allocate(xbuf_temp(len_trip))
          xbuf_temp(1:len_trip)
     &         = xbuf1(off_trip+1:off_trip+len_trip)
          





          deallocate(xbuf_temp)
          print *,'Here I go again'
          stop
        endif

      else !if(trim(r12_apprx).eq."A'")then
        ! Point to information arrays.
        igas_restr => str_info%igas_restr
        mostnd => orb_info%mostnd
        igamorb => orb_info%igamorb
        idx_gas => orb_info%idx_gas
        ngas_hpv => orb_info%ngas_hpv
        ngam = orb_info%nsym
          
        ! Point to actual operators and their me-lists.
        me_bmat => me_special(1)%mel
        nbmat = me_bmat%len_op

        me_xmat => me_special(2)%mel
        nxmat = me_xmat%len_op

        me_fmat => me_special(3)%mel
        nfmat = me_fmat%len_op
        f_op => me_fmat%op

        ! Extract the B and X matrices.
c dbg
        ! Perhaps need better memory allocation.
c dbg
        allocate(b_elem((4*orb_info%ntoob**2)**2))
        bnblk = ndisblk_mel(me_bmat)*10
        allocate(b_off(bnblk))
        call two_from_op(b_elem,b_off,bnblk,me_bmat,orb_info)

        allocate(x_elem((4*orb_info%ntoob**2)**2))
        xnblk = ndisblk_mel(me_xmat)*10
        allocate(x_off(xnblk))
        call two_from_op(x_elem,x_off,xnblk,me_xmat,orb_info)


        ! Find necessary block of the Fock operator, f_ij.
        f_hpvx(1,1:2)=2
        f_hpvx(2:ngastp,1:2)=0
        iocc_cls=iblk_occ(f_hpvx,.false.,f_op)
        idx_graph => me_fmat%idx_graph(1:ngastp,1:2,iocc_cls)

        ncidx=f_op%ica_occ(1,iocc_cls)
        naidx=f_op%ica_occ(2,iocc_cls)
          
        ! Limits of spin.
        mscmax = ncidx
        msamax = naidx
        maxidx=max(ncidx,naidx)
        allocate(idxorb(maxidx),idxspn(maxidx),idxdss(maxidx))

        ! Allocate buffer space.
        len_buf=me_fmat%len_op_occ(iocc_cls)
        allocate(buffer(len_buf))
        buffer(1:len_buf)=0d0

        ! Retrieve the actual elements of the Fock matrix -> f_dia
        allocate(f_dia(2*orb_info%ntoob))
        call onedia_from_op(f_dia,me_fmat,orb_info)
        
        ! Loop with idim=1 gets the necessary dimensions of xsum.
        ! Loop with idim=2 does the actual work.
        maxstrbuf=0
        idim_loop:do idim=1,2

          ! Loop over Ms of annihilator string.
          idxmsa = 0
          ioffbuf0=0
          msa_loop: do msa = msamax, -msamax, -2
            msc = msa + me_fmat%mst

            if(abs(msc).gt.mscmax) cycle msa_loop
            idxmsa = idxmsa+1

            ! Loop over irrep of annihilator string.
            gam_loop: do igama = 1, ngam

              igamc=multd2h(igama,me_fmat%gamt)
              if (me_fmat%len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)
     &             .eq.0) cycle gam_loop
             
              idis = 0
              first =.true.
              ! Loop over strings.
              distr_loop: do
                if(.not.next_msgamdist(first,
     &               msa,msc,igama,igamc,f_hpvx,ngam,msdst,igamdst))
     &               exit distr_loop
              
                first=.false.
                idis=idis+1

                if(me_fmat%len_op_gmox(iocc_cls)%
     &               d_gam_ms(idis,igama,idxmsa).eq.0) cycle distr_loop


                nloop=0
                istr=0
                nstrbuf=0
                ihpv = 1
                  
                fac=1d0
                nloop=nloop+1
                ioff_xsum(nloop)=istr
                nidx=2
                igrph=idx_graph(ihpv,1)
                ms_str=msdst(ihpv,1)
                igamstr=igamdst(ihpv,1)
                  
                if(idim.eq.1)then
                  nstrbuf=nstrbuf+
     &                 str_info%g(igrph)%lenstr_gm(igamstr,idxmsa)
                else
                    
                  first_str=.true.
                  ! Loop over all possible strings in the operator
                  ! block.
                  str_loop: do
                    if(.not.next_string(idxorb,idxspn,idxdss,
     &                 nidx,ms_str,igamstr,first_str,
     &                 igas_restr(1,1,1,igrph),
     &                 mostnd(1,1,idx_gas(ihpv)),igamorb,
     &                 ngam,ngas_hpv(ihpv))) exit str_loop

                    first_str = .false.
                    istr = istr + 1
                    xsum(istr)=0d0
                    do idx = 1,nidx
                      xsum(istr) = xsum(istr)+fac*f_dia(idxorb(idx))
                    enddo
                  enddo str_loop
                  ! Extract the correct elements of B and X.
                  idxms = (2-ms_str)/2+1
                  idx = (idxms-1)*ngam + igamstr
                  ilen = str_info%g(igrph)%lenstr_gm(igamstr,idxms)
                  ioff = b_off(idx)
                  bsum(1:ilen**2) = b_elem(ioff+1:ioff+ilen**2)
                  ioff = x_off(idx)
                  xmet(1:ilen**2) = x_elem(ioff+1:ioff+ilen**2)

                  nstr(nloop)=istr-ioff_xsum(nloop)
                  ! Disregard empty strings.
                  if(nstr(nloop).eq.0) nloop=nloop-1
                endif

                if(idim.eq.1)then
                  maxstrbuf=max(maxstrbuf,nstrbuf)
                else
                  ! No contribution?
                  if(nloop.eq.0)cycle distr_loop

                  ! Only one loop.
                  if(nloop.eq.1)then
c                    buffer(ioffbuf0+1:ioffbuf0+nstr(1)) =
c     &                   xsum(1:nstr(1))
                    do idx = 1, nstr(1)
                      idxbuf = 0
                      xsum_outer = xsum(ioff_xsum(nloop)+idx)
                      do jdx = 1, nstr(1)**2
                        idxbuf = idxbuf+1
                        buffer(idxbuf) =
     &                       bsum(jdx)-xsum_outer*xmet(jdx)
                      enddo

c                      ! Check the eigenvalues.
c                      if(nstr(1).gt.0)then
c                        allocate(eigen_vec(nstr(1),nstr(1)),
c     &                       eigen_val(nstr(1)))
c                        call jacobi(buffer,nstr(1),nstr(1),eigen_val,
c     &                       eigen_vec,nrot)
c
c                        print *,'eigen',eigen_val(1:nstr(1))
c                        deallocate(eigen_val,eigen_vec)
c                      endif


                      ! Invert the buffer.
                      call gaussj(buffer,nstr(1),nstr(1))

                      ! Multiply by the correct gradient elements.
                      temp(1:maxstrbuf**2)=0d0
                      do jdx = 1,nstr(1)
                        do kdx = 1,nstr(1)
                          temp(jdx) = temp(jdx) +
     &                         buffer(nstr(1)*(jdx-1)+kdx)*
     &                         xbuf1(nstr(1)*(idx-1)+ioff+kdx)
                        enddo
                        xbuf3(ioff+nstr(1)*(idx-1)+jdx)=
     &                       temp(jdx)
                      enddo
                      
                    enddo
                    ioffbuf0 = ioffbuf0+idxbuf

                    cycle distr_loop
                  endif

                  ! More than one loop.
                  nouter = nloop-2
                  n_inner2 = nstr(nloop)
                  n_inner1 = nstr(nloop-1)
                  lenblk = n_inner1*n_inner2
                  do iouter = nouter, 1, -1
                    nincr(iouter) = lenblk
                    lenblk = lenblk*nstr(iouter)
                  enddo

                  ! assemble the actual diagonal:
                  ioffbuf = ioffbuf0
                  ioff_inner1 = ioff_xsum(nouter+1)
                  ioff_inner2 = ioff_xsum(nouter+2)
                  do while(ioffbuf.lt.ioffbuf0+lenblk)
                    ! collect contributions from formal outer loops 
                    ! (over more than two strings):
                    idxbuf = ioffbuf+1
                    xsum_outer = 0d0
                    do iouter = 1, nouter
                      idx = ioff_xsum(iouter) + idxbuf/nincr(iouter)+1
                      idxbuf = mod(idxbuf,nincr(iouter)) 
                      xsum_outer = xsum_outer + xsum(idx)
                    enddo
                    ! explicit loops for the two innermost strings
                    do idx_inner1 = ioff_inner1+1, ioff_inner1+n_inner1
                      xsum_i1 = xsum(idx_inner1) + xsum_outer
                      buffer(ioffbuf+1:ioffbuf+n_inner2) =
     &                     xsum(ioff_inner2+1:ioff_inner2+n_inner2) +
     &                     xsum_i1
                      ioffbuf = ioffbuf+n_inner2
                    enddo
                  enddo
                  ioffbuf0 = ioffbuf
        
                endif
              enddo distr_loop
            
            enddo gam_loop

            ! Allocation of space for xsum.
            if (msa.eq.-msamax.and.idim.eq.1)then
              allocate(xsum(maxstrbuf),bsum(maxstrbuf**2),
     &             xmet(maxstrbuf**2),temp(maxstrbuf**2))
              xsum(1:maxstrbuf)=0d0
              bsum(1:maxstrbuf)=0d0
              xmet(1:maxstrbuf)=0d0
            endif
          
          enddo msa_loop
        enddo idim_loop

        deallocate(temp,xmet,bsum)
        deallocate(xsum,f_dia,buffer)
        deallocate(idxdss,idxspn,idxorb)
        deallocate(x_off)
        deallocate(x_elem)
        deallocate(b_off)
        deallocate(b_elem)

      endif

      return
      end


