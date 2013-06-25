      subroutine pseudoinv(mel_inp,mel_inv,nocc_cls,
     &     op_info,orb_info)
*----------------------------------------------------------------------*
*
*     Routine to invert simple matrices.
*     Allows matrix to be singular, will give pseudoinverse in this case
*     Currently restricted to symmetric matrices
*
*     matthias, may 2013
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_orbinf.h'
      include 'mdef_operator_info.h'
      include 'ifc_memman.h'
      include 'routes.h'

      integer, parameter ::
     &     ntest = 00

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(me_list), intent(in) ::
     &     mel_inp, mel_inv
      integer, intent(in) ::
     &     nocc_cls

      logical ::
     &     bufin, bufout
      integer ::
     &     ifree, nbuff, idxmsa, iocc_cls,
     &     msmax, msa, igama, idx, jdx, ngam,
     &     len_gam_ms, ioff, lwrk, info
      real(8) ::
     &     trace
      real(8), pointer ::
     &     buffer_in(:), buffer_out(:), scratch(:,:), eigen_val(:),
     &     scratch2(:,:), wrk(:), dum1,dum2

      type(filinf), pointer ::
     &     ffinp, ffinv
      type(operator), pointer ::
     &     op_inv, op_inp

      integer, external ::
     &     idx_oplist2

      ffinp => mel_inp%fhand
      ffinv => mel_inv%fhand
      op_inp => mel_inp%op
      op_inv => mel_inv%op

      ! Check whether files are buffered.
      bufin = .false.
      if(ffinp%buffered) bufin = .true.
      bufout = .false.
      if(ffinv%buffered) bufout = .true.

      ifree = mem_setmark('pseudoinv')

      ! Number of irreps in symmetry group.
      ngam = orb_info%nsym

      ! Allocations made to maximum block length to save time.
      if(.not.bufin)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          if(op_inp%formal_blk(iocc_cls))
     &         cycle

          nbuff = nbuff+mel_inp%len_op_occ(iocc_cls)
        enddo
        ifree = mem_alloc_real(buffer_in,nbuff,'buffer_in')
        call get_vec(ffinp,buffer_in,1,nbuff)
      else
        if(ntest.ge.100)
     &       write(luout,*)'Invert: input not incore'
        buffer_in => ffinp%buffer(1:)
      endif

      if(.not.bufout)then
        nbuff = 0
        do iocc_cls = 1, nocc_cls
          nbuff = nbuff + mel_inv%len_op_occ(iocc_cls)
        enddo
        ifree= mem_alloc_real(buffer_out,nbuff,'buffer_out')
        buffer_out(1:nbuff) = 0d0
      else
        if(ntest.ge.100)
     &       write(luout,*)'Symmetrise: output not incore'
        buffer_out => ffinv%buffer(1:)
      endif

      ! Loop over occupation class.
      iocc_loop: do iocc_cls = 1, nocc_cls
        if(op_inp%formal_blk(iocc_cls)) cycle iocc_loop
c dbg
c       print *,'iocc:',iocc_cls
c dbgend

        ! Loop over Ms of annihilator string.
        idxmsa = 0
        msmax = op_inp%ica_occ(1,iocc_cls)
        msa_loop : do msa = msmax, -msmax, -2

          idxmsa = idxmsa+1
c dbg
c       print *,'msa,idxmsa:',msa,idxmsa
c dbgend
      
          ! Loop over Irrep of annihilator string.
          igama_loop: do igama =1, ngam

            len_gam_ms = int(sqrt(dble(mel_inv%
     &         len_op_gmo(iocc_cls)%gam_ms(igama,idxmsa))))
            if (len_gam_ms.eq.0) cycle

            ioff = mel_inv%off_op_gmo(iocc_cls)%gam_ms(igama,idxmsa)
c dbg
c     print *,'igama,ioff, len_gam_ms:',igama,ioff,len_gam_ms
c dbgend

            lwrk = max(1024,len_gam_ms**2)
            allocate(scratch(len_gam_ms,len_gam_ms),
     &               scratch2(len_gam_ms,len_gam_ms),
     &               eigen_val(len_gam_ms),
     &               wrk(lwrk))
            trace = 0d0
            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms
c dbg
c     print *,'idx,jdx,el:',idx,jdx,ioff+(idx-1)*len_gam_ms+jdx
c     print *,'value:',buffer_in(ioff+(idx-1)*len_gam_ms+jdx)
c dbgend
                scratch(idx,jdx) =
     &               buffer_in(ioff+(idx-1)*len_gam_ms+jdx)
              enddo
              trace = trace + scratch(idx,idx)
            enddo

            ! if trace is negative: change sign
            if (trace.lt.0d0) scratch(1:len_gam_ms,1:len_gam_ms)
     &            = -scratch(1:len_gam_ms,1:len_gam_ms)

c dbg
c            write(luout,*) 'matrix:'
c            call wrtmat2(scratch,len_gam_ms,len_gam_ms,
c     &                   len_gam_ms,len_gam_ms)
c cbgend

            ! check if matrix is symmetric
            ! routine could be generalized with SVD in general sense
            do idx = 2, len_gam_ms
              do jdx = 1,idx-1
                if(abs(scratch(idx,jdx)-scratch(jdx,idx)).gt.1d-10)then
                  write(luout,*) 'idx,jdx:',idx,jdx
                  write(luout,*) 'scratch(idx,jdx): ',scratch(idx,jdx)
                  write(luout,*) 'scratch(jdx,idx): ',scratch(jdx,idx)
                  call quit(1,'pseudoinv',
     &                 'currently only for symmetric matrix!')
                end if
              end do
            end do

            ! Perform singular value decomposition, ...
            info = 0
            call dgesvd('O','N',len_gam_ms,len_gam_ms,
     &                  scratch,len_gam_ms,eigen_val,
     &                  dum1,1,dum2,1,wrk,lwrk,info)
            if (info.ne.0) then
              write(luout,*) 'WARNING in invsqrt_mat: SVD in trouble'
              call quit(1,'pseudoinv','cannot handle this matrix')
            end if
            ! ... compute (pseudo)eigenvalues to the -1/2, ...
            do idx = 1, len_gam_ms
              if (abs(eigen_val(idx)).le.sv_thresh) then
                eigen_val(idx) = 0d0
              else
                eigen_val(idx) = 1d0/sqrt(eigen_val(idx))
              end if
              ! ..., multiply them with the columns of U, ...
              scratch2(1:len_gam_ms,idx)
     &                = eigen_val(idx)*scratch(1:len_gam_ms,idx)
            end do
            ! ... and compute U*inv*U^+
            call dgemm('n','t',len_gam_ms,len_gam_ms,len_gam_ms,
     &                 1d0,scratch2,len_gam_ms,
     &                     scratch2,len_gam_ms,
     &                 0d0,scratch,len_gam_ms)

            ! if trace was negativ: change back sign
            if (trace.lt.0d0) scratch(1:len_gam_ms,1:len_gam_ms)
     &            = -scratch(1:len_gam_ms,1:len_gam_ms)

            do idx = 1,len_gam_ms
              do jdx = 1,len_gam_ms
                buffer_out((idx-1)*len_gam_ms+jdx+ioff) =
     &               scratch(idx,jdx)                
              enddo
            enddo

            deallocate(scratch,scratch2,eigen_val,wrk)
                    
          enddo igama_loop
        
        enddo msa_loop

      enddo iocc_loop


      if(.not.bufout)then
        call put_vec(ffinv,buffer_out,1,nbuff)
      endif  

      ifree = mem_flushmark('pseudoinv')

      return
      end
