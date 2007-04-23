*----------------------------------------------------------------------*
      subroutine wrt_op_buf(luout,level,bufop,op,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for wrt_op: operator elements on incore buffer
*     see below for further info
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, level, iblkst, iblknd
      real(8), intent(in), target ::
     &     bufop(*)
      type(operator), intent(in) ::
     &     op
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      type(filinf) ::
     &     ffdum

      call wrt_op(luout,level,.true.,bufop,ffdum,op,iblkst,iblknd,
     &     str_info,orb_info)

      return
      end

*----------------------------------------------------------------------*
      subroutine wrt_op_file(luout,level,ffop,op,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for wrt_op: operator elements on file
*     see below for further info
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, level, iblkst, iblknd
      type(filinf), intent(in) ::
     &     ffop
      type(operator), intent(in) ::
     &     op
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      real(8) ::
     &     xdum

      call wrt_op(luout,level,.false.,xdum,ffop,op,iblkst,iblknd,
     &     str_info,orb_info)

      return
      end

*----------------------------------------------------------------------*
      subroutine wrt_op(luout,level,incore,bufop,ffop,op,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     given: an operator definition (op) and a file handle
*         or buffer containing the matrix elements
*     print out operator blocks iblkst to iblknd 
*     (details acc. to level:
*       0 - only total norm
*       1 - +norm of blocks
*       2 - +norm of  Ms/Gamma distributions
*       3 - +elements (up to 100)
*       4 - all elements 
*       5 - all elements + indices )
*     simple version -- load complete symmetry block into memory
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, level, iblkst, iblknd
      logical, intent(in) ::
     &     incore
      real(8), intent(in), target ::
     &     bufop(*)
      type(filinf), intent(in) ::
     &     ffop
      type(operator), intent(in) ::
     &     op
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      logical ::
     &     first, close_again, blk_buf, scalar
      integer ::
     &     idxoff, idxoff_blk, iblk, lenblk, lenprt, ifree, mmax,
     &     msmax, idxms, ms, igam, idx_dis, ndis, nwarn, did, idum, nel
      integer ::
     &     msd(ngastp,2), igamd(ngastp,2)
      real(8) ::
     &     xnrm, xnrm_tot
      real(8), pointer ::
     &     buffer(:), curblk(:)

      logical, external ::
     &     next_msgamdist
      real(8), external ::
     &     ddot

      ifree = mem_setmark('wrt_op')

      if (.not.incore) then
        
        mmax = 0
        do iblk = iblkst, iblknd
          msmax = op%ica_occ(2,iblk)
          idxms = 0
          do ms = msmax, -msmax, -2
            idxms = idxms+1
            do igam = 1, orb_info%nsym
              mmax = max(mmax,op%len_op_gmo(iblk)%gam_ms(igam,idxms))
            end do
          end do
        end do

        ifree = mem_alloc_real(buffer,mmax,'buffer')

        close_again = .false.
        if (ffop%unit.le.0) then
          close_again = .true.
          call file_open(ffop)
        end if
      end if

      xnrm_tot = 0d0
      do iblk = iblkst, iblknd
        
        if (.not.incore) then
          blk_buf = ffop%buffered
          if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0
        end if

        if (level.ge.1) then
          if (level.ge.2) write(luout,'("+",77("="),"+")')
          write(luout,'(2x,a,i4,a,i12)') 'block no. ',iblk,' len = ',
     &         op%len_op_occ(iblk)
          if (level.ge.2) call wrt_occ(luout,op%ihpvca_occ(1,1,iblk))
          if (level.ge.2) write(luout,'("+",77("="),"+")')
        end if

        scalar = min(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        nwarn = 0
        msmax = op%ica_occ(2,iblk)
        nel = msmax+op%ica_occ(1,iblk)
        idxms = 0
        idxoff = 0
        idxoff_blk = 0
        do ms = msmax, -msmax, -2
          idxms = idxms+1
          do igam = 1, orb_info%nsym

            ! block offest and length:
            idxoff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lenblk = op%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle

            ! get current block
            if (.not.incore.and..not.blk_buf) then
              call get_vec(ffop,buffer,idxoff+1,idxoff+lenblk)
              curblk => buffer(1:lenblk)
            else if (.not.incore) then
c              ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
              ! currently: idxoff should be valid here, as well
              curblk => ffop%buffer(idxoff+1:idxoff+lenblk)
            else
              curblk => bufop(idxoff+1:idxoff+lenblk)
            end if

            ! reset offset within current buffer
            idxoff_blk = 0
        
            xnrm = sqrt(ddot(lenblk,curblk,1,curblk,1))
            xnrm_tot = xnrm_tot + xnrm*xnrm

            if (level.gt.0) then
              lenprt = lenblk
              if (level.lt.4) lenprt = min(maxprt,lenblk)
              write(luout,'(2x,a,i3,a,i2,a,i12,a,g12.6)')
     &           'Ms(A) = ',ms,'/2  IRREP(A) = ',igam,'  len = ',lenblk,
     &           '  norm = ',xnrm
              if (level.ge.2) write(luout,'("+",77("-"),"+")')
            end if
            if (level.ge.2) then
              ! print out distributions
              first = .true.
              ndis = op%off_op_gmox(iblk)%ndis(igam,idxms)
              if (ndis.eq.0) then
                write(luout,*) 'WARNING:'
                write(luout,*)
     &               ' op_info indicates no distribution at all'
                write(luout,*) ' skipping to next block'
                idxoff_blk = idxoff_blk+lenblk
                nwarn = nwarn+1
                cycle
              end if
              if (ndis.eq.1) then
                write(luout,*) 'block contains only single distribution'
                if (level.ge.3) then
                  if (level.ge.5)
     &                 write(luout,*) 'index of first element:',idxoff+1
                  write(luout,'("+",77("."),"+")')
                  if (level.ge.5) then
                    if (scalar) then
                      write(luout,'(2x,6f12.7)')
     &                   curblk(idxoff_blk+1)
                    else
                      call wrt_op_blk_wi(luout,curblk(idxoff_blk+1),
     &                   op,iblk,igam,idxms,1,
     &                   nel,str_info,orb_info)
                    end if
                  else
                    write(luout,'(2x,6f12.7)')
     &                   curblk(idxoff_blk+1:idxoff_blk+lenprt)
                  end if
                end if
                write(luout,'("+",77("-"),"+")')
                idxoff_blk = idxoff_blk+lenblk
                cycle
              end if
              distr_loop: do idx_dis = 1, ndis
                  
                did = op%off_op_gmox(iblk)%did(idx_dis,igam,idxms)

                call did2msgm(msd,igamd,did,
     &               op%ihpvca_occ(1,1,iblk),orb_info%nsym)

                lenblk =
     &               op%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
                if (lenblk.eq.0) cycle
                xnrm =
     &               sqrt(ddot(lenblk,curblk(idxoff_blk+1),1,
     &                                curblk(idxoff_blk+1),1))
                lenprt = lenblk
                if (level.lt.4) lenprt = min(maxprt/ndis,lenblk)
              
                write(luout,'(2x,a,i12,a,g12.6)'),
     &             ' Ms-Dst     Gamma-Dst   len = ',lenblk,
     &             '  norm = ',xnrm
                if (ngastp.eq.2) then
                  write(luout,'(2x,"/",2i3,"\/",2i3,"\")')
     &                 msd(1:ngastp,1), igamd(1:ngastp,1)
                  write(luout,'(2x,"\",2i3,"/\",2i3,"/")')
     &                 msd(1:ngastp,2), igamd(1:ngastp,2)
                else if (ngastp.eq.3) then
                  write(luout,'(2x,"/",3i3,"\/",3i3,"\")')
     &                 msd(1:ngastp,1), igamd(1:ngastp,1)
                  write(luout,'(2x,"\",3i3,"/\",3i3,"/")')
     &                 msd(1:ngastp,2), igamd(1:ngastp,2)
                else if (ngastp.eq.4) then
                  write(luout,'(2x,"/",4i3,"\/",4i3,"\")')
     &                 msd(1:ngastp,1), igamd(1:ngastp,1)
                  write(luout,'(2x,"\",4i3,"/\",4i3,"/")')
     &                 msd(1:ngastp,2), igamd(1:ngastp,2)
                end if
                if (level.ge.3) then
                  if (level.ge.5)
     &                 write(luout,*) 'index of first element:',idxoff+1
                  write(luout,'("+",77("."),"+")')
                  if (level.ge.5) then
                    call wrt_op_blk_wi(luout,curblk(idxoff_blk+1),
     &                   op,iblk,igam,idxms,idx_dis,
     &                   nel,str_info,orb_info)
                  else
                    write(luout,'(2x,6f12.7)')
     &                 curblk(idxoff_blk+1:idxoff_blk+lenprt)
                  end if
                end if
                idxoff_blk = idxoff_blk+lenblk
                if (idx_dis.ne.ndis)
     &               write(luout,'("+",77("."),"+")')
              end do distr_loop
            end if
            if (level.ge.2) write(luout,'("+",77("-"),"+")')
          end do
        end do
      end do

      write(luout,*) 'total norm = ',sqrt(xnrm_tot)

      if (nwarn.gt.0) then
        write(luout,*) '!!! There were ',nwarn,' warnings !!!'
        write(luout,*) 'look for "WARNING" in previous output!'
      end if

      if (.not.incore.and.close_again) call file_close_keep(ffop)

      ifree = mem_flushmark()

      return
      end
