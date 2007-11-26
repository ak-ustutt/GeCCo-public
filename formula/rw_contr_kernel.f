*----------------------------------------------------------------------*
      subroutine rw_contr_kernel(irw,lu,contr)
*----------------------------------------------------------------------*
*     read/write contraction from/to lu and expand/pack
*     long form on contr
*----------------------------------------------------------------------*

      implicit none

      integer, parameter ::
     &     ntest = 00

      include 'opdim.h'
      include 'stdunit.h'
      include 'def_contraction.h'

      integer, intent(in) ::
     &     lu, irw
      type(contraction), intent(inout) ::
     &     contr

      integer, parameter ::
     &     lbuf = 1024
      integer(2) ::
     &     buffer(lbuf)
      integer ::
     &     idx, ierr, ii, ica, igastp

      if (irw.gt.0) then

        read(lu,end=100) contr%fac,idx,buffer(1:idx)
        ! if we make it up to here in this case:
        if (idx.gt.lbuf)
     &     call quit(0,'rd_contr','too long record')

        contr%iblk_res = buffer(1)
        contr%nvtx = buffer(2)
        contr%nsupvtx = buffer(3)
        contr%narc = buffer(4) 
        contr%nfac = buffer(5) 

        ! (re)allocate the necessary space
        call resize_contr(contr,contr%nvtx,contr%narc,contr%nfac)
c        if (contr%nvtx.gt.contr%mxvtx) then
c          if (contr%mxvtx.gt.0) deallocate(contr%vertex)
c          contr%mxvtx = contr%nvtx
c          allocate(contr%vertex(contr%mxvtx))
c        end if
c        if (contr%narc.gt.contr%mxarc) then
c          if (contr%mxarc.gt.0) deallocate(contr%arc)
c          contr%mxarc = contr%narc
c          allocate(contr%arc(contr%mxarc))
c        end if
c        if (contr%nfac.gt.contr%mxfac) then
c          if (contr%mxfac.gt.0) deallocate(contr%inffac)
c          contr%mxfac = contr%nfac
c          allocate(contr%inffac(ld_inffac,contr%mxfac))
c        end if

        idx = 5

        do ii = 1, contr%nvtx
          contr%vertex(ii)%idx_op  = buffer(idx+1)
          contr%vertex(ii)%iblk_op = buffer(idx+2)
          idx = idx+2
        end do

        do ii = 1, contr%nvtx
          contr%svertex(ii) = buffer(idx+1)
          idx = idx+1
        end do

        ! set up "joined" array from this info
        call update_svtx4contr(contr)

        do ii = 1, contr%narc
          contr%arc(ii)%link(1)      = buffer(idx+1)
          contr%arc(ii)%link(2)      = buffer(idx+2)
          idx = idx+2
          do ica = 1,2
            do igastp = 1, ngastp
              idx = idx+1
              contr%arc(ii)%occ_cnt(igastp,ica) = buffer(idx)
            end do
          end do
        end do

        do ii = 1, contr%nfac
          contr%inffac(1,ii) = buffer(idx+1)
          contr%inffac(2,ii) = buffer(idx+2)
          contr%inffac(3,ii) = buffer(idx+3)
          contr%inffac(4,ii) = buffer(idx+4)
          contr%inffac(5,ii) = buffer(idx+5)
          idx = idx+5
        end do

        return

      else

        ! add contraction record
        ierr = 0
        buffer(1) = contr%iblk_res
        buffer(2) = contr%nvtx
        buffer(3) = contr%nsupvtx
        buffer(4) = contr%narc
        buffer(5) = contr%nfac
        if (buffer(1).ne.contr%iblk_res) ierr = ierr+1
        if (buffer(2).ne.contr%nvtx) ierr = ierr+1
        if (buffer(3).ne.contr%nsupvtx)ierr = ierr+1
        if (buffer(4).ne.contr%narc) ierr = ierr+1
        if (buffer(5).ne.contr%nfac) ierr = ierr+1
        idx = 5
        if (ierr.gt.0) goto 101
        if (idx+contr%nvtx*2.gt.lbuf) goto 103
        
        do ii = 1, contr%nvtx
          buffer(idx+1) = contr%vertex(ii)%idx_op
          if (buffer(idx+1).ne.contr%vertex(ii)%idx_op) ierr = ierr+1
          buffer(idx+2) = contr%vertex(ii)%iblk_op
          if (buffer(idx+2).ne.contr%vertex(ii)%iblk_op) ierr = ierr+1
          idx = idx+2
        end do
        if (ierr.gt.0) goto 102
        
        if (idx+contr%nvtx.gt.lbuf) goto 103

        do ii = 1, contr%nvtx
          buffer(idx+1) = contr%svertex(ii)
          idx = idx+1
        end do

        ! the "joined" array needs not be saved (see above)

        if (idx+contr%narc*8.gt.lbuf) goto 103
        do ii = 1, contr%narc
          buffer(idx+1) = contr%arc(ii)%link(1)
          buffer(idx+2) = contr%arc(ii)%link(2)
          idx = idx+2
          do ica = 1, 2
            do igastp = 1, ngastp
              idx = idx+1
              buffer(idx) = contr%arc(ii)%occ_cnt(igastp,ica)
            end do
          end do
        end do
        
        if (idx+contr%nfac*3.gt.lbuf) goto 103
        do ii = 1, contr%nfac
          buffer(idx+1) = contr%inffac(1,ii)
          buffer(idx+2) = contr%inffac(2,ii)
          buffer(idx+3) = contr%inffac(3,ii)
          buffer(idx+4) = contr%inffac(4,ii)
          buffer(idx+5) = contr%inffac(5,ii)
          idx = idx+5
        end do
        
        write(lu) contr%fac,idx,buffer(1:idx)
        
        if (ntest.ge.100) then
          write(luout,*) 'wrote ',8+4+idx,' bytes'
        end if
        
        return

      end if

      ! error handling
      ! unexpected EOF on read
 100  call quit(1,'rw_contr_kernel',
     &     'unexpected EOF while reading formula file')

      ! write errors:
 101  write(luout,*) 'too large nvtx, narc'
      goto 200
 102  write(luout,*)
     &     'too large numbers in vertex description'
      goto 200
 103  write(luout,*) 'insufficient buffer size'

 200  call quit(1,'wrt_contr','error writing contraction')

      end
