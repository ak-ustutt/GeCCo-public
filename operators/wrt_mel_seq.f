*----------------------------------------------------------------------*
      subroutine wrt_mel_seq(luout,mel,iblkst,iblknd,
     &     str_info,orb_info)
*----------------------------------------------------------------------*
*     wrapper for wrt_mel: matrix-elements on file
*     see below for further info
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout, iblkst, iblknd
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info

      real(8) ::
     &     xdum

      call wrt_mel2(luout,.false.,xdum,mel,iblkst,iblknd,
     &     str_info,orb_info)

      return
      end

*----------------------------------------------------------------------*
      subroutine wrt_mel2(luout,incore,bufmel,mel,iblkst,
     &     iblknd,str_info,orb_info)
*----------------------------------------------------------------------*
*     given: an operator definition (op) and a file handle
*         or buffer containing the matrix elements
*     copy the elements to a sequential file including integral value, 
*     indices and spins.
*     Initially used for printing Z and P-intermediates to enable frozen
*     core calculations.
*     Modified from wrt_mel.f GWR August 2008
*----------------------------------------------------------------------*
      implicit none
      include 'opdim.h'
      include 'def_filinf.h'
      include 'def_operator.h'
      include 'def_me_list.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'ifc_operators.h'
      include 'ifc_memman.h'
      include 'par_opnames_gen.h'

      integer, parameter ::
     &     maxprt = 100

      integer, intent(in) ::
     &     luout,  iblkst, iblknd
      logical, intent(in) ::
     &     incore
      real(8), intent(in), target ::
     &     bufmel(*)
      type(me_list), intent(in) ::
     &     mel
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in) ::
     &     orb_info
      
      logical ::
     &     first, close_again, blk_buf, scalar, dagger
      integer ::
     &     idoff, idxoff, idxoff_blk, iblk, lenblk, lenprt, ifree, mmax,
     &     msamax, mscmax, idxms, ms, igam, idx_dis, ndis, nwarn, did,
     &     idum, nel, mst,
     &     idxoff0, njoined, idx_occ, luwrt
      integer ::
     &     msd(ngastp,2,mel%op%njoined), igamd(ngastp,2,mel%op%njoined),
     &     scr(ngastp,2,2*mel%op%njoined), occ(ngastp,2,mel%op%njoined)
      real(8), pointer ::
     &     buffer(:), curblk(:)

      integer(8) ::
     &     nindex, maxlen

c dbg
      integer(2), pointer ::
     &     spins(:,:), indices(:,:)
      real(8),pointer ::
     &     val(:)
      integer ::
     &     idxstr
c dbg

      type(operator), pointer ::
     &     op
      type(filinf), pointer ::
     &     ffop
      type(filinf) ::
     &     fftemp
      integer, external ::
     &     iopen_nus
      logical, external ::
     &     next_msgamdist
      real(8), external ::
     &     ddot

      ifree = mem_setmark('wrt_op_seq')

      op => mel%op
      ffop => mel%fhand
      mst = mel%mst
      dagger = mel%op%dagger

      nwarn = 0

      ! Get the actual list (from file or memory)
      if (.not.incore) then
        
        mmax = 0
        do iblk = iblkst, iblknd
          if(op%formal_blk(iblk))cycle
          mscmax = op%ica_occ(1,iblk)
          msamax = op%ica_occ(2,iblk)
          idxms = 0
          do ms = msamax, -msamax, -2
            if (abs(mst+ms).gt.mscmax) cycle
            idxms = idxms+1
            do igam = 1, orb_info%nsym
              mmax = max(mmax,mel%len_op_gmo(iblk)%gam_ms(igam,idxms))
            end do
          end do
        end do
        ifree = mem_alloc_real(buffer,mmax,'buffer')

        close_again = .false.
        if (ffop%unit.le.0) then
          close_again = .true.
          call file_open(ffop)
        end if
      else
        ! incore: we assume that bufmel starts with first block
        !         to be displayed
        idxoff0 = mel%off_op_occ(iblkst)
      end if

      ! Open the file used to store the elements/indices.
c      fftemp%name = trim(op%name)
c      fftemp%type = 2
c      fftemp%unit = -1
      call file_init(fftemp,trim(op%name),ftyp_sq_unf,0)
      call file_open(fftemp)
      luwrt = fftemp%unit
      ! Write initial info on maximum block length and number of indices.
      maxlen = 1024
      nindex = 4
      if(trim(op%name).eq.op_z_inter) nindex = 6
      write(luwrt) nindex,maxlen

      njoined = mel%op%njoined

      idx_occ = (iblkst-1)*njoined+1
      do iblk = iblkst, iblknd
        if(op%formal_blk(iblk))cycle

        occ = op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        if (dagger) occ = iocc_dagger_n(occ,njoined)

        if (.not.incore) then
          blk_buf = ffop%buffered
          if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0
        end if

        scalar = max(op%ica_occ(1,iblk),op%ica_occ(2,iblk)).eq.0

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        nel = msamax+op%ica_occ(1,iblk)
        idxms = 0
        idxoff = 0
        idxoff_blk = 0
        do ms = msamax, -msamax, -2
          if (abs(ms+mst).gt.mscmax) cycle
          idxms = idxms+1
          do igam = 1, orb_info%nsym

            ! block offest and length:
            idxoff = mel%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle

            ! get current block
            if (.not.incore.and..not.blk_buf) then
              idoff = ffop%length_of_record*(ffop%current_record-1)
              call get_vec(ffop,buffer,
     &             idoff+idxoff+1,idoff+idxoff+lenblk)
              curblk => buffer(1:lenblk)
            else if (.not.incore) then
c              ioff = op%off_op_gmo(iblk)%gam_ms(igam,idxms)
              ! currently: idxoff should be valid here, as well
              curblk => ffop%buffer(idxoff+1:idxoff+lenblk)
            else
              curblk => bufmel(idxoff-idxoff0+1:idxoff-idxoff0+lenblk)
            end if

            ! reset offset within current buffer
            idxoff_blk = 0
        
            ! print out distributions
            first = .true.
            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
            if (ndis.eq.0) then
              write(luout,*) 'WARNING:'
              write(luout,*)
     &             ' op_info indicates no distribution at all'
              write(luout,*) ' skipping to next block'
              idxoff_blk = idxoff_blk+lenblk
              nwarn = nwarn+1
              cycle
            end if
            if (ndis.eq.1) then
c              write(luout,*) 'block contains only single distribution'
              if (scalar) then
                write(luout,'(2x,6f12.7)')
     &               curblk(idxoff_blk+1)
              else
                call wrt_mel_blk_seq(luwrt,curblk(idxoff_blk+1),
     &               mel,iblk,igam,idxms,1,
     &               nel,nindex,maxlen,str_info,orb_info)
              end if
              idxoff_blk = idxoff_blk+lenblk
              cycle
            end if
            distr_loop: do idx_dis = 1, ndis
            
              did = mel%off_op_gmox(iblk)%did(idx_dis,igam,idxms)

              call did2msgm(msd,igamd,did,
     &             op%ihpvca_occ(1,1,idx_occ),orb_info%nsym,njoined)
              if (.not.dagger) then
                scr(1:ngastp,1:2,1:njoined) =
     &               msd(1:ngastp,1:2,1:njoined)
                scr(1:ngastp,1:2,njoined+1:2*njoined) =
     &               igamd(1:ngastp,1:2,1:njoined)
              else
                scr(1:ngastp,1:2,1:njoined) =
     &               iocc_dagger_n(msd,njoined)
                  scr(1:ngastp,1:2,njoined+1:2*njoined) =
     &               iocc_dagger_n(igamd,njoined)
              endif

              lenblk =
     &             mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
              if (lenblk.eq.0) cycle
              lenprt = lenblk
              
              call wrt_mel_blk_seq(luwrt,curblk(idxoff_blk+1),
     &             mel,iblk,igam,idxms,idx_dis,
     &             nel,nindex,maxlen,str_info,orb_info)
              idxoff_blk = idxoff_blk+lenblk
            end do distr_loop
          end do                ! gam
        end do                  ! ms
        idx_occ = idx_occ+njoined
      end do

      ! Mark end of file with 0.
      write(luwrt) 0 

c dbg
c      ! Test to check that the integrals have been written properly.
c      print *,'TESTING'
c      nindex = 0
c      maxlen = 0
c      rewind luwrt
c      read(luwrt) nindex,maxlen
c      print *,'nindex, maxlen', nindex,maxlen
c
c      allocate(spins(nindex,maxlen),indices(nindex,maxlen),val(maxlen))
c
c      do 
c        read(luwrt) idxstr,indices(1:nindex,1:idxstr),
c     &       spins(1:nindex,1:idxstr),val(1:idxstr)
c        print *,'idxstr = ',idxstr
c        if(idxstr.eq.0) exit
c
c        do idum = 1,idxstr
c          write(luout,*)indices(1:nindex,idum),spins(1:nindex,idum),
c     &         val(idum)
c        enddo
c      enddo
c dbg

      if (nwarn.gt.0) then
        write(luout,*) '!!! There were ',nwarn,' warnings !!!'
        write(luout,*) 'look for "WARNING" in previous output!'
      end if

      if (.not.incore.and.close_again) call file_close_keep(ffop)
      call file_close_keep(fftemp)

      ifree = mem_flushmark()

      return
      end
