      subroutine import_intm_fc(mel,name,str_info,orb_info)

      implicit none

      integer, parameter ::
     &     ntest =  00

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_orbinf.h'
      include 'par_opnames_gen.h'
      include 'mdef_operator_info.h'
      include 'ifc_memman.h'

      type(me_list), intent(inout) ::
     &     mel
      character(*), intent(in) ::
     &     name
      type(strinf), intent(in) ::
     &     str_info
      type(orbinf), intent(in), target ::
     &     orb_info
      
      type(filinf) ::
     &     fftemp
      type(filinf), pointer ::
     &     ffop
      type(operator), pointer ::
     &     op

      integer(8) ::
     &     nindex,maxlen,lenbuf,nel8

      integer ::
     &     lutemp, mst, njoined, idxstr, iblk, idx_occ,
     &     mscmax, msamax, nel, idxms, idxoff, idxoff_blk, ms, lenblk,
     &     ndis, igam, idx_dis, lendis, did, ngam, ngas, idoff, mmax,
     &     ifree, idxbuf
      logical ::
     &     fexist, close_again, first, blk_buf

      integer ::
     &     msd(ngastp,2,mel%op%njoined),igamd(ngastp,2,mel%op%njoined),
     &     lexlscr(8,3), idorb(8), idorb2(8), idspn(8), idspc(8)
      integer, pointer ::
     &     occ(:,:,:), idx_graph(:,:,:), mostnd(:,:,:), igamorb(:),
     &     ngas_hpv(:), idx_gas(:), igas_restr(:,:,:,:,:)

      integer(2), pointer ::
     &     spins(:,:), indices(:,:)
      real(8),pointer ::
     &     val(:), curblk(:), buffer_reo(:)

      real(8) ::
     &     cpu0, sys0, wall0, cpu, sys, wall

      call atim_csw(cpu0,sys0,wall0)

      ngam = orb_info%nsym
      ngas = orb_info%ngas
      mostnd => orb_info%mostnd
      igamorb => orb_info%igamorb
      ngas_hpv => orb_info%ngas_hpv
      idx_gas => orb_info%idx_gas
      ffop => mel%fhand
      mst = mel%mst
      op => mel%op
      njoined = op%njoined
      igas_restr => str_info%igas_restr

      if(ntest.ge.100)then
        write(lulog,*) 'Frozen-core import of F12-intermediates'
        write(lulog,*) 'Intermediate = ',trim(op%name)
      endif

      ! Find maximum length needed for reorder buffer.
      mmax = 0
      do iblk = 1, op%n_occ_cls
        if (op%formal_blk(iblk)) cycle
        idx_occ = (iblk-1)*njoined+1
        occ => op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)

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
      ifree = mem_alloc_real(buffer_reo,mmax,'buffer_reo')

      ! open the file to contain the ME-list
      close_again = .false.
      if (ffop%unit.le.0) then
        close_again = .true.
        call file_open(ffop)
      end if

      ! Open the file containing the integrals and indices.
      inquire(file=trim(name),exist=fexist)
      if(fexist)then
        call file_init(fftemp,trim(name),ftyp_sq_unf,0)
        call file_open(fftemp)
        lutemp = fftemp%unit
        rewind lutemp
      else
        call quit(1,'import_intm_fc','file nonexistent')
      endif

      ! Read in initial info about block sizes and allocate arrays.
      read(lutemp) nindex,maxlen
c dbg
c      print *,'indices: ',nindex
c      print *,'block length: ',maxlen
c dbg
      allocate(spins(nindex,maxlen),indices(nindex,maxlen),val(maxlen))
      idxbuf = maxlen+1 ! indicate that we have to read a buffer
      lenbuf = 0

c      if((trim(op%name).eq.op_p_inter.and.nindex.ne.4).or.
c     &     (trim(op%name).eq.op_z_inter.and.nindex.ne.6))
c     &     call quit(1,'import_intm_fc','wrong nindex')

      ! Loop over blocks of the intermediate.
      blk_loop: do iblk = 1, op%n_occ_cls
        if (op%formal_blk(iblk)) cycle
        idx_occ = (iblk-1)*njoined+1
        
        occ => op%ihpvca_occ(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)
        idx_graph =>
     &       mel%idx_graph(1:ngastp,1:2,idx_occ:idx_occ+njoined-1)

        blk_buf = ffop%buffered
        if (blk_buf) blk_buf = blk_buf.and.ffop%incore(iblk).gt.0

        mscmax = op%ica_occ(1,iblk)
        msamax = op%ica_occ(2,iblk)
        nel = msamax+op%ica_occ(1,iblk)
        if(nel.gt.nindex) call quit(1,'import_intm_fc','nel>nindex')

        idxms = 0
        idxoff = 0
        idxoff_blk = 0

        ms_loop: do ms = msamax, -msamax, -2
          if (abs(ms+mst).gt.mscmax) cycle
          idxms = idxms+1

          igam_loop: do igam = 1, orb_info%nsym
              
            ! block offest and length:
            idxoff = mel%off_op_gmo(iblk)%gam_ms(igam,idxms)
            lenblk = mel%len_op_gmo(iblk)%gam_ms(igam,idxms)
            if (lenblk.eq.0) cycle
            ! Point to block.
            if (.not.blk_buf) then
              curblk => buffer_reo
            else 
              curblk => ffop%buffer(idxoff+1:idxoff+lenblk)
            endif

            first = .true.
            ndis = mel%off_op_gmox(iblk)%ndis(igam,idxms)
            if (ndis.eq.0)
     &           call quit(1,'import_intm_fc','ndis=0?')

            ! Loop over distribution.
            distr_loop: do idx_dis = 1, ndis
                  
              lendis =
     &           mel%len_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
              if (lendis.eq.0) cycle
              idxoff_blk =
     &             mel%off_op_gmox(iblk)%d_gam_ms(idx_dis,igam,idxms)
     &             - idxoff
              did = mel%off_op_gmox(iblk)%did(idx_dis,igam,idxms)
              call did2msgm(msd,igamd,did,occ,ngam,njoined)

              call inner_loop(curblk(idxoff_blk+1:))

            enddo distr_loop
            idoff = ffop%length_of_record*(ffop%current_record-1)
            call put_vec(ffop,curblk,
     &           idoff+idxoff+1,idoff+idxoff+lenblk)

          enddo igam_loop
        enddo ms_loop
      enddo blk_loop

      call file_close_keep(fftemp)
      if(close_again) call file_close_keep(ffop)

      call atim_csw(cpu,sys,wall)

      if (iprlvl.ge.5) 
     &     call prtim(lulog,'time in intermediate import',
     &     cpu-cpu0,sys-sys0,wall-wall0)

      return

      contains

      subroutine inner_loop(curdisblk)

      implicit none 

      include 'hpvxseq.h'

      real(8), intent(inout) ::
     &     curdisblk(*)

      integer ::
     &     idxstr, jdx
      logical ::
     &     match

      logical, external ::
     &     next_tupel_ca

c dbg
c      do 
c        read(lutemp) lenblk,indices(1:nindex,1:lenblk),
c     &       spins(1:nindex,1:lenblk),val(1:lenblk)
c        print *,'idxstr = ',lenblk
c        if(lenblk.eq.0) exit
c
c        do idx = 1,lenblk
c          write(lulog,*)indices(1:nindex,idx),spins(1:nindex,idx),
c     &         val(idx)
c        enddo
c      enddo
c      print *,'occ',occ
c      print *,'ms dist',msd
c      print *,'gam dist',igamd
c      print *,'nel, njoined, ngam, ngas',nel,njoined,ngam,ngas
c      print *,'ihpvseq',hpvxseq
c dbg

      first = .true.
      idxstr = 0
      do while(next_tupel_ca(idorb,idspn,idspc,
     &     nel,njoined,occ,
     &     idx_graph,
     &     msd,igamd,first,
     &     igas_restr,
     &     mostnd,igamorb,
     &     ngam,ngas,
     &     ngas_hpv,idx_gas,
     &     hpvxseq,lexlscr))
        first = .false.
        idxstr = idxstr+1

        do

          ! read next buffer, if necessary
          if (idxbuf.gt.lenbuf) then
            nel8 = nindex
            read(lutemp)  lenbuf,indices(1:nel8,1:lenbuf),
     &         spins(1:nel8,1:lenbuf),val(1:lenbuf)
c            read(lutemp)  nel8,lenbuf,indices(1:nel8,1:lenbuf),
c     &         spins(1:nel8,1:lenbuf),val(1:lenbuf)
c dbg
c            print *,'read new block, length = ',nel8,lenbuf
c            print *,'idxstr = ',idxstr
c dbg
            if (lenbuf.le.0)
     &           call quit(1,'import_intm_fc',
     &           'arrived at end of file, but import seems'//
     &           ' not complete !?')
            idxbuf = 0
          end if
          idxbuf = idxbuf+1

          match = .true.
          ! Check indices.
          do jdx = 1, nel8
            match = match.and.idorb(jdx).eq.indices(jdx,idxbuf)
          enddo
          if(.not.match) cycle

          ! Check spins.
          do jdx = 1, nel8
            match = match.and.idspn(jdx).eq.spins(jdx,idxbuf)
          enddo
          if(.not.match)cycle

          ! If correct indices and spins, copy to the buffer.
          curdisblk(idxstr) = val(idxbuf)
          exit
        end do

      end do

      return
      end subroutine
      end
