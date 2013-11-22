*----------------------------------------------------------------------*
      subroutine import_2el_sort(ffout,ffpre,ffchain,len_list,len_bin)
*----------------------------------------------------------------------*
*
*     process chains from ffpre/ffchain and sort to ffout
*
*----------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'
      include 'ifc_memman.h'

      integer, parameter ::
     &     ntest = 00

      type(filinf), intent(inout) ::
     &     ffpre, ffchain, ffout
      integer, intent(in) ::
     &     len_list, len_bin

      integer ::
     &     lu_pre, lu_chain, ifree, idx, idxoff,
     &     nbin, ibatch, len_batch, ichain, idxrec
      real(8) ::
     &     val

      integer(8) ::
     &     reclen, idxbuf, lenchain, ii, len, np_rd, mc_rd
      integer ::
     &     n_per_rec, max_chain

      integer, pointer ::
     &     index(:), chain(:)
      real(8), pointer ::
     &     value(:), buffer(:)


      if (ntest.ge.100) then
        call write_title(lulog,wst_dbg_subr,'import_2el_sort')
      end if

      ifree = mem_setmark('2el_sort')

      lu_pre = ffpre%unit
      lu_chain = ffchain%unit
      rewind(lu_chain)
      read(lu_chain) reclen, np_rd, mc_rd
      n_per_rec = np_rd
      max_chain = mc_rd

      ifree = mem_alloc_real(buffer,len_bin,'buffer')
      ifree = mem_alloc_real(value,n_per_rec,'value')
      ifree = mem_alloc_int (index,n_per_rec,'index')
      ifree = mem_alloc_int (chain,max_chain,'chain')

      nbin = (len_list-1)/len_bin + 1

      idxoff = 0
      ! loop over chain-file = batches over target file
      do ibatch = 1, nbin
        ! initialize present batch of target file
        len_batch = min(len_bin,len_list-idxoff)
        buffer(1:len_batch) = 0d0

c dbg
c        print *,'ibatch: ',ibatch
c        print *,'  ',idxoff+1,'->',idxoff+len_batch
c dbg

        ! read new chain
        read(lu_chain) idxbuf,lenchain,chain(1:lenchain)

        if (idxbuf.ne.ibatch) then
          write(lulog,*) 'idxbuf = ',idxbuf
          write(lulog,*) 'ibatch = ',ibatch
          call quit(1,'import_2el_sort','how that?')
        end if

        ! process chain:
        do ichain = 1, lenchain

          ! read present record
          idxrec = chain(ichain)
          read(lu_pre,rec=idxrec) len,value(1:len),index(1:len)

          ! sort contents to target buffer
          do ii = 1, len
            val = value(ii)
            idx = index(ii)-idxoff
            buffer(idx) = buffer(idx) + val
          end do

        end do ! loop over chain

        ! write buffer to target file
        call put_vec(ffout,buffer,idxoff+1,idxoff+len_batch)

        idxoff = idxoff + len_batch
      end do ! loop over chain-file

      ifree = mem_flushmark('2el_sort')

      return
      end
