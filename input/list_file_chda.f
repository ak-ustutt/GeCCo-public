      subroutine list_file_chda(fflist,ffchain)
*
*     print contents of DA file with chained records
*     files are assumed open
*
      implicit none

      include 'stdunit.h'
      include 'def_filinf.h'

      type(filinf), intent(in) ::
     &     fflist, ffchain
      
      integer(8) ::
     &     len, iii, icnt, maxchain, maxint, reclen, idxbuf, lenchain
      integer ::
     &     ich, un, unch

      integer, pointer ::
     &     idx(:)
      integer, pointer ::
     &     chain(:)
      real(8), pointer ::
     &     val(:)

      unch = ffchain%unit
      un   = fflist%unit

      rewind(unch)
      read(unch) reclen, maxint, maxchain

      
      allocate(chain(maxchain))
      allocate(idx(maxint),val(maxint))

      icnt = 0
      do
        read(unch) idxbuf,lenchain,chain(1:lenchain)
        if (idxbuf.le.0) exit
        write(luout,*) 'chain # ',idxbuf
c dbg
        write(luout,*) 'length: ',lenchain
        write(luout,*) 'records: ',chain(1:lenchain)
c dbg

        do ich = 1, lenchain
c dbg
          if (chain(ich).eq.0) write(luout,*)
     &         'warning: ignoring 0 record'
          if (chain(ich).eq.0) cycle
c dbg
          read(un,rec=chain(ich)) len,val(1:len),idx(1:len)          
          write(luout,*) 'record #',chain(ich),' length = ',len

          do iii = 1, len
            write(luout,'("M(",i8,") = ",g20.14)') idx(iii),val(iii)
          end do
        end do
      end do

      deallocate(val,idx)
      deallocate(chain)
      return
      end
