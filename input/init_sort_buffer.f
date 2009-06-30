      subroutine init_sort_buffer(sbuffer,nbin,len_bin,max_chain)

      implicit none

      include 'def_sort_buffer.h'
      include 'ifc_memman.h'

      type(sort_buffer), intent(inout) ::
     &     sbuffer
      integer, intent(in) ::
     &     nbin, len_bin, max_chain

      type(sort_bin), pointer ::
     &     bin(:)
      integer ::
     &     ifree, ibin
      character(len=16) ::
     &     label

      sbuffer%nbins = nbin
      sbuffer%max_length = len_bin
      sbuffer%max_chain  = max_chain

      allocate(sbuffer%bin(nbin))

      bin => sbuffer%bin

      call mem_pushmark()
      
      ifree = mem_setmark('sbuffer')

      do ibin = 1, nbin
        
        bin(ibin)%length = 0
        write(label,'("val#",i6.6)') ibin
        ifree = mem_alloc_real(bin(ibin)%value,len_bin,label)
        write(label,'("idx#",i6.6)') ibin
        ifree = mem_alloc_int (bin(ibin)%index,len_bin,label)
        write(label,'("chn#",i6.6)') ibin
        ifree = mem_alloc_int (bin(ibin)%chain,max_chain+1,label)
        !  first index is for counting

      end do

      call mem_popmark()

      return
      end
