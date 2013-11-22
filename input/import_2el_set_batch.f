*----------------------------------------------------------------------*
      subroutine import_2el_set_batch(len_bin,n_rec,nbin,maxchain,
     &                 npass,nbin_per_pass,
     &                 len_list,max_mem,len_rec,
     &                 len_bin_min,len_bin_max)
*----------------------------------------------------------------------*
*     set up bin lengthes for presort algorithm
*----------------------------------------------------------------------*

      implicit none

      include 'stdunit.h'

      integer, parameter ::
     &     ntest = 100

      integer, intent(in) ::
     &     len_list, max_mem, len_rec,
     &     len_bin_min, len_bin_max
      integer, intent(out) ::
     &     len_bin, n_rec, nbin, maxchain, npass, nbin_per_pass

      integer ::
     &     error, need, maxbuffer

      len_bin = len_bin_min

      error = 0
      do
        n_rec = (len_rec-1)/2
        maxchain = (2*len_bin-1)/n_rec+1  ! 2: up to 2 values may contribute
        need = len_rec + maxchain

        nbin = (len_list-1)/len_bin + 1
        maxbuffer = max_mem/need

        if (ntest.ge.100) then
          write(lulog,*) 'current state:'          
          write(lulog,*) ' max_mem   = ',max_mem
          write(lulog,*) ' len_list  = ',len_list
          write(lulog,*) ' len_rec   = ',len_rec
          write(lulog,*) ' len_bin   = ',len_bin
          write(lulog,*) ' maxchain  = ',maxchain
          write(lulog,*) ' maxbuffer = ',maxbuffer
        end if          

        if (maxbuffer.eq.0.and.maxchain.lt.9*len_rec) then
          error = 1
          exit
        end if

        npass = (nbin-1)/maxbuffer + 1
        nbin_per_pass = min(maxbuffer,nbin)

        if (ntest.ge.100) then
          write(lulog,*) ' npass         = ',npass
          write(lulog,*) ' nbin_per_pass = ',nbin_per_pass
        end if          

        if (npass.eq.1.or.len_bin*2.gt.len_bin_max) exit

        len_bin = len_bin*2

      end do

      if (error.eq.1) then
        write(lulog,*) 'max_mem = ',max_mem
        write(lulog,*) 'need    = ',need
        call quit(1,'import_2el_set_batch',
     &       'unexpectedly small memory')
      end if

      return
      end
