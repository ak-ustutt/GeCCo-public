*----------------------------------------------------------------------*
      subroutine set_frequency2(mel,mel_freq,fac)
*----------------------------------------------------------------------*
*     set frequency assigned to ME-list
*     get value from mel_freq
*     (we do *not* use the "species" descriptor here, sign is set by
*     user)
*     matthias, 2008 -> alternative version: andreas, 2025
*----------------------------------------------------------------------*

      implicit none

      include 'def_operator.h'
      include 'def_filinf.h'
      include 'def_me_list.h'
      include 'stdunit.h'

      type(me_list), intent(inout) ::
     &     mel, mel_freq
      real(8), intent(in) ::
     &     fac

      integer, parameter ::
     &     ntest = 00

      integer ::
     &     iprint, idoff_frq
      real(8) ::
     &     freq
      logical ::
     &     open_close_frq

      iprint = max(iprlvl, ntest)

      ! read frequency from mel_freq:
      open_close_frq = mel_freq%fhand%unit.le.0
      if (open_close_frq) then
         call file_open(mel_freq%fhand)
      endif

      idoff_frq = mel_freq%fhand%length_of_record
     &            *(mel_freq%fhand%current_record-1)
      if (mel_freq%fhand%buffered) then
         if (idoff_frq.gt.0) call quit(1,'set_frequency2','check this!')
         freq = mel_freq%fhand%buffer(1)
      else
         call get_vec(mel_freq%fhand,freq,idoff_frq+1,idoff_frq+1)
      end if

      if (open_close_frq)
     &     call file_close_keep(mel_freq%fhand)

      mel%frequency = freq*fac

      if (iprint.ge.10) 
     &      write(lulog,*)
     &           'Frequency associated with ',trim(mel%label),
     &           ': ',mel%frequency

      return
      end
