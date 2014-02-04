      subroutine help()

      implicit none

      include 'stdunit.h'

      write(lulog,'(7(x,a,/))')
     &'usage: gecco.x [options] input_file',
     &'   options:',
     &'    -l <logfile>, --logfile <logfile>: ',
     &'          print log-output to file logfile',
     &'          (leads to neat output on stdout)',
     &'    -h,--help: ',
     &'          print this help screen (and exit)'
      ! plus more to come

      return
      end
