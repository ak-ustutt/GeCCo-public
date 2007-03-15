      subroutine help()

      implicit none

      include 'stdunit.h'

      write(luout,'(3(x,a,/))')
     &'usage: gecco.x [options] input_file',
     &'   options:',
     &'    -h,--help: print this help screen (and exit)'
      ! plus more to come

      return
      end
