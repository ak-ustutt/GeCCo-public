      interface
        function iocc_xdn(ixdn,iocc)
        include 'opdim.h'
        integer ::
     &     iocc_xdn(ngastp,2)
        integer, intent(in) ::
     &     ixdn,
     &     iocc(ngastp,2)
        end function iocc_xdn

        function iocc_overlap(iocc,dagi,jocc,dagj)
        include 'opdim.h'
        integer ::
     &       iocc_overlap(ngastp,2)
        logical, intent(in) ::
     &       dagi, dagj
        integer, intent(in) ::
     &       iocc(ngastp,2), jocc(ngastp,2)
        end function iocc_overlap

        function iocc_add(ifac,iocc,dagi,jfac,jocc,dagj)
        include 'opdim.h'
        integer :: iocc_add(ngastp,2)  
        logical, intent(in) ::
     &       dagi, dagj
        integer, intent(in) ::
     &       ifac, jfac,
     &       iocc(ngastp,2), jocc(ngastp,2)
        end function iocc_add

        function iocc_dagger(iocc_in)
        include 'opdim.h'
        integer ::
     &       iocc_dagger(ngastp,2)
        integer, intent(in) ::
     &       iocc_in(ngastp,2)
        end function iocc_dagger

        function iocc_dagger_n(iocc_in,njoined)
        implicit none
        include 'opdim.h'
        integer, intent(in) :: njoined,iocc_in(ngastp,2,njoined)
        integer :: iocc_dagger_n(ngastp,2,njoined)
        end function

        function irest_dagger_n(irest_in,njoined,ngas,nspin)
        implicit none
        integer, intent(in) :: njoined,ngas,nspin,
     &       irest_in(2,ngas,2,2,nspin,njoined)
        integer :: irest_dagger_n(2,ngas,2,2,nspin,njoined)
        end function

c        integer function idx_oplist(opname,ops,nops)
c        implicit none
c        include 'def_operator.h'
c        character, intent(in) ::
c     &       opname*(*)
c        integer, intent(in) ::
c     &       nops
c        type(operator), intent(in) ::
c     &       ops(nops)
c        end function

        logical function iocc_equal(iocc,dagi,jocc,dagj)
        implicit none
        include 'opdim.h'
        logical, intent(in) ::
     &      dagi, dagj
        integer, intent(in) ::
     &     iocc(ngastp,2), jocc(ngastp,2)
        end function

        logical function iocc_equal_n(iocc,dagi,jocc,dagj,njoined)
        implicit none
        include 'opdim.h'
        logical, intent(in) ::
     &      dagi, dagj
        integer, intent(in) ::
     &     njoined,
     &     iocc(ngastp,2,njoined), jocc(ngastp,2,njoined)
        end function

        logical function iocc_nonzero(iocc)
        implicit none
        include 'opdim.h'
        integer ::
     &     iocc(ngastp,2)
        end function

        logical function iocc_zero(iocc)
        implicit none
        include 'opdim.h'
        integer ::
     &     iocc(ngastp,2)
        end function

        logical function iocc_bound(cbound,iocc,dagi,jocc,dagj)      
        implicit none
        include 'opdim.h'
        character, intent(in) ::
     &       cbound*(*)
        logical, intent(in) ::
     &       dagi, dagj
        integer, intent(in) ::
     &       iocc(ngastp,2), jocc(ngastp,2)
        end function

c        integer function op_type(op)
c        implicit none
c        include 'opdim.h'
c        include 'def_operator.h'
c        type(operator), intent(in) ::
c     &       op
c	end function

      function irest_xdn(ixdn,irest,hpvxgas,ngas,nspin)
      implicit none
      include 'opdim.h'
      integer ::
     &     irest_xdn(2,ngas,2,2,nspin)
      integer, intent(in) ::
     &     ixdn, ngas, nspin, hpvxgas(ngas),
     &     irest(2,ngas,2,2,nspin)
      end function

      end interface
