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

        integer function idx_oplist(opname,ops,nops)
        implicit none
        include 'def_operator.h'
        character, intent(in) ::
     &       opname*(*)
        integer, intent(in) ::
     &       nops
        type(operator), intent(in) ::
     &       ops(nops)
        end function

        logical function iocc_equal(iocc,dagi,jocc,dagj)
        implicit none
        include 'opdim.h'
        logical, intent(in) ::
     &      dagi, dagj
        integer, intent(in) ::
     &     iocc(ngastp,2), jocc(ngastp,2)
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

      end interface
