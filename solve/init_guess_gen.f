
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine init_guess_gen(guess_gen, maxtrials,
     &     me_diag, me_trv, nopt,
     &     op_info,str_info,strmap_info,orb_info)
*----------------------------------------------------------------------*
      
      implicit none
      include 'stdunit.h'
      include "def_guess_gen.h"
      include 'mdef_operator_info.h'
      include 'def_graph.h'
      include 'def_strinf.h'
      include 'def_strmapinf.h'
      include 'def_orbinf.h'
      
      type(guess_generator),intent(out)::
     &     guess_gen
      integer,parameter::
     &     ntest = 1000
      integer,intent(in)::
     &     maxtrials, nopt
      type(me_list_array)::
     &     me_diag(nopt),me_trv(nopt)

      integer::
     &     iopt,isign(nopt),ntrials(nopt)
      integer,pointer::
     &     idxlist(:,:)
      real(8),pointer::
     &     xlist_all(:), xlist(:,:)
      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info

      print *, "initializing guess_gen"
      
      allocate(guess_gen%idxlist_all(2,maxtrials),
     &     guess_gen%idxlist_ba(maxtrials))
      allocate(idxlist(maxtrials,nopt))
      allocate(xlist(maxtrials,nopt),xlist_all(maxtrials))

      print *, "allocated"
      do iopt = 1,nopt
        call find_nmin_list(xlist(1:maxtrials,iopt),
     &       idxlist(1:maxtrials,iopt),
     &       maxtrials,me_diag(iopt)%mel)
        ntrials(iopt) = min(maxtrials,me_diag(iopt)%mel%len_op)
      end do
      if (nopt .gt.1 )then
        call merge_min_lists(
     &       xlist_all,guess_gen%idxlist_all,guess_gen%ntrials_all,
     &       xlist,idxlist,
     &       nopt,maxtrials,ntrials,0)
      else
        xlist_all(1:maxtrials) = xlist(1:maxtrials,1)
        guess_gen%idxlist_all(2,1:maxtrials) = idxlist(1:maxtrials,1)
        guess_gen%idxlist_all(1,1:maxtrials) = 1
      end if

      do iopt =1,nopt
        isign(iopt) = me_trv(iopt)%mel%absym
      end do
      call create_inverted_spin_indexlist_h(
     &     guess_gen%idxlist_ba,
     &     idxlist,guess_gen%idxlist_all,
     &     guess_gen%ntrials_all,
     &     maxtrials,ntrials,me_diag,isign,nopt,
     &     op_info,str_info,strmap_info,orb_info)
      deallocate( idxlist, xlist, xlist_all)
      ! And initialize iguess
      guess_gen%iguess =0
      contains
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine create_inverted_spin_indexlist_h(
     &     idxlist_ba,idxlist,idxlist_all, ntrials_all,
     &     maxtrials,ntrials,
     &     me_dia,isign,nopt,
     &     op_info, str_info, strmap_info, orb_info)
*----------------------------------------------------------------------*
      implicit none
      
      integer,intent(out)::
     &     idxlist_ba(ntrials_all)

      integer,intent(in)::
     &     maxtrials,nopt,ntrials_all,
     &     ntrials(nopt),
     &     isign(nopt),
     &     idxlist(maxtrials,nopt),
     &     idxlist_all(2,ntrials_all)

      type(me_list_array)::
     &     me_dia(nopt)

      type(orbinf), intent(in) ::
     &     orb_info
      type(operator_info), intent(in) ::
     &     op_info
      type(strinf), intent(in) ::
     &     str_info
      type(strmapinf), intent(in) ::
     &     strmap_info
      
      integer::
     &     iopt, idx, iguess
      integer,pointer::
     &     idxscr(:)
      print *, "creating inverted spin_indexlist"

      allocate(idxscr(maxtrials))
      do iopt = 1, nopt
        if (isign(iopt).eq.0) cycle
        
        call symidx_ab(idxscr,
     &       idxlist(1,iopt),ntrials(iopt),me_dia(iopt)%mel,
     &       op_info,str_info,strmap_info,orb_info)
        
        if (ntest.ge.200) then
          write(lulog,*) 'ntrials(iopt),maxtrials: ',
     &         ntrials(iopt),maxtrials
          write(lulog,*) 'idxscr = '
          write(lulog,'(1x,5i4,x,5i4)') idxscr(1:ntrials(iopt))
        end if
        
        idx = 0
        do iguess = 1, ntrials_all
          if (idxlist_all(1,iguess) .ne.iopt) cycle
          idx = idx+1
          idxlist_ba(iguess) = idxscr(idx)
        end do
      end do
      deallocate(idxscr)
      return
      end subroutine
      
      end subroutine


*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      function generate_guess(guess_gen,me_trv,iopt,nopt,choice)
*----------------------------------------------------------------------*
      implicit none
      include 'stdunit.h'
      include "def_guess_gen.h"
      include 'mdef_me_list.h'
      logical::
     &     generate_guess
      character(*),parameter::
     &     i_am = "generate_guess"
      integer,parameter::
     &     ntest=1000
      
      type(guess_generator),intent(inout)::
     &     guess_gen
      integer,intent(out)::
     &     iopt
      integer,intent(in)::
     &     nopt,choice
      type(me_list_array)::
     &     me_diag(nopt),me_trv(nopt)

      
      integer::
     &     iroot,isign(nopt)
      guess_gen%iguess =  guess_gen%iguess +1
      
      do iopt =1,nopt
        isign(iopt) = me_trv(iopt)%mel%absym
      end do
      
      if (choice .ne. 0)then
        iopt = choice
      else
        iopt = guess_gen%idxlist_all(1,guess_gen%iguess)
      end if
      
      iroot = me_trv(1)%mel%fhand%current_record
      do while(.not. generate_guess_list_h(
     &     me_trv,isign,iopt,nopt,
     &     iroot, 
     &     guess_gen%idxlist_all,
     &     guess_gen%idxlist_ba,
     &     guess_gen%iguess,
     &     guess_gen%ntrials_all) )
        print *, iopt,guess_gen%iguess
        if (guess_gen%iguess.ge.guess_gen%ntrials_all)then
          generate_guess = .false.
          return
        end if 
        
        guess_gen%iguess =  guess_gen%iguess +1
        
        if (choice .ne. 0)then
          iopt = choice
        else
          iopt = guess_gen%idxlist_all(1,guess_gen%iguess)
        end if
         
      end do
      generate_guess = .true.
      contains
*----------------------------------------------------------------------*
      function generate_guess_list_h(
     &     me_trv,isign,iopt,nopt,
     &     iroot,
     &     idxlist_all,idxlist_ba,iguess,ntrials_all)
*----------------------------------------------------------------------*
      implicit none

      
      
      integer,intent(in)::
     &     isign(nopt),iguess,nopt,iopt,
     &     iroot,
     &     ntrials_all,
     &     idxlist_all(2,ntrials_all),
     &     idxlist_ba(ntrials_all)

      logical::
     &     generate_guess_list_h

      type(me_list_array)::
     &     me_trv(nopt)
      
      integer::
     &     idxset(2),
     &     nset,
     &     jopt

      real(8)::
     &     valset(2)


      nset = 0
      
      if (isign(iopt) .ne. 0) then
        if (idxlist_all(2,iguess).lt.abs(idxlist_ba(iguess)))then
          nset = 2
          idxset(1) = abs(idxlist_all(2,iguess))
          idxset(2) = abs(idxlist_ba(iguess))
          valset(1) = 1d0/sqrt(2d0)
          valset(2) = dble(isign(iopt))*
     &         dble(sign(1,idxlist_ba(iguess)))
     &         /sqrt(2d0)
        else if (idxlist_all(2,iguess).eq.abs(idxlist_ba(iguess))
     &         .and.isign(iopt).eq.+1) then
          if (idxlist_ba(iguess).lt.0)
     &         call quit(1,i_am,'unexpected case')
! set a single element
          nset = 1
          idxset(1) = abs(idxlist_all(2,iguess))
          valset(1) = 1d0
        end if
      else
! set a single element
        nset = 1
        idxset(1) = abs(idxlist_all(2,iguess))
        valset(1) = 1d0
      end if

      generate_guess_list_h=.false.
      if (nset.eq.0) return
      generate_guess_list_h=.true.
      
      if (ntest.ge.100) then
        write(lulog,*) 'iopt:   ',iopt
        write(lulog,*) 'idxset: ',idxset(1:nset)
        write(lulog,*) 'valset: ',valset(1:nset)
        write(lulog,*) '=> switching to iroot = ',iroot
      end if
      
      do jopt = 1, nopt
        if (jopt.eq.iopt) then
          call set_list(me_trv(jopt)%mel,idxset,valset,nset)
        else
          call zeroop(me_trv(jopt)%mel)
        end if
      end do
      
      return
      end function
      end function
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
      subroutine del_guess_gen(guess_gen)
*----------------------------------------------------------------------*
      implicit none
      include "def_guess_gen.h"
      type(guess_generator),intent(inout)::
     &     guess_gen
      deallocate(guess_gen%idxlist_all,
     &     guess_gen%idxlist_ba)
      return
      end subroutine
