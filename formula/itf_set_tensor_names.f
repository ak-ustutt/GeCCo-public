*---------------------------------------------------------------------*
      subroutine itf_set_tensor_names(itf_names,rename_list,nlist)
*---------------------------------------------------------------------*
*     parse the list rename_list(1:nlist) and set replacement for
*     GeCCo operator names to ITF tensor names
*
*     nlist is expected to consist of an even number of pair entries
*     odd entries are GeCCo operator names, even entries are IFT tensor
*     names. The latter allow the following special syntax:
*
*     "HAMGC","<ITF-HAM>" replaces the blocks of operator HAMGC by the
*     canonical names for the tensors representing matrix elements of the
*     Hamiltonian (f, J, K).
*
*     "OPGC","OPITF<RANK>" generates the ITF tensor names OPITF1, OPITF2
*     etc., depending on the particle rank of the corresponding block of OPGC
*
*     if no "<ITF-HAM>" is found in the list, "H","<ITF-HAM>" is
*     automatically assumed
*
*     andreas, Feb 21
*---------------------------------------------------------------------*
      implicit none

      include 'stdunit.h'
      include 'opdim.h'
      include 'def_contraction.h'
      include 'def_itf_contr.h'
      include 'def_target.h'

      character(len=18), parameter ::
     &     i_am = 'itf_set_tensor_names'
      integer, parameter ::
     &     ntest = 100
      
      integer, intent(in) ::
     &     nlist
      character(len=len_command_par), intent(in) ::
     &     rename_list(nlist)
      type(tensor_names), intent(inout) ::
     &     itf_names

      logical ::
     &     found
      integer ::
     &     ii, jj, idx, nalloc
      
! sanity checks:
      if (mod(nlist,2)==1) then
        do ii = 1, nlist-1, 2
          write(lulog,'(1x,a,"->",a)') rename_list(ii:ii+1)
        end do
        write(lulog,'(1x,a)') rename_list(nlist)
        call quit(0,i_am,'expected even number of entries for RENAME')
      end if

! check if <ITF-HAM> is present
      found = .false.
      do ii = 1, nlist-1, 2
        found = found.or.rename_list(ii+1)(1:9)=="<ITF-HAM>"
      end do

      nalloc = nlist/2
      if (.not.found) nalloc = nalloc+1

      allocate(itf_names%gc_itf_rename(2,nalloc))
      allocate(itf_names%rename_type(nalloc))
      do ii = 1, nalloc
        itf_names%gc_itf_rename(1:2,ii)(1:MAXLEN_BC_LABEL) = " "
      end do
      
      jj=1
      if (.not.found) then
        itf_names%gc_itf_rename(1,jj) = "H"
        itf_names%gc_itf_rename(2,jj) = ""
        itf_names%rename_type(jj) = RENAME_HAMILTONIAN
        jj = jj+1
      end if

      do ii = 1, nlist-1, 2
        if (rename_list(ii+1)(1:9)=="<ITF-HAM>") then
          itf_names%gc_itf_rename(1,jj) = trim(rename_list(ii))
          itf_names%gc_itf_rename(2,jj) = ""
          itf_names%rename_type(jj) = RENAME_HAMILTONIAN
          jj = jj+1
        else if (index(rename_list(ii+1),'<')>0) then
          idx = index(rename_list(ii+1),'<')
          if (rename_list(ii+1)(idx:idx+5)=="<RANK>") then
            itf_names%gc_itf_rename(1,jj) = trim(rename_list(ii))
            itf_names%gc_itf_rename(2,jj) = rename_list(ii+1)(1:idx-1)
            itf_names%rename_type(jj) = RENAME_ADD_RANK
            jj = jj+1
          else
            write(lulog,'(1x,a,"->",a)') rename_list(ii:ii+1)
            call quit(0,i_am,'Could not intepret content after "<"')
          end if
        else
          itf_names%gc_itf_rename(1,jj) = trim(rename_list(ii))
          itf_names%gc_itf_rename(2,jj) = trim(rename_list(ii+1))
          itf_names%rename_type(jj) = RENAME_BASIC
          jj = jj+1
        end if

      end do

      itf_names%nrename = nalloc

      if (ntest.ge.100) then
        write(lulog,'(" At the end of ",a)') i_am
        write(lulog,'(" Number of cases: ",i5)') nalloc
        do jj = 1, nalloc
          select case(itf_names%rename_type(jj))
          case(RENAME_BASIC)
            write(lulog,'(" Rename ",a12," -> ",a12)')
     &           itf_names%gc_itf_rename(1,jj),
     &           itf_names%gc_itf_rename(2,jj)
          case(RENAME_HAMILTONIAN)
            write(lulog,'(" Rename ",a12," -> f, J, K")')
     &           itf_names%gc_itf_rename(1,jj)
          case(RENAME_ADD_RANK)
            write(lulog,'(" Rename ",a12," -> ",a12, " and add rank")')
     &           itf_names%gc_itf_rename(1,jj),
     &           itf_names%gc_itf_rename(2,jj)
          end select
        end do
      end if
      
      end subroutine
