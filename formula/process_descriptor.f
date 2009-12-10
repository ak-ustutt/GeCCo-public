      subroutine process_descriptor(vtx1,vtx2,occ_list,nlist,
     &                              dstring,maxlist)

      implicit none

      include 'stdunit.h'
      include 'def_word_list.h'

      integer, parameter ::
     &     ntest = 00
      character(len=18) ::
     &     i_am = 'process_descriptor'
      

      integer, intent(in) ::
     &     maxlist
      character(len=*), intent(in) ::
     &     dstring
      integer, intent(out) ::
     &     vtx1, vtx2, nlist, occ_list(4,2,maxlist)

      type(word_list) ::
     &     wlist
      type(word_list_entry) ::
     &     wl_pnt

      character(len=maxlen_word) ::
     &     cur_entry
      logical ::
     &     next
      integer ::
     &     nl_c, nl_a, icnt, il_a, il_c,
     &     occ_list_c(4,maxlist),
     &     occ_list_a(4,maxlist)

      logical, external ::
     &     next_word_list_entry

      if (ntest.ge.100) then
        call write_title(luout,wst_dbg_subr,i_am)
        write(luout,*) 'string: "',trim(dstring),'"'
      end if

      call init_word_list(wlist)
      call lex_line(wlist,dstring,' ,',' ')

      if (ntest.ge.100) then
        write(luout,*) 'lex''ed string:'
        call print_word_list(luout,wlist)
      end if

      nl_c = 1
      nl_a = 1
      occ_list_c = 0
      occ_list_a = 0

      icnt = 0
      do
        next = next_word_list_entry(cur_entry,wlist)
        icnt = icnt+1

        select case (icnt)
        case(1)
          read(cur_entry,*,err=100) vtx1
        case(2)
          if (len_trim(cur_entry).eq.0.or.cur_entry.eq.'-') then
            vtx2 = 0
          else
            read(cur_entry,*,err=100) vtx2
          end if
        case(3)
          call parse_hpvx_descr(occ_list_c,nl_c,cur_entry,maxlist)
        case(4)          
          call parse_hpvx_descr(occ_list_a,nl_a,cur_entry,maxlist)
        case default
          exit
        end select 

        if (.not.next) exit
      end do

      if (nl_c*nl_a.gt.maxlist)
     &     call quit(1,'process_descriptor','maxlen too small')

      nlist = 0
      do il_c = 1, nl_c
        do il_a = 1, nl_a
          nlist = nlist + 1
          occ_list(1:4,1,nlist) = occ_list_c(1:4,il_c)
          occ_list(1:4,2,nlist) = occ_list_a(1:4,il_a)
        end do
      end do

      return

 100  call quit(1,'process_descriptor','could not read integer from "'//
     &     trim(cur_entry)//'"')
      
      end subroutine

