      integer, parameter ::
     &     form_maxlen_comment = 256,
     &     form_maxlen_label = 32

      type formula
        character ::
     &     label*(form_maxlen_label),
     &     comment*(form_maxlen_comment)
        type(filinf) ::
     &     fhand
      end type formula
