      subroutine menulist
      character*60 promptlist,prompttype
      common /quizinfo/ nqst(600),lqst(600),iqfst(600),iqlst(600),
     -  maxq,init,nprompttype,promptlist(600),prompttype(600)
      character*55 submenu1(20)
      character*55 submenue(20)
      character*35 submenua(20)
      data submenu1 /
     -  '                                                       ',
     -  'optimization type                                      ',
     -  '                                                       ',
     -  'structure file format conversion type                  ',
     -  'trajectory file format conversion type                 ',
     -  'name conversion type                                   ',
     -  'trajectory/configuration stack packing/unpacking       ',
     -  'conformation manipulation                              ',
     -  'miscellaneous file creation                            ',
     -  'configuration analysis                                 ',
     -  'clustering source                                      ',
     -  9*'                                                       '/
      data submenua /
     -  'topology/geometry analysis         ',
     -  'bond tracking                      ',
     -  'atomic property calculation        ',
     -  'molecular property calculation     ',
     -  'RMSD calculation                   ',
     -  'distance analysis                  ',
     -  14*'                                   '/
      data submenue /
     -  'selecting option                                       ',
     -  '(additional) conformation transformation type          ',
     -  '                                                       ',
     -  'modification type                                      ',
     -  16*'                                                       '/
      call findmenu('run type',8,ix1)
      ix10=ix1
      ix1=ix1+1
      nn=0
      do while (promptlist(ix1)(1:3) .ne. '<Q>')
        write (6,1000) '   ',promptlist(ix1)
        call lastchar(submenu1(ix1-ix10),lc1,55)
        if (lc1 .gt. 1) then
          call findmenu(submenu1(ix1-ix10),lc1,ix2)
          ix20=ix2
          ix2=ix2+1
          do while (promptlist(ix2)(1:4) .ne. '****' .and.
     -              promptlist(ix2)(1:3) .ne. '<Q>')
            write (6,1000) '     ',promptlist(ix2)
            if (submenu1(ix1-ix10)(1:lc1) .eq.
     -          'configuration analysis') then
              call lastchar(submenua(ix2-ix20),lc2,35)
              if (lc2 .gt. 1) then
                call findmenu(submenua(ix2-ix20),lc2,ix3)
                ix3=ix3+1
                do while (promptlist(ix3)(1:4) .ne. '****' .and.
     -                    promptlist(ix3)(1:3) .ne. '<Q>')
                  write (6,1000) '       ',promptlist(ix3)
                  ix3=ix3+1
                end do
              end if
            else if (submenu1(ix1-ix10)(1:lc1) .eq.
     -          'conformation manipulation') then
c             Configuration edit submenus
              call lastchar(submenue(ix2-ix20),lc2,55)
              if (lc2 .gt. 1) then
                call findmenu(submenue(ix2-ix20),lc2,ix3)
                ix3=ix3+1
                do while (promptlist(ix3)(1:4) .ne. '****' .and.
     -                    promptlist(ix3)(1:3) .ne. '<Q>')
                  write (6,1000) '       ',promptlist(ix3)
                  ix3=ix3+1
                end do
              end if
            end if
            ix2=ix2+1
            if (ix2-ix20 .gt. 20) then
              write (6,2000)
              stop
            end if
          end do
        end if
        ix1=ix1+1
      end do
      return
1000  format(a,a)
2000  format(' PROGRAM ERROR: submenua or submenue array is too short',
     -  ' in subroutine menulist')
      end
