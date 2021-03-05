      subroutine regularpdb(inp,reg,itofrom)
      character*4 inp,reg,regg,inpl
      regg='    '
      call leftadjust4(inp,inpl)
c     write (06,*) 'REGPDB inp=',inp,' inpl=',inpl,' itofrom=',itofrom
      if (itofrom .eq. 1) then
c       To PDB
        if (inpl(1:1) .eq. 'C' .or. inpl(1:1) .eq. 'O' .or.
     -    inpl(1:1) .eq. 'N' .or. inpl(1:1) .eq. 'S' .or.
     -    inpl(1:1) .eq. 'P') then
          regg(2:4)=inpl(1:3)
          if (inpl(4:4) .ne. ' ') then
            if (idigit(inpl(4:4),1) .eq. 1) then
              regg(1:1)= inpl(4:4)
            else if (idigit(inpl(3:3),1) .eq. 1) then
              regg(1:1)= inpl(3:3)
              regg(4:4)= inpl(4:4)
            else if (idigit(inpl(2:2),1) .eq. 1) then
              regg(1:1)= inpl(3:3)
              regg(3:3)= inpl(3:3)
              regg(4:4)= inpl(4:4)
            else
              write (6,1000) ' ',inpl
            end if
          end if
        else if (inpl(1:1) .eq. 'H' .or. inpl(1:1) .eq. 'D') then
          if (inpl .eq. 'H   ') then
            regg=' H  '
          else if (inpl .eq. 'D   ') then
            regg=' D  '
          else if (inpl(3:4) .eq. '  ') then
c           Two-character H name
            regg(2:3)=inpl(1:2)
          else if (inpl(4:4) .eq. ' ') then
c           Three-character H name
            regg(1:1)=' '
            regg(2:4)=inpl(1:3)
          else
c           Four-character H name
            if (idigit(inpl(4:4),1) .eq. 1) then
c             HXXd -> dHXX
              regg(1:1)=inpl(4:4)
              regg(2:4)=inpl(1:3)
            else if (idigit(inpl(2:2),1) .eq. 1 .and.
     -               idigit(inpl(3:3),1) .eq. 1) then
c             HddX -> dHXd
              regg(1:1)=inpl(2:2)
              regg(2:2)=inpl(1:1)
              regg(3:3)=inpl(4:4)
              regg(4:4)=inpl(3:3)
            else if (idigit(inpl(2:2),1) .eq. 1) then
c             HdXX -> dHXX
              regg(1:1)=inpl(2:2)
              regg(2:2)=inpl(1:1)
              regg(3:4)=inpl(3:4)
            else if (idigit(inpl(3:3),1) .eq. 1) then
c             HXdX -> dHXX
              regg(1:1)=inpl(3:3)
              regg(2:2)=inpl(1:1)
              regg(3:3)=inpl(2:2)
              regg(4:4)=inpl(4:4)
            else
              write (6,1000) ' H ',inpl
              regg=inpl
            end if
          end if
        else if (idigit(inpl(1:1),1) .eq. 1 .and.
     -           inpl(4:4) .eq. ' ') then
           regg(2:3)=inpl(2:3)
           if (inpl(3:4) .eq. '  ') then
             regg(3:3)=inpl(1:1)
           else
             regg(4:4)=inpl(1:1)
           end if
        else
          regg=inpl
        end if
      else
c       From PDB
c       Just move the number away from the first position, if any
        if (idigit(inpl(1:1),1) .eq. 1) then
          if (inpl(3:4) .eq. '  ' .or. inpl(3:4) .eq. "''" .or.
     -        inpl(3:4) .eq. '**' .or. inpl(3:4) .eq. "' " .or.
     -        inpl(3:4) .eq. '* ') then
c           Two-character H name
            regg(1:1)=inpl(2:2)
            regg(2:2)=inpl(1:1)
          else if (inpl(4:4) .eq. ' ' .or. inpl(4:4) .eq. '*' .or.
     -             inpl(4:4) .eq. "'") then
c           Three-character H name
            regg(1:2)=inpl(2:3)
            regg(3:3)=inpl(1:1)
          else
c           Four-character H name
c           Permute cyclicly
            regg(1:3)=inpl(2:4)
            regg(4:4)=inpl(1:1)
          end if
        else
          regg=inpl
        end if
      end if
      reg=regg
c     write (06,*) 'REGPDB reg=',reg
      return
1000  format(' ERROR: can not regularize 4-character',a,
     -  'name withot a digit:',a)
      end
