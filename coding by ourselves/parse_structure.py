
def parse_structure(filename):
    """
    Parse PDB entry file
    to get the CA, C, O, N and H from the backbone chain
    Returns: the corresponding coordinates for each residue
    """
    
    do ia = n1, nslt
            atnam(1: lnam) = line(index(ia))(inamcol1: inamcol2)
            if (lnam .eq. 4) atnam(5: 8) = '    '
            if (lnam .gt. 4) call leftadjustn(atnam, atnam, lnam)
            if (atnam(1: 4) .eq. 'C   ' . or . atnam(1: 4) .eq. ' C   ') then
              icfound = 1
              ixc(nres+1) = ia
            else if (atnam(1: 4) .eq. 'O   ' . or . atnam(1: 4) .eq. ' O   ')
         -           then
              iofound=1
              ixo(nres+1)=ia
            else if (atnam(1:4) .eq. 'N   ' .or. atnam(1:4) .eq. ' N   ')
         -           then
              infound=1
              ixn(nres+1)=ia
            else if (atnam(1:4) .eq. 'CA  ' .or. atnam(1:4) .eq. ' CA  ')
         -           then
              iafound=1
              ixa(nres+1)=ia
            else if (atnam(1:4) .eq. 'H   ' .or. atnam(1:4) .eq. ' H   '
         -    .or. atnam(1:4) .eq. ' D  ' .or. atnam(1:4) .eq. 'D   ' .or.
         -    atnam(1:4) .eq. 'HN  ' .or. atnam(1:4) .eq. ' HN ') then
              ihfound=1
              call trnsfr(ch(1,nres+1),c(1,ia),3)
            end if
            if (ia .eq. nslt) then
              iresn=-1
