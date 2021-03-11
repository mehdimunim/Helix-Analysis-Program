# Will be possibly higly modified
# Need to interpret the code
def (q, lenqin, iyn, idefans, ians, ihelp, itip):
    """

    ***
    Parameters:

    q:
    lenqin:
    iyn:
    idefans:
    ians:
    ihelp:
    itip


    """

    ians = 0
    if (idefans == - 1:
        defans=' [n] '
        lendef=5
    elif(idefans == 1):
        defans=" [y] "
        lendef=5
    else:
        defans=" "
        lendef=1
    lenq=lenqin
    pline[0:1]=" "
    pline[1:lenq + 1]=q[0:lenq]
    pline[lenq+2:lenq+6]=" \(y//n"
    lenq += 6
    if (ihelp > 0):
        pline(lenq+1: lenq+2)=" ?"
        lenq += 2
    if (itip > 0):
        pline[lenq+1:lenq+2]=" +"
        lenq += 2
    lenq += 1
    pline[lenq:lenq+1]="\)"
    pline[lenq+1:lenq+lendef]=defans[:lendef]
    # 100   call writeline(6, pline, 1, lenq+lendef, 1)
    ans=" "
    read(5, 1001, end=99, err=99) ans
    if (ans == "?"):
        if (ihelp == 0):
            print("Sorry, no help on this")
            # go to 100
        else:
            call explanation(ihelp, 0)
            # go to 100
    if (ans == "+ "):
        if (itip == 0):
            print("Sorry, no tip on this")
            # go to 100
        else:
            # call explanation(0,itip)
            #  go to 100
    if (ans.lower() != "n" and ans.lower() != "y"):
        # call lastchar(ans,ilc,1)
        if (ilc == 1:
            print("Invalid answer - please, answer yes or no")
            # go to 100
        elif(idefans != 1 and idefans != -1):
            print("Pls answer y or no")
            # go to 100
        else:
            if (idefans == -1):
                ans="n"
            elif(idefans == 1):
                ans="y"
    if(logfile > 0):
        with open(logfile, "w+") as f:
            f.write(ans)
    if (ans.lower() == "y"):
        ians=1
    if (iyn == 0):
        ians=1 - ians
    return ians
