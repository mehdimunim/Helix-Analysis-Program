      subroutine getindex(nlist,ip,iv,i1,i2)
c     Generate the vertex indices corresponding the (nlist+ip)th list element
      iv=nlist+ip
      i1=nlist+mod(ip,3)+1
      i2=nlist+mod(ip+1,3)+1
      return
      end
