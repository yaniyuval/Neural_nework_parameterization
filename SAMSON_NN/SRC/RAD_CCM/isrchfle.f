      integer function isrchfle(n, array, inc, target)
c
c $Id: isrchfle.f,v 1.1.1.1 2007/01/11 22:58:18 cw Exp $
c $Author: cw $
c
      integer n,inc
      real array(*),target
      integer i
      if (n.le.0) then
         isrchfle = 0
         return
      end if
      ind = 1
      do i=1,n
         if (array(ind).le.target) then
            isrchfle = i
            return
         else
            ind = ind + inc
         end if
      end do
      isrchfle = n + 1
      return
      end

 
