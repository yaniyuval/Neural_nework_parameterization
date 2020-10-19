      subroutine whenfgt(n,array,inc,target,index,nval)
c
c $Id: whenfgt.f,v 1.1.1.1 2007/01/11 22:58:18 cw Exp $
c $Author: cw $
c
	dimension index(*), array(*)
	ina=1
	nval=0
	if(inc .lt. 0) ina=(-inc)*(n-1)+1
	do 100 i=1,n
	    if(array(ina) .gt. target) then
	    nval=nval+1
	    index(nval)=i
	    end if
	    ina=ina+inc
 100    continue
      return
      end           
 
