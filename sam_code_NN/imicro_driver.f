      subroutine imicro_driver(nk,tabsz,presz,qmixz,tssz,pssz,rhoz,
     &     dostats)

      implicit none

      COMMON/VTERM/VTS,VTG,VTR
      COMMON/MVAR/T,PRESR,QMIX,TSS,PSS,RHOK

      real vts, vtr, vtg, t, presr, qmix(6), tss(7), pss(26), rhok

C     inputs/outputs
      integer nk
      real    tabsz(nk)
      real    presz(nk)
      real    qmixz(6,nk)
      real    tssz(7,nk)
      real    pssz(35)
      real    rhoz(nk)
      logical dostats
      
C     local variables
      integer k, i

      if (dostats) then
         do i = 1,35
            pssz(i) = 0.
         end do
      end if

      do k = 1,nk
         
         t = tabsz(k)
         rhok = rhoz(k)
         presr = presz(k)
         
         do i = 1,6
            qmix(i) = qmixz(i,k)
         end do

         do i = 1,7
            tssz(i,k) = 0.
         end do

         if (dostats) then
            do i = 1,26
               pss(i) = 0.
            end do
         end if

         call imicro()

         do i = 1,7
            tssz(i,k) = tss(i)
         end do

         if (dostats) then
            do i = 1,26
               pssz(i) = pssz(i) + pss(i)
            end do
            if (qmix(6).gt.0.) then
               pssz(27) = pssz(27) + vtr 
               pssz(28) = pssz(28) + 1.
            end if
            if (qmix(5).gt.0.) then
               pssz(29) = pssz(29) + vtg
               pssz(30) = pssz(30) + 1.
            end if
            if (qmix(4).gt.0.) then
               pssz(31) = pssz(31) + vts
               pssz(32) = pssz(32) + 1.
            end if
            if (qmix(3).gt.0.) then
c no ice fall velocity in Kreugers scheme
c              pssz(33) = pssz(33) + vts
               pssz(34) = pssz(34) + 1.
            end if
            if (qmix(2).gt.0.) then
               pssz(35) = pssz(35) + 1.
            end if
         end if
         
      end do

      return
      end
