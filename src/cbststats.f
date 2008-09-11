c     function for calculating the full max t-statistic
      subroutine tmaxo(n,x,tss,sx,iseg,ostat,ibin)
      integer n,iseg(2)
      double precision x(n),tss,sx(n),ostat
      logical ibin

      integer i,j, sxmxi, tmaxi, tmaxj, nmj, nmji
c      double precision xsum,sx2,x1,x2,rij,rw,tij,xvar
      double precision rn,rj,absx,sxmx,bssmx,tmax

c     max and min partial sums and their locations
      integer ipsmin, ipsmax
      double precision psmin, psmax, psdiff, bsslim, rnjov1

      rn = dfloat(n)
      ostat = -0.5
      sx(1) = x(1)
      psmin = sx(1)
      ipsmin = 1
      psmax = sx(1)
      ipsmax = 1
      do 20 i = 2,n
         sx(i) = sx(i-1) + x(i)
         if (sx(i) .lt. psmin) then 
            psmin = sx(i)
            ipsmin = i
         endif
         if (sx(i) .gt. psmax) then 
            psmax = sx(i)
            ipsmax = i
         endif
 20   continue
      psdiff = psmax - psmin

      tmax = -0.5
c     compute the max statistic for segments of length j
      do 40 j = 2,(n-1)/2
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         if (ibin) then
            bsslim = rnjov1*(psdiff-0.5)**2
         else
            bsslim = rnjov1*psdiff**2
         endif
c         call dblepr("bss limit", 9, bsslim, 1)
         if (bsslim .lt. tmax) go to 60
c     sxmx is the maximum of abs partial sums x[i+1] + ... + x[i+j-1]
         sxmx = -0.5
         do 30 i = 1,n-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) then 
               sxmx = absx
c     i+1 is the start point for the current maximum
               sxmxi = i
            endif
 30      continue
         nmj = n - j
         do 35 i = 1, j
            nmji = nmj + i
            absx = abs(sx(nmji) - sx(i))
            if (sxmx.lt.absx) then 
               sxmx = absx
c     i+1 is the start point for the current maximum
               sxmxi = nmji
            endif
 35      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (tmax.lt.bssmx) then 
            tmax = bssmx
            tmaxi = sxmxi
            tmaxj = j
         endif
 40   continue
c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         rj = dfloat(j)
         sxmx = abs(sx(j))
         sxmxi = 0
         do 50 i = 1,n-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) then 
               sxmx = absx
c     i+1 is the start point for the current maximum
               sxmxi = i
            endif
 50      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmax.lt.bssmx) then 
            tmax = bssmx
            tmaxi = sxmxi
            tmaxj = j
         endif
      endif
c     convert statistic to t^2 form
 60   if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmax = tmax/(tss/rn)
      else
         if (tss.le.tmax+0.0001) tss = tmax + 1.0
         tmax = tmax/((tss-tmax)/(rn-2.0))
      endif

      ostat = tmax
      if (tmaxi+tmaxj.le.n) then
         iseg(1) = tmaxi
         iseg(2) = tmaxi+tmaxj
      else
         iseg(1) = tmaxi+tmaxj-n
         iseg(2) = tmaxi
      endif

      return
      end

c     function for calculating the full max t-statistic on permuted data
      double precision function tmaxp(n,tss,px,sx,ibin)
      integer n
      double precision tss,px(n),sx(n)
      logical ibin

      integer i, j, nmj, nmji
      double precision rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff,
     1     bsslim, rnjov1

      rn = dfloat(n)
      sx(1) = px(1)
      psmin = sx(1)
      psmax = sx(1)
      do 20 i = 2,n
         sx(i) = sx(i-1) + px(i)
         psmin = min(psmin, sx(i))
         psmax = max(psmax, sx(i))
 20   continue
      psdiff = psmax - psmin

      tmaxp = 0.0
c     compute the max statistic for segments of length j
      do 50 j = 2,k
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         if (ibin) then
            bsslim = rnjov1*(psdiff-0.5)**2
         else
            bsslim = rnjov1*psdiff**2
         endif
         if (bsslim .lt. tmaxp) go to 70
         sxmx = 0.0d0
         do 30 i = 1,n-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         nmj = n - j
         do 40 i = 1, j
            nmji = nmj + i
            absx = abs(sx(nmji) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 40      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (tmaxp.lt.bssmx) tmaxp = bssmx
 50   continue

c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         rj = dfloat(j)
         sxmx = 0.0
         do 60 i = 1,n/2
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 60      continue
         if (ibin) then
            bssmx = rn*(abs(sxmx)-0.5)**2/(rj*(rn-rj))
         else
            bssmx = rn*sxmx**2/(rj*(rn-rj))
         endif
         if (tmaxp.lt.bssmx) tmaxp = bssmx
      endif
 70   if (ibin) then
         if (tss.le.0.0001) tss = 1.0
         tmaxp = tmaxp/(tss/rn)
      else
         if (tss.le.tmaxp+0.0001) tss = tmaxp + 1.0
         tmaxp = tmaxp/((tss-tmaxp)/(rn-2.0))
      endif

      return
      end

c     function for the max (over small arcs) t-statistic on permuted data
      double precision function htmaxp(n,k,tss,px,sx,ibin)
      integer n,k
      double precision tss,px(n),sx(n)
      logical ibin

      integer i, j, nmj, nmji
      double precision rn, rj, absx, sxmx, bssmx, psmin, psmax, psdiff,
     1     bsslim, rnjov1

      rn = dfloat(n)
      sx(1) = px(1)
      psmin = sx(1)
      psmax = sx(1)
      do 20 i = 2,n
         sx(i) = sx(i-1) + px(i)
         psmin = min(psmin, sx(i))
         psmax = max(psmax, sx(i))
c         if (sx(i) .lt. psmin) psmin = sx(i)
c         if (sx(i) .gt. psmax) psmax = sx(i)
 20   continue
      psdiff = psmax - psmin

      htmaxp = 0.0d0
      do 50 j = 2,k
         rj = dfloat(j)
         rnjov1 = rn/(rj*(rn-rj))
         if (ibin) then
            bsslim = rnjov1*(psdiff-0.5)**2
         else
            bsslim = rnjov1*psdiff**2
         endif
         if (bsslim .lt. htmaxp) go to 60
         sxmx = 0.0d0
         do 30 i = 1,n-j
            absx = abs(sx(i+j) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 30      continue
         nmj = n - j
         do 40 i = 1, j
            nmji = nmj + i
            absx = abs(sx(nmji) - sx(i))
            if (sxmx.lt.absx) sxmx = absx
 40      continue
         if (ibin) then
            bssmx = rnjov1*(abs(sxmx)-0.5)**2
         else
            bssmx = rnjov1*sxmx**2
         endif
         if (htmaxp.lt.bssmx) htmaxp = bssmx
 50   continue
 60   if (ibin) then
         if (tss .le. 0.0001d0) tss = 1.0d0
         htmaxp = htmaxp/(tss/rn)
      else
         if (tss .le. htmaxp+0.0001d0) tss = htmaxp + 1.0d0
         htmaxp = htmaxp/((tss-htmaxp)/(rn-2.0d0))
      endif

      return
      end
