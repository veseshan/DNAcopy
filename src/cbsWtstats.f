c     these are the subroutines to do the weigthed version of CBS
c     which is useful in order to merge data from multiple platforms
c     --------------------------------------------------------------
c         This is relevant only for log-ratio not binary data
c     --------------------------------------------------------------
c     function for calculating the full max weighted t-statistic
      subroutine wtmaxo(n,x,wts,tss,sx,cwts,mncwt,iseg,ostat,al0)
      integer n,iseg(2),al0
      double precision x(n),wts(n),tss,sx(n),cwts(n),mncwt(n),ostat
c     wts is the vector of weights for each probe
c     cwts is the cumsum of wts divided by sqrt(cwts(n))
c     mncwt is the min seg weight for a given seg length

      integer i,j, tmaxi, tmaxj, nmj, ipj, ipnmj
c      double precision xsum,sx2,x1,x2,rij,rw,tij,xvar
      double precision rn, rj, bssmx, tmax, bssij

c     max and min partial sums and their locations
      integer ipsmin, ipsmax
      double precision psmin, psmax, psrj, psrnj, rnj

      rn = cwts(n)
      ostat = -0.5
      sx(1) = x(1)*wts(1)
      psmin = sx(1)
      ipsmin = 1
      psmax = sx(1)
      ipsmax = 1
      do 20 i = 2,n
         sx(i) = sx(i-1) + x(i)*wts(i)
         if (sx(i) .lt. psmin) then 
            psmin = sx(i)
            ipsmin = i
         endif
         if (sx(i) .gt. psmax) then 
            psmax = sx(i)
            ipsmax = i
         endif
 20   continue
      psrj = abs(cwts(ipsmax) - cwts(ipsmin))
      psrnj = psrj*(rn-psrj)
      psrj = min(psrj, rn-psrj)

      bssmx = (psmax - psmin)**2/psrnj
      if (ipsmax .gt. ipsmin) then
         tmaxi = ipsmin
         tmaxj = ipsmax
      else
         tmaxi = ipsmax
         tmaxj = ipsmin
      endif

c     compute the max bss for segments of length j
      do 40 j = al0,(n-1)/2
c     sxmx is the maximum of abs partial sums x[i+1] + ... + x[i+j-1]
         do 30 i = 1,n-j
            ipj = i + j
            rj = cwts(ipj) - cwts(i)
            mncwt(j) = min(mncwt(j), rj)
            rnj = rj*(rn-rj)
            if (rnj .lt. psrnj) then
               bssij = (sx(ipj) - sx(i))**2/rnj
               if (bssij .gt. bssmx) then 
                  bssmx = bssij
                  tmaxi = i
                  tmaxj = ipj
               endif
            endif
 30      continue
         nmj = n - j
         do 35 i = 1, j
            ipnmj = i + nmj
            rj = cwts(ipnmj) - cwts(i)
            mncwt(j) = min(mncwt(j), rn-rj)
            rnj = rj*(rn-rj)
            if (rnj .lt. psrnj) then
               bssij = (sx(ipnmj) - sx(i))**2/rnj
c     absx = abs(sx(ipnmj) - sx(i))
               if (bssij .gt. bssmx) then 
                  bssmx = bssij
                  tmaxi = i
                  tmaxj = ipnmj
               endif
            endif
 35      continue
         if (mncwt(j) .ge. psrj) go to 60
 40   continue
c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         do 50 i = 1,n-j
            ipj = i + j
            rj = cwts(ipj) - cwts(i)
            rnj = rj*(rn-rj)
            if (rnj .lt. psrnj) then
               bssij = (sx(ipj) - sx(i))**2/rnj
               if (bssij .gt. bssmx) then 
                  bssmx = bssij
                  tmaxi = i
                  tmaxj = ipj
               endif
            endif
 50      continue
      endif
c     convert statistic to t^2 form
 60   if (tss .le. bssmx+0.0001) tss = bssmx + 1.0
      tmax = bssmx/((tss-bssmx)/(dfloat(n)-2.0))
      ostat = tmax
      iseg(1) = tmaxi
      iseg(2) = tmaxj

      return
      end

c     function for calculating the full max wtd t-statistic on permuted data
      double precision function wtmaxp(n,tss,px,wts,sx,cwts,mncwt,al0)
      integer n,al0
      double precision tss,px(n),wts(n),sx(n),cwts(n),mncwt(n)

      integer i, j, nmj, ipj, ipnmj
      double precision rn, rj, rnj, bssmax, bssij, psmin, psmax, psdiff,
     1     bsslim,ssq

      rn = cwts(n)
      sx(1) = px(1)*wts(1)
      ssq = wts(1)*px(1)**2
      psmin = sx(1)
      psmax = sx(1)
      do 20 i = 2,n
         sx(i) = sx(i-1) + px(i)*wts(i)
         ssq = ssq + wts(i)*px(i)**2
         psmin = min(psmin, sx(i))
         psmax = max(psmax, sx(i))
 20   continue
      psdiff = psmax - psmin
      tss = ssq - (sx(n)/rn)**2

      bssmax = 0.0
c     compute the max statistic for segments of length j
      do 50 j = al0,(n-1)/2         
         rj = mncwt(j)
c     since psdiff is the maximum sx(i+j) - sx(i) can be 
         bsslim = psdiff**2/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 70
         do 30 i = 1,n-j
            ipj = i+j
            rj = cwts(ipj) - cwts(i)
            rnj = rj*(rn-rj)
            bssij = (sx(ipj) - sx(i))**2/rnj
            if (bssij .gt. bssmax) bssmax = bssij
 30      continue
         nmj = n - j
         do 40 i = 1, j
            ipnmj = i + nmj
            rj = cwts(ipnmj) - cwts(i)
            rnj = rj*(rn-rj)
            bssij = (sx(ipnmj) - sx(i))**2/rnj
            if (bssij .gt. bssmax) bssmax = bssij
 40      continue
 50   continue

c     compute the max statistic for segments of length n/2 (if integer)
      if (n.eq.2*(n/2)) then
         j = n/2
         do 60 i = 1,n/2
            ipj = i + j
            rj = cwts(ipj) - cwts(i)
            rnj = rj*(rn-rj)
            bssij = (sx(ipj) - sx(i))**2/rnj
            if (bssij .gt. bssmax) bssmax = bssij
 60      continue
      endif
 70   if (tss .le. bssmax+0.0001) tss = bssmax + 1.0
      wtmaxp = bssmax/((tss-bssmax)/(dfloat(n)-2.0))

      return
      end

c     function for the max (over small arcs) wtd t-statistic on permuted data
      double precision function hwtmaxp(n,k,tss,px,wts,sx,cwts,mncwt,
     1     al0)
      integer n,k,al0
      double precision tss,px(n),wts(n),sx(n),cwts(n),mncwt(n)

      integer i, j, nmj, ipj, ipnmj
      double precision rn, rj, rnj, bssmax, bssij, psmin, psmax, psdiff,
     1     bsslim, ssq

      rn = cwts(n)
      sx(1) = px(1)*wts(1)
      ssq = wts(1)*px(1)**2
      psmin = sx(1)
      psmax = sx(1)
      do 20 i = 2,n
         sx(i) = sx(i-1) + px(i)*wts(i)
         ssq = ssq + wts(i)*px(i)**2
         psmin = min(psmin, sx(i))
         psmax = max(psmax, sx(i))
c         if (sx(i) .lt. psmin) psmin = sx(i)
c         if (sx(i) .gt. psmax) psmax = sx(i)
 20   continue
      psdiff = psmax - psmin
      tss = ssq - (sx(n)/rn)**2

      bssmax = 0.0d0
      do 50 j = al0,k
         rj = mncwt(j)
c     since psdiff is the maximum sx(i+j) - sx(i) can be 
         bsslim = psdiff**2/(rj*(rn-rj))
         if (bsslim .lt. bssmax) go to 60
         do 30 i = 1,n-j
            ipj = i+j
            rj = cwts(ipj) - cwts(i)
            rnj = rj*(rn-rj)
            bssij = (sx(ipj) - sx(i))**2/rnj
            if (bssij .gt. bssmax) bssmax = bssij
 30      continue
         nmj = n - j
         do 40 i = 1, j
            ipnmj = i + nmj
            rj = cwts(ipnmj) - cwts(i)
            rnj = rj*(rn-rj)
            bssij = (sx(ipnmj) - sx(i))**2/rnj
            if (bssij .gt. bssmax) bssmax = bssij
 40      continue
 50   continue
 60   if (tss .le. bssmax+0.0001d0) tss = bssmax + 1.0d0
      hwtmaxp = bssmax/((tss-bssmax)/(dfloat(n)-2.0d0))

      return
      end
