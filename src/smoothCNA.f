c     n = length(gdat) i.e. number of markers and log-ratio
c     k is the neighborhood width
c     oSD = how far the outlier is to the nearest observation
c     sSD = how close should it be moved

      subroutine smoothLR(n, gdat, sgdat, k, oSD, sSD)
      integer n, k
      double precision gdat(n), sgdat(n), oSD, sSD

      integer i, j, ilo, ihi, k1, k2, j1
      double precision mnnbd, mxnbd, xmed

c     temporary arrayd for finding median
      double precision, allocatable :: xnbhd(:)
c     order vector for calculating the median
      integer, allocatable :: loc(:)

      k2 = 2*k + 1
      allocate(xnbhd(k2), loc(k2))

      do 100 i = 1, n
c     range of the neighborhood
         ilo = max(1, i-k)
         ihi = min(n, i+k)
c     check if ith observation is an outlier
         if (i .eq. ilo) then
            mnnbd = gdat(ihi)
            mxnbd = gdat(ihi)
         else
            mnnbd = gdat(ilo)
            mxnbd = gdat(ilo)
         endif
         if ((i .le. k) .or. (i .gt. n-k)) then
            k1 = ihi - ilo + 1
            do 10 j = ilo, ihi
               if (j .ne. i) then
c     neighborhood maximum
                  if (gdat(j) .ge. mxnbd) mxnbd = gdat(j)
c     neighborhood minimum
                  if (gdat(j) .le. mnnbd) mnnbd = gdat(j)
               endif
 10         continue
         else
            k1 = k2
            do 15 j = 1, k
c     neighborhood maximum
               if (gdat(i-j) .ge. mxnbd) mxnbd = gdat(i-j)
               if (gdat(i+j) .ge. mxnbd) mxnbd = gdat(i+j)
c     neighborhood minimum
               if (gdat(i-j) .le. mnnbd) mnnbd = gdat(i-j)
               if (gdat(i+j) .le. mnnbd) mnnbd = gdat(i+j)
 15         continue
         endif
c     if it is bring it closer to the median
         if (gdat(i) .lt. mnnbd - oSD) then
            do 20 j = ilo, ihi
               j1 = j - ilo + 1
               xnbhd(j1) = gdat(j)
               loc(j1) = j1
 20         continue
            call qsort4(xnbhd, loc, 1, k1)
            if ((i .le. k) .or. (i .gt. n-k)) then
               j1 = k1/2
               if (k1 .eq. 2*j1) then
                  xmed = (xnbhd(j1) + xnbhd(j1+1))/2
               else
                  xmed = xnbhd(j1+1)
               endif
            else
               xmed = xnbhd(k+1)
            endif
            sgdat(i) = xmed - sSD
         elseif (gdat(i) .gt. mxnbd + oSD) then
            do 30 j = ilo, ihi
               j1 = j - ilo + 1
               xnbhd(j1) = gdat(j)
               loc(j1) = j1
 30         continue
            call qsort4(xnbhd, loc, 1, k1)
            if ((i .le. k) .or. (i .gt. n-k)) then
               j1 = k1/2
               if (k1 .eq. 2*j1) then
                  xmed = (xnbhd(j1) + xnbhd(j1+1))/2
               else
                  xmed = xnbhd(j1+1)
               endif
            else
               xmed = xnbhd(k+1)
            endif
            sgdat(i) = xmed + sSD
         else
            sgdat(i) = gdat(i)
         endif
 100  continue

      deallocate(xnbhd, loc)

      return
      end
