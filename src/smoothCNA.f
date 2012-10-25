c     n = length(gdat) i.e. number of markers and log-ratio
c     k is the neighborhood width
c     oSD = how far the outlier is to the nearest observation
c     sSD = how close should it be moved
c     nchr = number of chromosomes
c     cfrq = number of probes in respective chromosomes
      subroutine smoothLR(n, gdat, nchr, cfrq, sgdat, k, oSD, sSD)
      integer n, nchr, cfrq(nchr), k
      double precision gdat(n), sgdat(n), oSD, sSD

      integer i, j, ilo, ihi, k1, k2, j1, ic, cilo, cihi
      double precision mnnbd, mxnbd, distij, xmed

c     temporary array for finding median
      double precision, allocatable :: xnbhd(:)

      k2 = 2*k + 1
      allocate(xnbhd(k2))

c     initial values for start and end of chromosomes
      cilo = 1
      cihi = 0
c     loop over chromomsomes
      do 100 ic = 1, nchr
c     end of the current chromosome
         cihi = cihi + cfrq(ic)
         do 50 i = cilo, cihi
c     range of the neighborhood
            ilo = max(cilo, i-k)
            ihi = min(cihi, i+k)
c     check if ith observation is an outlier
c     initialize the distances to be large
            mxnbd = 100*oSD
            mnnbd = 100*oSD
            do 10 j = ilo, ihi
               if (j .ne. i) then
c     calculate distance from between ith and jth obsn
                  distij = gdat(i) - gdat(j)
c     if distance is less than oSD no smoothing necessary
                  if (abs(distij) .le. oSD) then
                     sgdat(i) = gdat(i)
                     go to 50
c     otherwise calculate distances from above and below
                  else
c     mxnbd is distance from above
                     if (distij .lt. mxnbd) mxnbd = distij
c     mnnbd is distance from below
                     if (-distij .lt. mnnbd) mnnbd = -distij
                  endif
               endif
 10         continue
c     if all the points in the nbhd are above than mxnbd will be negative
c     and mnnbd will be greater than oSD. Vice versa if all points below
c     
c     If both are negative then the ith point is singleton in the middle 
c     but distance oSD away from all points in the nbhd. No smoothing done.
            if ((mxnbd .le. 0) .and. (mnnbd .le. 0)) then
               sgdat(i) = gdat(i)
               go to 50
            else
c     calculate the median of the nbhd
c     number of points in the nbhd
               k1 = ihi - ilo + 1
c     get the data into temporary array
               do 20 j = ilo, ihi
                  xnbhd(j-ilo+1) = gdat(j)
 20            continue
c     sort the data
               call qsort3(xnbhd, 1, k1)
c     median is the midpoint if n is odd and average of the two if even
               j1 = k1/2
               if (k1 .eq. 2*j1) then
                  xmed = (xnbhd(j1) + xnbhd(j1+1))/2
               else
                  xmed = xnbhd(j1+1)
               endif
c     if point is above the nbhd bring it down
               if (mxnbd .gt. 0) sgdat(i) = xmed + sSD
c     if point is below the nbhd bring it up
               if (mnnbd .gt. 0) sgdat(i) = xmed - sSD
            endif
 50      continue
c     beginning of next chromosome
         cilo = cilo + cfrq(ic)
 100  continue
      deallocate(xnbhd)

      return
      end
