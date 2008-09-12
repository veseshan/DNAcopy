      subroutine prune(n,x,nseg,lseg,pcut,sx,ncpt,loc,loc1,pncpt)
      integer n, nseg, lseg(nseg), ncpt, loc(ncpt), loc1(2,ncpt), pncpt
      double precision x(n), pcut, sx(nseg)

      integer i, j, k, kmj
      double precision ssq, wssqk, wssq1, wssqj
      logical jleft

      double precision errssq
      external errssq

      ssq = 0.0
      do 10 i = 1,n
         ssq = ssq + x(i)**2
 10   continue
      k = 0
      do 15 i = 1,nseg
         sx(i) = 0
         do 14 j = 1,lseg(i)
            k = k + 1
            sx(i) = sx(i) + x(k)
 14      continue
 15   continue

      k = nseg - 1
      do 16 i = 1,k
         loc(i) = i
         loc1(2,i) = i
 16   continue
      wssqk = ssq - errssq(nseg,lseg,sx,k,loc)
      do 100 j = k-1, 1, -1
         kmj = k - j
         jleft = .TRUE.
         do 20 i = 1,j
            loc(i) = i
            loc1(1,i) = i
 20      continue
         wssqj = ssq - errssq(nseg,lseg,sx,j,loc)
         do 30 while(jleft) 
            call combn(j, kmj, loc, jleft)
            wssq1 = ssq - errssq(nseg,lseg,sx,j,loc)
            if (wssq1 .le. wssqj) then
               wssqj = wssq1
               do 25 i = 1,j
                  loc1(1,i) = loc(i)
 25            continue
            endif
 30      continue
         if (wssqj/wssqk .gt. 1+pcut) then
            pncpt = j+1
            do 35 i = 1,pncpt
               loc(i) = loc1(2,i)
 35         continue
            return
         else
            do 40 i = 1,j
               loc1(2,i) = loc1(1,i)
 40         continue
         endif
 100  continue
      pncpt = 0
      return
      end

      double precision function errssq(nseg,lseg,sx,k,loc)
      integer nseg, lseg(nseg),k,loc(k)
      double precision sx(nseg)

      double precision segsx
      integer segnx, i, j

      errssq = 0.0
      segsx = 0.0
      segnx = 0
      do 10 i = 1,loc(1)
         segsx = segsx + sx(i)
         segnx = segnx + lseg(i)
 10   continue
      errssq = errssq + segsx**2/dfloat(segnx)
      do 20 j = 2,k
         segsx = 0.0
         segnx = 0
         do 15 i = loc(j-1)+1,loc(j)
            segsx = segsx + sx(i)
            segnx = segnx + lseg(i)
 15      continue
         errssq = errssq + segsx**2/dfloat(segnx)
 20   continue
      segsx = 0.0
      segnx = 0
      do 25 i = loc(k)+1,nseg
         segsx = segsx + sx(i)
         segnx = segnx + lseg(i)
 25   continue
      errssq = errssq + segsx**2/dfloat(segnx)

      return
      end
c
c     This program generates Choose(n,r) combinations one at a time
c     Adapted from Algorithm AS 88  Appl. Statist. (1975) Vol.24, No. 3
c     
      subroutine combn(r, nmr, loc, rleft)
      integer r, nmr, loc(r)
      logical rleft

      integer i,j

      i = r
      do 10 while (loc(i) .eq. nmr+i)
         i = i-1
 10   continue
      loc(i) = loc(i) + 1
      do 20 j = i+1,r
         loc(j) = loc(j-1)+1
 20   continue
      if (loc(1) .eq. nmr+1) rleft = .FALSE.

      return
      end
