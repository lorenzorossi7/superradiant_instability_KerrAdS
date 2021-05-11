c----------------------------------------------------------------------
c routine to select chrbdy mask, which is not ex only at the outermost points
c that we use for first order radial extrapolation
c
c If we use derivatives to define near boundary quantities, we will only define 
c them at points between is and ie (js and je, ks and ke).  
c Therefore, for extrapolation, we can only select near boundary points whose  
c neighbors used for extrapolation are within that range
c We also need to make sure that those neighbours are not excised.
c The condition (chrbdy2(i+1,j,k).ne.ex) makes sure that (i,j,k) is the 
c outmost point satisfying the conditions of the previous for-loop, which 
c sets chrbdy2 as well as chrbdy. 
c In other words, if there's an outer point w.r.t. (i,j,k) that satisfies 
c those conditions, then we don't want to use (i,j,k) for extrapolation, 
c but we will use that other point. 
c----------------------------------------------------------------------

        subroutine firstord_chrbdy_radextrap(
     &                  chrbdy,
     &                  chrbdy2,
     &                  is,ie,js,je,ks,ke,
     &                  i,j,k,
     &                  xp1,yp1,zp1,rhop1,chip1,xip1,
     &                  chr,ex,Nx,Ny,Nz)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

              if ((abs(xp1).gt.abs(yp1)).and.
     &            (abs(xp1).gt.abs(zp1))) then !(i.e., |xp1|>|yp1|,|zp1|, so xp1 cannot be 0)
            
                if (xp1.gt.0) then !(i.e., we are in the upper part of 
                                    !the spatial grid, called "a" sector)

                  ipa=i-1
                  ipb=i-1
                  ipc=i-1
                  ipd=i-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i+1).le.Nx) then 
                    if (chrbdy2(i+1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                else !(i.e., xp1<0, i.e., we are in the lower part of 
                     !the spatial grid, called "b" sector)

                  ipa=i+1
                  ipb=i+1
                  ipc=i+1
                  ipd=i+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((i-1).ge.1) then 
                    if (chrbdy2(i-1,j,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(yp1).gt.abs(zp1)) then !(i.e., |yp1|>|zp1|, so yp1 cannot be 0)
                    if (yp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIIa)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |zp1|>=|yp1|, it's possible that zp1=yp1=0)

                    if (zp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IVa)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(yp1).gt.abs(zp1))

                end if !closes condition on (xp1.gt.0)











              else if ((abs(yp1).gt.abs(zp1)).and.
     &                 (abs(yp1).gt.abs(xp1))) then !(i.e., |yp1|>|zp1|,|xp1|, so yp1 cannot be 0)


                if (yp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IVa or Q.IVb)

                  jpa=j-1
                  jpb=j-1
                  jpc=j-1
                  jpd=j-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j+1).le.Ny) then 
                    if (chrbdy2(i,j+1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IVa or Q.IVb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1| - xp1=zp1=0 is not an issue)

                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IVa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (zp1.gt.0) then !(i.e., quadrant Ia)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (zp1.gt.0) then !(i.e., quadrant Ib)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IVb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                else !(i.e., yp1<0, so we are in either Q.IIa or Q.IIb or Q.IIIa or Q.IIIb)

                  jpa=j+1
                  jpb=j+1
                  jpc=j+1
                  jpd=j+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((j-1).ge.1) then 
                    if (chrbdy2(i,j-1,k).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(zp1).gt.abs(xp1)) then !(i.e., |zp1|>|xp1|, so zp1 cannot be 0)
                    if (zp1.gt.0) then !(i.e., either quadrant IIa or IIb)

                      kpa=k
                      kpb=k
                      kpc=k-1
                      kpd=k-1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., zp1<=0, so Q.IIIa or IIIb)

                      kpa=k
                      kpb=k
                      kpc=k+1
                      kpd=k+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |xp1|>=|zp1|, it's possible that xp1=zp1=0)

                    if (xp1.gt.0) then !(i.e., either quadrant IIa or IIIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (zp1.gt.0) then !(i.e., quadrant IIa)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIa)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIb or Q.IIIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (zp1.gt.0) then !(i.e., quadrant IIb)

                        kpa=k
                        kpb=k-1
                        kpc=k-1
                        kpd=k


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., zp1<=0: quadrant IIIb)

                        kpa=k
                        kpb=k+1
                        kpc=k+1
                        kpd=k

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(zp1).gt.abs(xp1))

                end if !closes condition on (yp1.gt.0)











              else !(i.e., |zp1|>=|xp1|,|yp1|, so zp1 cannot be 0)

                if (zp1.gt.0) then !(i.e., we are in either Q.Ia or Q.Ib or Q.IIa or Q.IIb)

                  kpa=k-1
                  kpb=k-1
                  kpc=k-1
                  kpd=k-1

                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k+1).le.ke) then 
                    if (chrbdy2(i,j,k+1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant Ia or IIa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (yp1.gt.0) then !(i.e., quadrant Ia)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.Ib or IIb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (yp1.gt.0) then !(i.e., quadrant Ib)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1| - yp1=xp1=0 is not an issue)

                    if (yp1.gt.0) then !(i.e., either quadrant Ia or Ib)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (xp1.gt.0) then !(i.e., quadrant Ia)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant Ib)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIa or IIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k-1).lt.ks)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))


                else !(i.e., zp1<0, so we are in either Q.IIIa or Q.IIIb or Q.IVa or Q.IVb)

                  kpa=k+1
                  kpb=k+1
                  kpc=k+1
                  kpd=k+1


                  !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation 
                  !if it is not the outermost point in the range in which we look for candidates
                  if ((k-1).ge.ks) then 
                    if (chrbdy2(i,j,k-1).ne.ex) then
                      chrbdy(i,j,k)=ex
                    end if
                  end if

                  if (abs(xp1).gt.abs(yp1)) then !(i.e., |xp1|>|yp1|, so xp1 cannot be 0)
                    if (xp1.gt.0) then !(i.e., either quadrant IIIa or IVa)

                      ipa=i
                      ipb=i
                      ipc=i-1
                      ipd=i-1

                      if (yp1.gt.0) then !(i.e., quadrant IVa)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIa)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., xp1<=0, so Q.IIIb or IVb)

                      ipa=i
                      ipb=i
                      ipc=i+1
                      ipd=i+1

                      if (yp1.gt.0) then !(i.e., quadrant IVb)

                        jpa=j
                        jpb=j-1
                        jpc=j-1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., yp1<=0: quadrant IIIb)

                        jpa=j
                        jpb=j+1
                        jpc=j+1
                        jpd=j

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  else !(i.e., |yp1|>=|xp1|, it's possible that yp1=xp1=0)

                    if (yp1.gt.0) then !(i.e., either quadrant IVa or IVb)

                      jpa=j
                      jpb=j
                      jpc=j-1
                      jpd=j-1

                      if (xp1.gt.0) then !(i.e., quadrant IVa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IVb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j-1).lt.js).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    else !(i.e., yp1<=0, so Q.IIIa or IIIb)

                      jpa=j
                      jpb=j
                      jpc=j+1
                      jpd=j+1

                      if (xp1.gt.0) then !(i.e., quadrant IIIa)

                        ipa=i
                        ipb=i-1
                        ipc=i-1
                        ipd=i


                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i-1).lt.is).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if


                      else !(i.e., xp1<=0: quadrant IIIb)

                        ipa=i
                        ipb=i+1
                        ipc=i+1
                        ipd=i

                        !eliminates the grid point p1 from the list of candidate points to use for radial extrapolation
                        !if the grid points needed to extrapolate the value of the function at the second point used for radial extrapolation p2 (p2 is not a grid point in general) are not available 
                        if (((i+1).gt.ie).or.
     &                      ((j+1).gt.je).or.
     &                      ((k+1).gt.ke)) then 
                          chrbdy(i,j,k)=ex
                        else if ((chr(ipa,jpa,kpa).eq.ex).or.
     &                           (chr(ipb,jpb,kpb).eq.ex).or. 
     &                           (chr(ipc,jpc,kpc).eq.ex).or. 
     &                           (chr(ipd,jpd,kpd).eq.ex) 
     &                                          ) then
                          chrbdy(i,j,k)=ex
                        end if

                      end if

                    end if

                  end if !closes condition on (abs(xp1).gt.abs(yp1))

                end if !closes condition on (zp1.gt.0)






              end if !closes condition on ((abs(xp1).gt.abs(yp1)).and.(abs(xp1).gt.abs(zp1)))

        return
        end


c------------------------------------------------------------------------------------------------------