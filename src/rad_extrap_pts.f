c----------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for setting the mask chrbdy to a value different 
c from ex at the outermost grid points that will then be used 
c for radial extrapolation.  
c
c The extrapolation technique goes as follows:
c 1. define a range in which we look for grid points  
c to use for extrapolation. Consider the outermost  
c points in this range. 
c 2. Take one of such points, denoted by p1.  
c This will be the vertex of a cube (the other vertices  
c are other grid points). Consider the radial direction  
c connecting the point to the origin. This direction  
c intersects a face of the cube at a point p2. This will be  
c the second point used for extrapolation. The angular coords  
c of p2 are the same as p1, the radial coord must be obtained  
c and it depends on the quadrant of the 3-dimensional grid in  
c which p1 is.
c 3. To obtain the value of the function at p2, we use bilinear  
c extrapolation from the 4 vertices of the face of the cube (if  
c these are not accessible from the current process, we discard  
c the point p1 and move on to the next one).
c----------------------------------------------------------------------

        subroutine nexttobdypoints_radextrap(
     &                  chrbdy,
     &                  numbdypoints,
     &                  bdy_extrap_order,
     &                  currentres_ratio_Lhighres_Llowres,
     &                  half_steps_from_bdy_ext,
     &                  half_steps_from_bdy_int,
     &                  x,y,z,chr,L,ex,Nx,Ny,Nz,phys_bdy,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz),chrbdy2(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke
        integer ipa,jpa,kpa
        integer ipb,jpb,kpb
        integer ipc,jpc,kpc
        integer ipd,jpd,kpd

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1,rhop1,chip1,xip1
        real*8 maxxyzp1
        integer numbdypoints
        integer half_steps_from_bdy_ext
        integer half_steps_from_bdy_int
        real*8 currentres_ratio_Lhighres_Llowres
        integer bdy_extrap_order

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        numbdypoints=0

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

         do i=1,Nx
          do j=1,Ny
           do k=1,Nz
           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if

!select only points that have radius smaller than (1-5*dx/2).
!For each resolution, points between rhobdy=1 and 1-dx/2 are excised, points between rhobdy=1 and 1-3*dx/2 use forward/backward stencils (so when we cannot expect convergence at these points because the stencils used are different for different resolutions), points between rhobdy=1 and 1-5*dx/2 use points that are set by forward/backward stencils, so we cannot expect convergence. If we want to check convergence at the grid points used for extrapolation at the boundary, we need to pick points that have rho<1-5*dx/2. However, we see that we get good convergence only if we use points with rho<1-9*dx/2. This should be understood in the future. The condition rhop1>currentres_ratio_Lhighres_Llowres ensures that the first suitable point for extrapolation is not too far from the AdS boundary.
            if ((rhop1.lt.(1-half_steps_from_bdy_ext*dx/2))
     &          .and.(rhop1.gt.
     &                (1-half_steps_from_bdy_int
     &                  *currentres_ratio_Lhighres_Llowres*dx/2))
     &          .and.(chr(i,j,k).ne.ex)
     &         ) then
              chrbdy(i,j,k) =ex-1.0d0
              chrbdy2(i,j,k)=ex-1.0d0
            else
              chrbdy(i,j,k) =ex
              chrbdy2(i,j,k)=ex
            end if


            if ((i.lt.is).or.(i.gt.ie).or.
     &          (j.lt.js).or.(j.gt.je).or.
     &          (k.lt.ks).or.(k.gt.ke)) then

               chrbdy(i,j,k) =ex
               chrbdy2(i,j,k)=ex
            end if


!! eliminate troublesome points
!!removing points with abs(x)=abs(y)=abs(z): even if the grid function quasiset_tracell converges at these points, the extrapolated boundary quantity quasisettrace, obtained using these points, is not converging. It would be interesting to understand why. For now, we just eliminate these points from the set of points used for extrapolation
!            if (
!     &          (abs(abs(xp1)-abs(yp1)).lt.10.0d0**(-10)).and.
!     &          (abs(abs(xp1)-abs(zp1)).lt.10.0d0**(-10))
!     &         ) then
!                  chrbdy(i,j,k)=ex
!                  chrbdy2(i,j,k)=ex
!            end if




! eliminate troublesome points
!y=z=0 (i.e. the x axis, where chi=0 or 1) are troublesome points where the xi coordinate is not defined  
!We will fill and impose regularity at these points in post-processing in Mathematica.
            if (
     &          (abs(yp1).lt.10.0d0**(-10)).and.
     &          (abs(zp1).lt.10.0d0**(-10)) 
     &         ) then
             chrbdy(i,j,k)=ex
             chrbdy2(i,j,k)=ex
            end if

           end do
          end do
         end do


         do i=is,ie
          do j=js,je
           do k=ks,ke

           xp1=x(i)
           yp1=y(j)
           zp1=z(k)
           rhop1=sqrt(xp1**2+yp1**2+zp1**2)
           chip1=(1/PI)*acos(xp1/rhop1)
           if (zp1.lt.0) then
              xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
           else
              xip1=(1/(2*PI))*atan2(zp1,yp1)
           end if


           if (chrbdy(i,j,k).ne.ex) then


            if (bdy_extrap_order.eq.1) then


              call firstord_chrbdy_radextrap(
     &                  chrbdy,
     &                  chrbdy2,
     &                  is,ie,js,je,ks,ke,
     &                  i,j,k,
     &                  xp1,yp1,zp1,rhop1,chip1,xip1,
     &                  chr,ex,Nx,Ny,Nz)
           

            end if !closes condition on bdy_extrap_order.eq.1


          end if !closes condition on chrbdy(i,j,k).ne.ex





          if (chrbdy(i,j,k).ne.ex) then
            numbdypoints=numbdypoints+1
          end if


         end do
        end do
      end do

        return
        end


c------------------------------------------------------------------------------------------------------


c------------------------------------------------------------------------------------------------------
c in Cartesian coordinates t,x,y,z for x/y/z in [-1,1]
c
c routine for computing the Cartesian coordinates of each AdS boundary point pbdy where we want to extrapolate the value of functions 
c AND the coordinates of the outermost point pout used for radial extrapolation for each of the boundary points
c-------------------------------------------------------------------------------------------------------------------------

        subroutine xyz_bdy_out_radextrap(
     &                  xpbdy,ypbdy,zpbdy,
     &                  xpout,ypout,zpout,
     &                  chrbdy,
     &                  numbdypoints,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,ghost_width)

!----------------------------------------------------------------------

        implicit none

        integer Nx,Ny,Nz
        integer phys_bdy(6),ghost_width(6)
        real*8 L
        real*8 x(Nx),y(Ny),z(Nz),dt,chr(Nx,Ny,Nz),ex

        real*8 chrbdy(Nx,Ny,Nz)

        integer i,j,k,is,ie,js,je,ks,ke

        integer ix1,ix2,ix3,jy1,jy2,jy3,kz1,kz2,kz3
        integer ix,jy,kz
        integer lind

        real*8 x0,y0,z0,rho0,q

        real*8 dx,dy,dz
        real*8 xp1,yp1,zp1
        real*8 rhop1,chip1,xip1
        real*8 rhopbdy,chipbdy,xipbdy
        real*8 maxxyzp1
        integer numbdypoints
        real*8 xpbdy(numbdypoints)
        real*8 ypbdy(numbdypoints)
        real*8 zpbdy(numbdypoints)
        real*8 xpout(numbdypoints)
        real*8 ypout(numbdypoints)
        real*8 zpout(numbdypoints)

        real*8 PI
        parameter (PI=3.141592653589793d0)

!----------------------------------------------------------------------

        dx=x(2)-x(1)
        dy=y(2)-y(1)
        dz=z(2)-z(1)

        ! set index bounds for main loop
        is=2
        ie=Nx-1
        js=2
        je=Ny-1
        ks=2
        ke=Nz-1

        ! adjust index bounds to compensate for ghost_width
        if (ghost_width(1).gt.0) is=is+ghost_width(1)-1
        if (ghost_width(2).gt.0) ie=ie-(ghost_width(2)-1)
        if (ghost_width(3).gt.0) js=js+ghost_width(3)-1
        if (ghost_width(4).gt.0) je=je-(ghost_width(4)-1)
        if (ghost_width(5).gt.0) ks=ks+ghost_width(5)-1
        if (ghost_width(6).gt.0) ke=ke-(ghost_width(6)-1)

        lind=0
        do i=is,ie
         do j=js,je
          do k=ks,ke

            if (chrbdy(i,j,k).ne.ex) then
              lind=lind+1

              xp1=x(i)
              yp1=y(j)
              zp1=z(k)
              rhop1=sqrt(xp1**2+yp1**2+zp1**2)
              chip1=(1/PI)*acos(xp1/rhop1)
              if (zp1.lt.0) then
                xip1=(1/(2*PI))*(atan2(zp1,yp1)+2*PI)
              else
                xip1=(1/(2*PI))*atan2(zp1,yp1)
              end if

              xpout(lind)=xp1
              ypout(lind)=yp1
              zpout(lind)=zp1

              rhopbdy=1
              chipbdy=chip1
              xipbdy=xip1

              xpbdy(lind)=rhopbdy*cos(PI*chipbdy)
              ypbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *cos(2*PI*xipbdy)
              zpbdy(lind)=rhopbdy*sin(PI*chipbdy)
     &              *sin(2*PI*xipbdy)


            end if !closes condition on (chrbdy(i,j,k).ne.ex)
          end do
         end do
        end do

        return
        end
c--------------------------------------------------------------------------------------