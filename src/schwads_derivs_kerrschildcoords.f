c----------------------------------------------------------------------
c calculates Kerr-schwads metric and its first and second derivatives in the version of 
c the Kerr-Schild coordinates (horizon-penetrating) that is non-rotating at the boundary.
c
c r0 is the radius parameter, i.e. r0=2*M0, where M0 is the BHmass 
c (r0 has no physical meaning, it is NOT the horizon radius)
c----------------------------------------------------------------------
        subroutine schwads_derivs_kerrschildcoords(
     &                  gschwads_ll,gschwads_uu,gschwads_ll_x,
     &                  gschwads_uu_x,gschwads_ll_xx,
     &                  Hschwads_l,
     &                  gammaschwads_ull,
     &                  phi1schwads,
     &                  phi1schwads_x,
     &                  x,y,z,dt,chr,L,ex,Nx,Ny,Nz,i,j,k,
     &                  ief_bh_r0)

        implicit none

        integer Nx,Ny,Nz
        integer i,j,k
        real*8  ief_bh_r0,a_rot,M0,M0_min

        real*8 chr(Nx,Ny,Nz),ex
        real*8 x(Nx),y(Ny),z(Nz),dt,L

        integer a,b,c,d,e,f,g,h
        real*8 dx,dy,dz
        real*8 x0,y0,z0
        real*8 rho0,theta0,phi0
        real*8 f0

        real*8 PI
        parameter (PI=3.141592653589793d0)

        !--------------------------------------------------------------
        ! variables for tensor manipulations
        !(indices are t,x,w,y,z)
        !--------------------------------------------------------------

        real*8 gschwads_ll(4,4),gschwads_uu(4,4)
        real*8 gschwads_ll_x(4,4,4),gschwads_uu_x(4,4,4)
        real*8 gschwads_ll_xx(4,4,4,4)
        real*8 Hschwads_l(4)
        real*8 phi10_x(4),phi10_xx(4,4)
        real*8 gschwads_ll_sph(4,4),gschwads_uu_sph(4,4)
        real*8 gschwads_ll_sph_x(4,4,4),gschwads_uu_sph_x(4,4,4)
        real*8 gschwads_ll_sph_xx(4,4,4,4)
        real*8 dxsph_dxcar(4,4)
        real*8 d2xsph_dxcardxcar(4,4,4)
        real*8 d3xsph_dxcardxcardxcar(4,4,4,4)
        real*8 gammaschwads_ull(4,4,4)
        real*8 boxschwadsx_u(4)
        real*8 phi1schwads, phi1schwads_x(4)
        real*8 phi1schwads_xx(4,4)
        real*8 gammaschwads_ull_x(4,4,4,4)
        real*8 riemannschwads_ulll(4,4,4,4)
        real*8 riccischwads_ll(4,4),riccischwads_lu(4,4)
        real*8 riccischwads
        real*8 einsteinschwads_ll(4,4),setschwads_ll(4,4)
        real*8 grad_phi1schwads_sq

        real*8 g0_tt_schwads_sph_t
        real*8 g0_tt_schwads_sph_tt
        real*8 g0_tt_schwads_sph_trho
        real*8 g0_tt_schwads_sph_ttheta
        real*8 g0_tt_schwads_sph_tphi
        real*8 g0_tt_schwads_sph_rho
        real*8 g0_tt_schwads_sph_rhorho
        real*8 g0_tt_schwads_sph_rhotheta
        real*8 g0_tt_schwads_sph_theta
        real*8 g0_tt_schwads_sph_thetatheta
        real*8 g0_tt_schwads_sph_phi
        real*8 g0_tt_schwads_sph_rhophi
        real*8 g0_tt_schwads_sph_thetaphi
        real*8 g0_tt_schwads_sph_phiphi
        real*8 g0_trho_schwads_sph_t
        real*8 g0_trho_schwads_sph_tt
        real*8 g0_trho_schwads_sph_trho
        real*8 g0_trho_schwads_sph_ttheta
        real*8 g0_trho_schwads_sph_tphi
        real*8 g0_trho_schwads_sph_rho
        real*8 g0_trho_schwads_sph_rhorho
        real*8 g0_trho_schwads_sph_rhotheta
        real*8 g0_trho_schwads_sph_theta
        real*8 g0_trho_schwads_sph_thetatheta
        real*8 g0_trho_schwads_sph_phi
        real*8 g0_trho_schwads_sph_rhophi
        real*8 g0_trho_schwads_sph_thetaphi
        real*8 g0_trho_schwads_sph_phiphi
        real*8 g0_ttheta_schwads_sph_t
        real*8 g0_ttheta_schwads_sph_tt
        real*8 g0_ttheta_schwads_sph_trho
        real*8 g0_ttheta_schwads_sph_ttheta
        real*8 g0_ttheta_schwads_sph_tphi
        real*8 g0_ttheta_schwads_sph_rho
        real*8 g0_ttheta_schwads_sph_rhorho
        real*8 g0_ttheta_schwads_sph_rhotheta
        real*8 g0_ttheta_schwads_sph_theta
        real*8 g0_ttheta_schwads_sph_thetatheta
        real*8 g0_ttheta_schwads_sph_phi
        real*8 g0_ttheta_schwads_sph_rhophi
        real*8 g0_ttheta_schwads_sph_thetaphi
        real*8 g0_ttheta_schwads_sph_phiphi
        real*8 g0_tphi_schwads_sph_t
        real*8 g0_tphi_schwads_sph_tt
        real*8 g0_tphi_schwads_sph_trho
        real*8 g0_tphi_schwads_sph_ttheta
        real*8 g0_tphi_schwads_sph_tphi
        real*8 g0_tphi_schwads_sph_rho
        real*8 g0_tphi_schwads_sph_rhorho
        real*8 g0_tphi_schwads_sph_rhotheta
        real*8 g0_tphi_schwads_sph_theta
        real*8 g0_tphi_schwads_sph_thetatheta
        real*8 g0_tphi_schwads_sph_phi
        real*8 g0_tphi_schwads_sph_rhophi
        real*8 g0_tphi_schwads_sph_thetaphi
        real*8 g0_tphi_schwads_sph_phiphi
        real*8 g0_rhorho_schwads_sph_t
        real*8 g0_rhorho_schwads_sph_tt
        real*8 g0_rhorho_schwads_sph_trho
        real*8 g0_rhorho_schwads_sph_ttheta
        real*8 g0_rhorho_schwads_sph_tphi
        real*8 g0_rhorho_schwads_sph_rho
        real*8 g0_rhorho_schwads_sph_rhorho
        real*8 g0_rhorho_schwads_sph_rhotheta
        real*8 g0_rhorho_schwads_sph_theta
        real*8 g0_rhorho_schwads_sph_thetatheta
        real*8 g0_rhorho_schwads_sph_phi
        real*8 g0_rhorho_schwads_sph_rhophi
        real*8 g0_rhorho_schwads_sph_thetaphi
        real*8 g0_rhorho_schwads_sph_phiphi
        real*8 g0_rhotheta_schwads_sph_t
        real*8 g0_rhotheta_schwads_sph_tt
        real*8 g0_rhotheta_schwads_sph_trho
        real*8 g0_rhotheta_schwads_sph_ttheta
        real*8 g0_rhotheta_schwads_sph_tphi
        real*8 g0_rhotheta_schwads_sph_rho
        real*8 g0_rhotheta_schwads_sph_rhorho
        real*8 g0_rhotheta_schwads_sph_rhotheta
        real*8 g0_rhotheta_schwads_sph_theta
        real*8 g0_rhotheta_schwads_sph_thetatheta
        real*8 g0_rhotheta_schwads_sph_phi
        real*8 g0_rhotheta_schwads_sph_rhophi
        real*8 g0_rhotheta_schwads_sph_thetaphi
        real*8 g0_rhotheta_schwads_sph_phiphi
        real*8 g0_rhophi_schwads_sph_t
        real*8 g0_rhophi_schwads_sph_tt
        real*8 g0_rhophi_schwads_sph_trho
        real*8 g0_rhophi_schwads_sph_ttheta
        real*8 g0_rhophi_schwads_sph_tphi
        real*8 g0_rhophi_schwads_sph_rho
        real*8 g0_rhophi_schwads_sph_rhorho
        real*8 g0_rhophi_schwads_sph_rhotheta
        real*8 g0_rhophi_schwads_sph_theta
        real*8 g0_rhophi_schwads_sph_thetatheta
        real*8 g0_rhophi_schwads_sph_phi
        real*8 g0_rhophi_schwads_sph_rhophi
        real*8 g0_rhophi_schwads_sph_thetaphi
        real*8 g0_rhophi_schwads_sph_phiphi
        real*8 g0_thetatheta_schwads_sph_t
        real*8 g0_thetatheta_schwads_sph_tt
        real*8 g0_thetatheta_schwads_sph_trho
        real*8 g0_thetatheta_schwads_sph_ttheta
        real*8 g0_thetatheta_schwads_sph_tphi
        real*8 g0_thetatheta_schwads_sph_rho
        real*8 g0_thetatheta_schwads_sph_rhorho
        real*8 g0_thetatheta_schwads_sph_rhotheta
        real*8 g0_thetatheta_schwads_sph_theta
        real*8 g0_thetatheta_schwads_sph_thetatheta
        real*8 g0_thetatheta_schwads_sph_phi
        real*8 g0_thetatheta_schwads_sph_rhophi
        real*8 g0_thetatheta_schwads_sph_thetaphi
        real*8 g0_thetatheta_schwads_sph_phiphi
        real*8 g0_thetaphi_schwads_sph_t
        real*8 g0_thetaphi_schwads_sph_tt
        real*8 g0_thetaphi_schwads_sph_trho
        real*8 g0_thetaphi_schwads_sph_ttheta
        real*8 g0_thetaphi_schwads_sph_tphi
        real*8 g0_thetaphi_schwads_sph_rho
        real*8 g0_thetaphi_schwads_sph_rhorho
        real*8 g0_thetaphi_schwads_sph_rhotheta
        real*8 g0_thetaphi_schwads_sph_theta
        real*8 g0_thetaphi_schwads_sph_thetatheta
        real*8 g0_thetaphi_schwads_sph_phi
        real*8 g0_thetaphi_schwads_sph_rhophi
        real*8 g0_thetaphi_schwads_sph_thetaphi
        real*8 g0_thetaphi_schwads_sph_phiphi
        real*8 g0_phiphi_schwads_sph_t
        real*8 g0_phiphi_schwads_sph_tt
        real*8 g0_phiphi_schwads_sph_trho
        real*8 g0_phiphi_schwads_sph_ttheta
        real*8 g0_phiphi_schwads_sph_tphi
        real*8 g0_phiphi_schwads_sph_rho
        real*8 g0_phiphi_schwads_sph_rhorho
        real*8 g0_phiphi_schwads_sph_rhotheta
        real*8 g0_phiphi_schwads_sph_theta
        real*8 g0_phiphi_schwads_sph_thetatheta
        real*8 g0_phiphi_schwads_sph_phi
        real*8 g0_phiphi_schwads_sph_rhophi
        real*8 g0_phiphi_schwads_sph_thetaphi
        real*8 g0_phiphi_schwads_sph_phiphi

        real*8 g0_tt_schwads_sph0,g0_rhorho_schwads_sph0
        real*8 g0_trho_schwads_sph0,g0_ttheta_schwads_sph0
        real*8 g0_tphi_schwads_sph0
        real*8 g0_rhotheta_schwads_sph0
        real*8 g0_thetatheta_schwads_sph0
        real*8 g0_phiphi_schwads_sph0
        real*8 g0_rhophi_schwads_sph0,g0_thetaphi_schwads_sph0

        real*8 detg0_schwads_sph0
        real*8 g0u_tt_schwads_sph0,g0u_rhorho_schwads_sph0
        real*8 g0u_trho_schwads_sph0,g0u_ttheta_schwads_sph0
        real*8 g0u_tphi_schwads_sph0
        real*8 g0u_rhotheta_schwads_sph0
        real*8 g0u_thetatheta_schwads_sph0
        real*8 g0u_phiphi_schwads_sph0
        real*8 g0u_rhophi_schwads_sph0
        real*8 g0u_thetaphi_schwads_sph0

        real*8 detg0_schwads0

!!!!!!!!!!DEBUG DERIVATIVE STENCILS!!!!!!!!!!!
        real*8 testf1(Nx,Ny,Nz),testf2(Nx,Ny,Nz),testf3(Nx,Ny,Nz)
        real*8 testf1_t,testf1_x,testf1_y,testf1_z
        real*8 testf2_t,testf2_x,testf2_y,testf2_z
        real*8 testf3_t,testf3_x,testf3_y,testf3_z
        real*8 testf1_tt,testf1_tx,testf1_ty
        real*8 testf1_xx,testf1_xy,testf1_yy
        real*8 testf1_tz,testf1_xz,testf1_yz,testf1_zz
        real*8 testf2_tt,testf2_tx,testf2_ty
        real*8 testf2_xx,testf2_xy,testf2_yy
        real*8 testf2_tz,testf2_xz,testf2_yz,testf2_zz
        real*8 testf3_tt,testf3_tx,testf3_ty
        real*8 testf3_xx,testf3_xy,testf3_yy
        real*8 testf3_tz,testf3_xz,testf3_yz,testf3_zz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------------
        
        dx=(x(2)-x(1))
        dy=(y(2)-y(1))
        dz=(z(2)-z(1))

        x0=x(i)
        y0=y(j)
        z0=z(k)
        rho0=sqrt(x0**2+y0**2+z0**2)
        if (rho0.ne.0.0d0) then
          theta0=acos(x0/rho0)
        end if
        if ((y0.ne.0.0d0).or.(z0.ne.0.0d0)) then
          phi0=atan2(z0,y0)
          if (phi0.lt.0) phi0=phi0+2*PI
        end if

      ! Black hole mass
        M0=ief_bh_r0/2


        g0_tt_schwads_sph0 =
     -   M0/rho0 - M0*rho0 - 
     -  (4*rho0**2 + (-1 + rho0**2)**2)/
     -   ((-1 + rho0)**2*(1 + rho0)**2)
        g0_trho_schwads_sph0 =
     -    (2*M0*(1 - 2*rho0**2 + rho0**4))/
     -  (rho0*Sqrt((-1 + rho0**2)**2)*(1 + rho0**2))
        g0_ttheta_schwads_sph0 =
     -   0
        g0_tphi_schwads_sph0 =
     -   0
        g0_rhorho_schwads_sph0 =
     -   (-4*M0*(-1 + rho0**2))/
     -   (rho0*(1 + rho0**2)**2) + 
     -  (4*(1 + rho0**2)**2)/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2))
        g0_rhotheta_schwads_sph0 =
     -   0
        g0_rhophi_schwads_sph0 =
     -   0
        g0_thetatheta_schwads_sph0 =
     -    (4*rho0**2)/(-1 + rho0**2)**2
        g0_thetaphi_schwads_sph0 =
     -   0
        g0_phiphi_schwads_sph0 =
     -   (4*rho0**2*Sin(theta0)**2)/(-1 + rho0**2)**2

             ! give values to the metric inverse
        call calc_g0uu(g0_tt_schwads_sph0,
     &         g0_trho_schwads_sph0,
     &         g0_ttheta_schwads_sph0,
     &         g0_tphi_schwads_sph0,
     &         g0_rhorho_schwads_sph0,
     &         g0_rhotheta_schwads_sph0,
     &         g0_rhophi_schwads_sph0,
     &         g0_thetatheta_schwads_sph0,
     &         g0_thetaphi_schwads_sph0,
     &         g0_phiphi_schwads_sph0,
     &         g0u_tt_schwads_sph0,
     &         g0u_trho_schwads_sph0,
     &         g0u_ttheta_schwads_sph0,
     &         g0u_tphi_schwads_sph0,
     &         g0u_rhorho_schwads_sph0,
     &         g0u_rhotheta_schwads_sph0,
     &         g0u_rhophi_schwads_sph0,
     &         g0u_thetatheta_schwads_sph0,
     &         g0u_thetaphi_schwads_sph0,
     &         g0u_phiphi_schwads_sph0,
     &         detg0_schwads_sph0)


        g0_tt_schwads_sph_t  =
     -   0
        g0_tt_schwads_sph_rho  =
     -   -M0 - M0/rho0**2 - 
     -  (8*rho0 + 4*rho0*(-1 + rho0**2))/
     -   ((-1 + rho0)**2*(1 + rho0)**2) + 
     -  (2*(4*rho0**2 + (-1 + rho0**2)**2))/
     -   ((-1 + rho0)**2*(1 + rho0)**3) + 
     -  (2*(4*rho0**2 + (-1 + rho0**2)**2))/
     -   ((-1 + rho0)**3*(1 + rho0)**2)
        g0_tt_schwads_sph_theta  =
     -   0
        g0_tt_schwads_sph_phi  =
     -   0
        g0_tt_schwads_sph_tt =
     -   0
        g0_tt_schwads_sph_trho =
     -   0
        g0_tt_schwads_sph_ttheta =
     -   0
        g0_tt_schwads_sph_tphi =
     -   0
        g0_tt_schwads_sph_rhorho =
     -    (2*M0)/rho0**3 - 
     -  (8 + 8*rho0**2 + 4*(-1 + rho0**2))/
     -   ((-1 + rho0)**2*(1 + rho0)**2) + 
     -  (4*(8*rho0 + 4*rho0*(-1 + rho0**2)))/
     -   ((-1 + rho0)**2*(1 + rho0)**3) + 
     -  (4*(8*rho0 + 4*rho0*(-1 + rho0**2)))/
     -   ((-1 + rho0)**3*(1 + rho0)**2) - 
     -  (6*(4*rho0**2 + (-1 + rho0**2)**2))/
     -   ((-1 + rho0)**2*(1 + rho0)**4) - 
     -  (8*(4*rho0**2 + (-1 + rho0**2)**2))/
     -   ((-1 + rho0)**3*(1 + rho0)**3) - 
     -  (6*(4*rho0**2 + (-1 + rho0**2)**2))/
     -   ((-1 + rho0)**4*(1 + rho0)**2)
        g0_tt_schwads_sph_rhotheta =
     -   0
        g0_tt_schwads_sph_rhophi =
     -   0
        g0_tt_schwads_sph_thetatheta =
     -   0
        g0_tt_schwads_sph_thetaphi =
     -   0
        g0_tt_schwads_sph_phiphi =
     -   0

        g0_trho_schwads_sph_t  =
     -   0
        g0_trho_schwads_sph_rho  =
     -   (2*M0*(-4*rho0 + 4*rho0**3))/
     -   (rho0*Sqrt((-1 + rho0**2)**2)*(1 + rho0**2))
     -   - (4*M0*(1 - 2*rho0**2 + rho0**4))/
     -   (Sqrt((-1 + rho0**2)**2)*(1 + rho0**2)**2) - 
     -  (4*M0*(-1 + rho0**2)*
     -     (1 - 2*rho0**2 + rho0**4))/
     -   (((-1 + rho0**2)**2)**1.5*(1 + rho0**2)) - 
     -  (2*M0*(1 - 2*rho0**2 + rho0**4))/
     -   (rho0**2*Sqrt((-1 + rho0**2)**2)*
     -     (1 + rho0**2))
        g0_trho_schwads_sph_theta  =
     -   0
        g0_trho_schwads_sph_phi  =
     -   0
        g0_trho_schwads_sph_tt =
     -   0
        g0_trho_schwads_sph_trho =
     -   0
        g0_trho_schwads_sph_ttheta =
     -   0
        g0_trho_schwads_sph_tphi =
     -   0
        g0_trho_schwads_sph_rhorho =
     -   (2*M0*(-4 + 12*rho0**2))/
     -   (rho0*Sqrt((-1 + rho0**2)**2)*(1 + rho0**2))
     -   - (8*M0*(-4*rho0 + 4*rho0**3))/
     -   (Sqrt((-1 + rho0**2)**2)*(1 + rho0**2)**2) - 
     -  (8*M0*(-1 + rho0**2)*(-4*rho0 + 4*rho0**3))/
     -   (((-1 + rho0**2)**2)**1.5*(1 + rho0**2)) - 
     -  (4*M0*(-4*rho0 + 4*rho0**3))/
     -   (rho0**2*Sqrt((-1 + rho0**2)**2)*
     -     (1 + rho0**2)) + 
     -  (16*M0*rho0*(1 - 2*rho0**2 + rho0**4))/
     -   (Sqrt((-1 + rho0**2)**2)*(1 + rho0**2)**3) + 
     -  (16*M0*rho0*(-1 + rho0**2)*
     -     (1 - 2*rho0**2 + rho0**4))/
     -   (((-1 + rho0**2)**2)**1.5*(1 + rho0**2)**2) + 
     -  (4*M0*(1 - 2*rho0**2 + rho0**4))/
     -   (rho0*Sqrt((-1 + rho0**2)**2)*
     -     (1 + rho0**2)**2) + 
     -  (16*M0*rho0*(1 - 2*rho0**2 + rho0**4))/
     -   (((-1 + rho0**2)**2)**1.5*(1 + rho0**2)) + 
     -  (4*M0*(-1 + rho0**2)*
     -     (1 - 2*rho0**2 + rho0**4))/
     -   (rho0*((-1 + rho0**2)**2)**1.5*(1 + rho0**2))
     -   + (4*M0*(1 - 2*rho0**2 + rho0**4))/
     -   (rho0**3*Sqrt((-1 + rho0**2)**2)*
     -     (1 + rho0**2))
        g0_trho_schwads_sph_rhotheta =
     -   0
        g0_trho_schwads_sph_rhophi =
     -   0
        g0_trho_schwads_sph_thetatheta =
     -   0
        g0_trho_schwads_sph_thetaphi =
     -   0
        g0_trho_schwads_sph_phiphi =
     -   0

        g0_ttheta_schwads_sph_t  =
     -   0
        g0_ttheta_schwads_sph_rho  =
     -   0
        g0_ttheta_schwads_sph_theta  =
     -   0
        g0_ttheta_schwads_sph_phi  =
     -   0
        g0_ttheta_schwads_sph_tt =
     -   0
        g0_ttheta_schwads_sph_trho =
     -   0
        g0_ttheta_schwads_sph_ttheta =
     -   0
        g0_ttheta_schwads_sph_tphi =
     -   0
        g0_ttheta_schwads_sph_rhorho =
     -   0
        g0_ttheta_schwads_sph_rhotheta =
     -   0
        g0_ttheta_schwads_sph_rhophi =
     -   0
        g0_ttheta_schwads_sph_thetatheta =
     -   0
        g0_ttheta_schwads_sph_thetaphi =
     -   0
        g0_ttheta_schwads_sph_phiphi =
     -   0

        g0_tphi_schwads_sph_t  =
     -   0
        g0_tphi_schwads_sph_rho  =
     -   0
        g0_tphi_schwads_sph_theta  =
     -   0
        g0_tphi_schwads_sph_phi  =
     -   0
        g0_tphi_schwads_sph_tt =
     -   0
        g0_tphi_schwads_sph_trho =
     -   0
        g0_tphi_schwads_sph_ttheta =
     -   0
        g0_tphi_schwads_sph_tphi =
     -   0
        g0_tphi_schwads_sph_rhorho =
     -   0
        g0_tphi_schwads_sph_rhotheta =
     -   0
        g0_tphi_schwads_sph_rhophi =
     -   0
        g0_tphi_schwads_sph_thetatheta =
     -   0
        g0_tphi_schwads_sph_thetaphi =
     -   0
        g0_tphi_schwads_sph_phiphi =
     -   0

        g0_rhorho_schwads_sph_t  =
     -   0
        g0_rhorho_schwads_sph_rho  =
     -   (16*M0*(-1 + rho0**2))/(1 + rho0**2)**3 - 
     -  (8*M0)/(1 + rho0**2)**2 + 
     -  (4*M0*(-1 + rho0**2))/
     -   (rho0**2*(1 + rho0**2)**2) - 
     -  (4*(1 + rho0**2)**2*
     -     (8*rho0 + 4*rho0*(-1 + rho0**2)))/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)**2) + 
     -  (16*rho0*(1 + rho0**2))/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)) - 
     -  (16*rho0*(1 + rho0**2)**2)/
     -   ((-1 + rho0**2)**3*
     -     (4*rho0**2 + (-1 + rho0**2)**2))
        g0_rhorho_schwads_sph_theta  =
     -   0
        g0_rhorho_schwads_sph_phi  =
     -   0
        g0_rhorho_schwads_sph_tt =
     -   0
        g0_rhorho_schwads_sph_trho =
     -   0
        g0_rhorho_schwads_sph_ttheta =
     -   0
        g0_rhorho_schwads_sph_tphi =
     -   0
        g0_rhorho_schwads_sph_rhorho =
     -   (-96*M0*rho0*(-1 + rho0**2))/
     -   (1 + rho0**2)**4 + 
     -  (64*M0*rho0)/(1 + rho0**2)**3 - 
     -  (16*M0*(-1 + rho0**2))/
     -   (rho0*(1 + rho0**2)**3) + 
     -  (8*M0)/(rho0*(1 + rho0**2)**2) - 
     -  (8*M0*(-1 + rho0**2))/
     -   (rho0**3*(1 + rho0**2)**2) + 
     -  (8*(1 + rho0**2)**2*
     -     (8*rho0 + 4*rho0*(-1 + rho0**2))**2)/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)**3) - 
     -  (4*(1 + rho0**2)**2*
     -     (8 + 8*rho0**2 + 4*(-1 + rho0**2)))/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)**2) - 
     -  (32*rho0*(1 + rho0**2)*
     -     (8*rho0 + 4*rho0*(-1 + rho0**2)))/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)**2) + 
     -  (32*rho0*(1 + rho0**2)**2*
     -     (8*rho0 + 4*rho0*(-1 + rho0**2)))/
     -   ((-1 + rho0**2)**3*
     -     (4*rho0**2 + (-1 + rho0**2)**2)**2) + 
     -  (32*rho0**2)/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)) - 
     -  (128*rho0**2*(1 + rho0**2))/
     -   ((-1 + rho0**2)**3*
     -     (4*rho0**2 + (-1 + rho0**2)**2)) + 
     -  (16*(1 + rho0**2))/
     -   ((-1 + rho0**2)**2*
     -     (4*rho0**2 + (-1 + rho0**2)**2)) + 
     -  (96*rho0**2*(1 + rho0**2)**2)/
     -   ((-1 + rho0**2)**4*
     -     (4*rho0**2 + (-1 + rho0**2)**2)) - 
     -  (16*(1 + rho0**2)**2)/
     -   ((-1 + rho0**2)**3*
     -     (4*rho0**2 + (-1 + rho0**2)**2))
        g0_rhorho_schwads_sph_rhotheta =
     -   0
        g0_rhorho_schwads_sph_rhophi =
     -   0
        g0_rhorho_schwads_sph_thetatheta =
     -   0
        g0_rhorho_schwads_sph_thetaphi =
     -   0
        g0_rhorho_schwads_sph_phiphi =
     -   0

        g0_rhotheta_schwads_sph_t  =
     -   0
        g0_rhotheta_schwads_sph_rho  =
     -   0
        g0_rhotheta_schwads_sph_theta  =
     -   0
        g0_rhotheta_schwads_sph_phi  =
     -   0
        g0_rhotheta_schwads_sph_tt =
     -   0
        g0_rhotheta_schwads_sph_trho =
     -   0
        g0_rhotheta_schwads_sph_ttheta =
     -   0
        g0_rhotheta_schwads_sph_tphi =
     -   0
        g0_rhotheta_schwads_sph_rhorho =
     -   0
        g0_rhotheta_schwads_sph_rhotheta =
     -   0
        g0_rhotheta_schwads_sph_rhophi =
     -   0
        g0_rhotheta_schwads_sph_thetatheta =
     -   0
        g0_rhotheta_schwads_sph_thetaphi =
     -   0
        g0_rhotheta_schwads_sph_phiphi =
     -   0

        g0_rhophi_schwads_sph_t  =
     -   0
        g0_rhophi_schwads_sph_rho  =
     -   0
        g0_rhophi_schwads_sph_theta  =
     -   0
        g0_rhophi_schwads_sph_phi  =
     -   0
        g0_rhophi_schwads_sph_tt =
     -   0
        g0_rhophi_schwads_sph_trho =
     -   0
        g0_rhophi_schwads_sph_ttheta =
     -   0
        g0_rhophi_schwads_sph_tphi =
     -   0
        g0_rhophi_schwads_sph_rhorho =
     -   0
        g0_rhophi_schwads_sph_rhotheta =
     -   0
        g0_rhophi_schwads_sph_rhophi =
     -   0
        g0_rhophi_schwads_sph_thetatheta =
     -   0
        g0_rhophi_schwads_sph_thetaphi =
     -   0
        g0_rhophi_schwads_sph_phiphi =
     -   0

        g0_thetatheta_schwads_sph_t  =
     -   0
        g0_thetatheta_schwads_sph_rho  =
     -   (-16*rho0**3)/(-1 + rho0**2)**3 + 
     -  (8*rho0)/(-1 + rho0**2)**2
        g0_thetatheta_schwads_sph_theta  =
     -   0
        g0_thetatheta_schwads_sph_phi  =
     -   0
        g0_thetatheta_schwads_sph_tt =
     -   0
        g0_thetatheta_schwads_sph_trho =
     -   0
        g0_thetatheta_schwads_sph_ttheta =
     -   0
        g0_thetatheta_schwads_sph_tphi =
     -   0
        g0_thetatheta_schwads_sph_rhorho =
     -   (96*rho0**4)/(-1 + rho0**2)**4 - 
     -  (80*rho0**2)/(-1 + rho0**2)**3 + 
     -  8/(-1 + rho0**2)**2
        g0_thetatheta_schwads_sph_rhotheta =
     -   0
        g0_thetatheta_schwads_sph_rhophi =
     -   0
        g0_thetatheta_schwads_sph_thetatheta =
     -   0
        g0_thetatheta_schwads_sph_thetaphi =
     -   0
        g0_thetatheta_schwads_sph_phiphi =
     -   0

        g0_thetaphi_schwads_sph_t  =
     -   0
        g0_thetaphi_schwads_sph_rho  =
     -   0
        g0_thetaphi_schwads_sph_theta  =
     -   0
        g0_thetaphi_schwads_sph_phi  =
     -   0
        g0_thetaphi_schwads_sph_tt =
     -   0
        g0_thetaphi_schwads_sph_trho =
     -   0
        g0_thetaphi_schwads_sph_ttheta =
     -   0
        g0_thetaphi_schwads_sph_tphi =
     -   0
        g0_thetaphi_schwads_sph_rhorho =
     -   0
        g0_thetaphi_schwads_sph_rhotheta =
     -   0
        g0_thetaphi_schwads_sph_rhophi =
     -   0
        g0_thetaphi_schwads_sph_thetatheta =
     -   0
        g0_thetaphi_schwads_sph_thetaphi =
     -   0
        g0_thetaphi_schwads_sph_phiphi =
     -   0

        g0_phiphi_schwads_sph_t  =
     -   0
        g0_phiphi_schwads_sph_rho  =
     -   (-16*rho0**3*Sin(theta0)**2)/
     -   (-1 + rho0**2)**3 + 
     -  (8*rho0*Sin(theta0)**2)/(-1 + rho0**2)**2
        g0_phiphi_schwads_sph_theta  =
     -   (8*rho0**2*Cos(theta0)*Sin(theta0))/
     -      (-1 + rho0**2)**2
        g0_phiphi_schwads_sph_phi  =
     -   0
        g0_phiphi_schwads_sph_tt =
     -   0
        g0_phiphi_schwads_sph_trho =
     -   0
        g0_phiphi_schwads_sph_ttheta =
     -   0
        g0_phiphi_schwads_sph_tphi =
     -   0
        g0_phiphi_schwads_sph_rhorho =
     -   (96*rho0**4*Sin(theta0)**2)/
     -   (-1 + rho0**2)**4 - 
     -  (80*rho0**2*Sin(theta0)**2)/
     -   (-1 + rho0**2)**3 + 
     -  (8*Sin(theta0)**2)/(-1 + rho0**2)**2
        g0_phiphi_schwads_sph_rhotheta =
     -   (-32*rho0**3*Cos(theta0)*Sin(theta0))/
     -   (-1 + rho0**2)**3 + 
     -  (16*rho0*Cos(theta0)*Sin(theta0))/
     -   (-1 + rho0**2)**2
        g0_phiphi_schwads_sph_rhophi =
     -   0
        g0_phiphi_schwads_sph_thetatheta =
     -   (8*rho0**2*Cos(theta0)**2)/
     -      (-1 + rho0**2)**2 - 
     -   (8*rho0**2*Sin(theta0)**2)/
     -      (-1 + rho0**2)**2
        g0_phiphi_schwads_sph_thetaphi =
     -   0
        g0_phiphi_schwads_sph_phiphi =
     -   0     


        ! give values to the schwads metric
        gschwads_ll_sph(1,1)=g0_tt_schwads_sph0
        gschwads_ll_sph(1,2)=g0_trho_schwads_sph0
        gschwads_ll_sph(1,3)=g0_ttheta_schwads_sph0
        gschwads_ll_sph(1,4)=g0_tphi_schwads_sph0
        gschwads_ll_sph(2,2)=g0_rhorho_schwads_sph0
        gschwads_ll_sph(2,3)=g0_rhotheta_schwads_sph0
        gschwads_ll_sph(2,4)=g0_rhophi_schwads_sph0
        gschwads_ll_sph(3,3)=g0_thetatheta_schwads_sph0
        gschwads_ll_sph(3,4)=g0_thetaphi_schwads_sph0
        gschwads_ll_sph(4,4)=g0_phiphi_schwads_sph0


        gschwads_uu_sph(1,1)=g0u_tt_schwads_sph0
        gschwads_uu_sph(1,2)=g0u_trho_schwads_sph0
        gschwads_uu_sph(1,3)=g0u_ttheta_schwads_sph0
        gschwads_uu_sph(1,4)=g0u_tphi_schwads_sph0
        gschwads_uu_sph(2,2)=g0u_rhorho_schwads_sph0
        gschwads_uu_sph(2,3)=g0u_rhotheta_schwads_sph0
        gschwads_uu_sph(2,4)=g0u_rhophi_schwads_sph0
        gschwads_uu_sph(3,3)=g0u_thetatheta_schwads_sph0
        gschwads_uu_sph(3,4)=g0u_thetaphi_schwads_sph0
        gschwads_uu_sph(4,4)=g0u_phiphi_schwads_sph0



        gschwads_ll_sph_x(1,1,1)   =g0_tt_schwads_sph_t
        gschwads_ll_sph_x(1,1,2)   =g0_tt_schwads_sph_rho
        gschwads_ll_sph_x(1,1,3)   =g0_tt_schwads_sph_theta
        gschwads_ll_sph_x(1,1,4)   =g0_tt_schwads_sph_phi
        gschwads_ll_sph_xx(1,1,1,1)=g0_tt_schwads_sph_tt
        gschwads_ll_sph_xx(1,1,1,2)=g0_tt_schwads_sph_trho
        gschwads_ll_sph_xx(1,1,1,3)=g0_tt_schwads_sph_ttheta
        gschwads_ll_sph_xx(1,1,1,4)=g0_tt_schwads_sph_tphi
        gschwads_ll_sph_xx(1,1,2,2)=g0_tt_schwads_sph_rhorho
        gschwads_ll_sph_xx(1,1,2,3)=g0_tt_schwads_sph_rhotheta
        gschwads_ll_sph_xx(1,1,2,4)=g0_tt_schwads_sph_rhophi
        gschwads_ll_sph_xx(1,1,3,3)=g0_tt_schwads_sph_thetatheta
        gschwads_ll_sph_xx(1,1,3,4)=g0_tt_schwads_sph_thetaphi
        gschwads_ll_sph_xx(1,1,4,4)=g0_tt_schwads_sph_phiphi

        gschwads_ll_sph_x(1,2,1)   =g0_trho_schwads_sph_t
        gschwads_ll_sph_x(1,2,2)   =g0_trho_schwads_sph_rho
        gschwads_ll_sph_x(1,2,3)   =g0_trho_schwads_sph_theta
        gschwads_ll_sph_x(1,2,4)   =g0_trho_schwads_sph_phi
        gschwads_ll_sph_xx(1,2,1,1)=g0_trho_schwads_sph_tt
        gschwads_ll_sph_xx(1,2,1,2)=g0_trho_schwads_sph_trho
        gschwads_ll_sph_xx(1,2,1,3)=g0_trho_schwads_sph_ttheta
        gschwads_ll_sph_xx(1,2,1,4)=g0_trho_schwads_sph_tphi
        gschwads_ll_sph_xx(1,2,2,2)=g0_trho_schwads_sph_rhorho
        gschwads_ll_sph_xx(1,2,2,3)=g0_trho_schwads_sph_rhotheta
        gschwads_ll_sph_xx(1,2,2,4)=g0_trho_schwads_sph_rhophi
        gschwads_ll_sph_xx(1,2,3,3)=g0_trho_schwads_sph_thetatheta
        gschwads_ll_sph_xx(1,2,3,4)=g0_trho_schwads_sph_thetaphi
        gschwads_ll_sph_xx(1,2,4,4)=g0_trho_schwads_sph_phiphi

        gschwads_ll_sph_x(1,3,1)   =g0_ttheta_schwads_sph_t
        gschwads_ll_sph_x(1,3,2)   =g0_ttheta_schwads_sph_rho
        gschwads_ll_sph_x(1,3,3)   =g0_ttheta_schwads_sph_theta
        gschwads_ll_sph_x(1,3,4)   =g0_ttheta_schwads_sph_phi
        gschwads_ll_sph_xx(1,3,1,1)=g0_ttheta_schwads_sph_tt
        gschwads_ll_sph_xx(1,3,1,2)=g0_ttheta_schwads_sph_trho
        gschwads_ll_sph_xx(1,3,1,3)=g0_ttheta_schwads_sph_ttheta
        gschwads_ll_sph_xx(1,3,1,4)=g0_ttheta_schwads_sph_tphi
        gschwads_ll_sph_xx(1,3,2,2)=g0_ttheta_schwads_sph_rhorho
        gschwads_ll_sph_xx(1,3,2,3)=g0_ttheta_schwads_sph_rhotheta
        gschwads_ll_sph_xx(1,3,2,4)=g0_ttheta_schwads_sph_rhophi
        gschwads_ll_sph_xx(1,3,3,3)=g0_ttheta_schwads_sph_thetatheta
        gschwads_ll_sph_xx(1,3,3,4)=g0_ttheta_schwads_sph_thetaphi
        gschwads_ll_sph_xx(1,3,4,4)=g0_ttheta_schwads_sph_phiphi

        gschwads_ll_sph_x(1,4,1)   =g0_tphi_schwads_sph_t
        gschwads_ll_sph_x(1,4,2)   =g0_tphi_schwads_sph_rho
        gschwads_ll_sph_x(1,4,3)   =g0_tphi_schwads_sph_theta
        gschwads_ll_sph_x(1,4,4)   =g0_tphi_schwads_sph_phi
        gschwads_ll_sph_xx(1,4,1,1)=g0_tphi_schwads_sph_tt
        gschwads_ll_sph_xx(1,4,1,2)=g0_tphi_schwads_sph_trho
        gschwads_ll_sph_xx(1,4,1,3)=g0_tphi_schwads_sph_ttheta
        gschwads_ll_sph_xx(1,4,1,4)=g0_tphi_schwads_sph_tphi
        gschwads_ll_sph_xx(1,4,2,2)=g0_tphi_schwads_sph_rhorho
        gschwads_ll_sph_xx(1,4,2,3)=g0_tphi_schwads_sph_rhotheta
        gschwads_ll_sph_xx(1,4,2,4)=g0_tphi_schwads_sph_rhophi
        gschwads_ll_sph_xx(1,4,3,3)=g0_tphi_schwads_sph_thetatheta
        gschwads_ll_sph_xx(1,4,3,4)=g0_tphi_schwads_sph_thetaphi
        gschwads_ll_sph_xx(1,4,4,4)=g0_tphi_schwads_sph_phiphi

        gschwads_ll_sph_x(2,2,1)   =g0_rhorho_schwads_sph_t
        gschwads_ll_sph_x(2,2,2)   =g0_rhorho_schwads_sph_rho
        gschwads_ll_sph_x(2,2,3)   =g0_rhorho_schwads_sph_theta
        gschwads_ll_sph_x(2,2,4)   =g0_rhorho_schwads_sph_phi
        gschwads_ll_sph_xx(2,2,1,1)=g0_rhorho_schwads_sph_tt
        gschwads_ll_sph_xx(2,2,1,2)=g0_rhorho_schwads_sph_trho
        gschwads_ll_sph_xx(2,2,1,3)=g0_rhorho_schwads_sph_ttheta
        gschwads_ll_sph_xx(2,2,1,4)=g0_rhorho_schwads_sph_tphi
        gschwads_ll_sph_xx(2,2,2,2)=g0_rhorho_schwads_sph_rhorho
        gschwads_ll_sph_xx(2,2,2,3)=g0_rhorho_schwads_sph_rhotheta
        gschwads_ll_sph_xx(2,2,2,4)=g0_rhorho_schwads_sph_rhophi
        gschwads_ll_sph_xx(2,2,3,3)=g0_rhorho_schwads_sph_thetatheta
        gschwads_ll_sph_xx(2,2,3,4)=g0_rhorho_schwads_sph_thetaphi
        gschwads_ll_sph_xx(2,2,4,4)=g0_rhorho_schwads_sph_phiphi

        gschwads_ll_sph_x(2,3,1)   =g0_rhotheta_schwads_sph_t
        gschwads_ll_sph_x(2,3,2)   =g0_rhotheta_schwads_sph_rho
        gschwads_ll_sph_x(2,3,3)   =g0_rhotheta_schwads_sph_theta
        gschwads_ll_sph_x(2,3,4)   =g0_rhotheta_schwads_sph_phi
        gschwads_ll_sph_xx(2,3,1,1)=g0_rhotheta_schwads_sph_tt
        gschwads_ll_sph_xx(2,3,1,2)=g0_rhotheta_schwads_sph_trho
        gschwads_ll_sph_xx(2,3,1,3)=g0_rhotheta_schwads_sph_ttheta
        gschwads_ll_sph_xx(2,3,1,4)=g0_rhotheta_schwads_sph_tphi
        gschwads_ll_sph_xx(2,3,2,2)=g0_rhotheta_schwads_sph_rhorho
        gschwads_ll_sph_xx(2,3,2,3)=g0_rhotheta_schwads_sph_rhotheta
        gschwads_ll_sph_xx(2,3,2,4)=g0_rhotheta_schwads_sph_rhophi
        gschwads_ll_sph_xx(2,3,3,3)=g0_rhotheta_schwads_sph_thetatheta
        gschwads_ll_sph_xx(2,3,3,4)=g0_rhotheta_schwads_sph_thetaphi
        gschwads_ll_sph_xx(2,3,4,4)=g0_rhotheta_schwads_sph_phiphi

        gschwads_ll_sph_x(2,4,1)   =g0_rhophi_schwads_sph_t
        gschwads_ll_sph_x(2,4,2)   =g0_rhophi_schwads_sph_rho
        gschwads_ll_sph_x(2,4,3)   =g0_rhophi_schwads_sph_theta
        gschwads_ll_sph_x(2,4,4)   =g0_rhophi_schwads_sph_phi
        gschwads_ll_sph_xx(2,4,1,1)=g0_rhophi_schwads_sph_tt
        gschwads_ll_sph_xx(2,4,1,2)=g0_rhophi_schwads_sph_trho
        gschwads_ll_sph_xx(2,4,1,3)=g0_rhophi_schwads_sph_ttheta
        gschwads_ll_sph_xx(2,4,1,4)=g0_rhophi_schwads_sph_tphi
        gschwads_ll_sph_xx(2,4,2,2)=g0_rhophi_schwads_sph_rhorho
        gschwads_ll_sph_xx(2,4,2,3)=g0_rhophi_schwads_sph_rhotheta
        gschwads_ll_sph_xx(2,4,2,4)=g0_rhophi_schwads_sph_rhophi
        gschwads_ll_sph_xx(2,4,3,3)=g0_rhophi_schwads_sph_thetatheta
        gschwads_ll_sph_xx(2,4,3,4)=g0_rhophi_schwads_sph_thetaphi
        gschwads_ll_sph_xx(2,4,4,4)=g0_rhophi_schwads_sph_phiphi

        gschwads_ll_sph_x(3,3,1)   =g0_thetatheta_schwads_sph_t
        gschwads_ll_sph_x(3,3,2)   =g0_thetatheta_schwads_sph_rho
        gschwads_ll_sph_x(3,3,3)   =g0_thetatheta_schwads_sph_theta
        gschwads_ll_sph_x(3,3,4)   =g0_thetatheta_schwads_sph_phi
        gschwads_ll_sph_xx(3,3,1,1)=g0_thetatheta_schwads_sph_tt
        gschwads_ll_sph_xx(3,3,1,2)=g0_thetatheta_schwads_sph_trho
        gschwads_ll_sph_xx(3,3,1,3)=g0_thetatheta_schwads_sph_ttheta
        gschwads_ll_sph_xx(3,3,1,4)=g0_thetatheta_schwads_sph_tphi
        gschwads_ll_sph_xx(3,3,2,2)=g0_thetatheta_schwads_sph_rhorho
        gschwads_ll_sph_xx(3,3,2,3)=g0_thetatheta_schwads_sph_rhotheta
        gschwads_ll_sph_xx(3,3,2,4)=g0_thetatheta_schwads_sph_rhophi
        gschwads_ll_sph_xx(3,3,3,3)=g0_thetatheta_schwads_sph_thetatheta
        gschwads_ll_sph_xx(3,3,3,4)=g0_thetatheta_schwads_sph_thetaphi
        gschwads_ll_sph_xx(3,3,4,4)=g0_thetatheta_schwads_sph_phiphi

        gschwads_ll_sph_x(3,4,1)   =g0_thetaphi_schwads_sph_t
        gschwads_ll_sph_x(3,4,2)   =g0_thetaphi_schwads_sph_rho
        gschwads_ll_sph_x(3,4,3)   =g0_thetaphi_schwads_sph_theta
        gschwads_ll_sph_x(3,4,4)   =g0_thetaphi_schwads_sph_phi
        gschwads_ll_sph_xx(3,4,1,1)=g0_thetaphi_schwads_sph_tt
        gschwads_ll_sph_xx(3,4,1,2)=g0_thetaphi_schwads_sph_trho
        gschwads_ll_sph_xx(3,4,1,3)=g0_thetaphi_schwads_sph_ttheta
        gschwads_ll_sph_xx(3,4,1,4)=g0_thetaphi_schwads_sph_tphi
        gschwads_ll_sph_xx(3,4,2,2)=g0_thetaphi_schwads_sph_rhorho
        gschwads_ll_sph_xx(3,4,2,3)=g0_thetaphi_schwads_sph_rhotheta
        gschwads_ll_sph_xx(3,4,2,4)=g0_thetaphi_schwads_sph_rhophi
        gschwads_ll_sph_xx(3,4,3,3)=g0_thetaphi_schwads_sph_thetatheta
        gschwads_ll_sph_xx(3,4,3,4)=g0_thetaphi_schwads_sph_thetaphi
        gschwads_ll_sph_xx(3,4,4,4)=g0_thetaphi_schwads_sph_phiphi

        gschwads_ll_sph_x(4,4,1)   =g0_phiphi_schwads_sph_t
        gschwads_ll_sph_x(4,4,2)   =g0_phiphi_schwads_sph_rho
        gschwads_ll_sph_x(4,4,3)   =g0_phiphi_schwads_sph_theta
        gschwads_ll_sph_x(4,4,4)   =g0_phiphi_schwads_sph_phi
        gschwads_ll_sph_xx(4,4,1,1)=g0_phiphi_schwads_sph_tt
        gschwads_ll_sph_xx(4,4,1,2)=g0_phiphi_schwads_sph_trho
        gschwads_ll_sph_xx(4,4,1,3)=g0_phiphi_schwads_sph_ttheta
        gschwads_ll_sph_xx(4,4,1,4)=g0_phiphi_schwads_sph_tphi
        gschwads_ll_sph_xx(4,4,2,2)=g0_phiphi_schwads_sph_rhorho
        gschwads_ll_sph_xx(4,4,2,3)=g0_phiphi_schwads_sph_rhotheta
        gschwads_ll_sph_xx(4,4,2,4)=g0_phiphi_schwads_sph_rhophi
        gschwads_ll_sph_xx(4,4,3,3)=g0_phiphi_schwads_sph_thetatheta
        gschwads_ll_sph_xx(4,4,3,4)=g0_phiphi_schwads_sph_thetaphi
        gschwads_ll_sph_xx(4,4,4,4)=g0_phiphi_schwads_sph_phiphi


        do a=1,3
          do b=a+1,4
            gschwads_ll_sph(b,a)=gschwads_ll_sph(a,b)
            gschwads_uu_sph(b,a)=gschwads_uu_sph(a,b)
            do c=1,4
              gschwads_ll_sph_x(b,a,c)=gschwads_ll_sph_x(a,b,c)
              do d=1,4
                gschwads_ll_sph_xx(b,a,c,d)=gschwads_ll_sph_xx(a,b,c,d)
                gschwads_ll_sph_xx(c,d,b,a)=gschwads_ll_sph_xx(c,d,a,b)
              end do
            end do
          end do
        end do

        !define transformation matrix between spherical to Cartesian coordinates, 
        !e.g. dxsph_dxcar(2,3)=drho/dtheta, d2xsph_dxcardxcar(2,3,4)=d^2rho/(dtheta dphi), 
        ! d3xsph_dxcardxcardxcar(2,3,4,2)=d^2rho/(dtheta dphi drho)

        dxsph_dxcar(1,1)=1
        dxsph_dxcar(1,2)=0
        dxsph_dxcar(1,3)=0
        dxsph_dxcar(1,4)=0

        dxsph_dxcar(2,1)=0
        dxsph_dxcar(2,2)=x0/rho0
        dxsph_dxcar(2,3)=y0/rho0
        dxsph_dxcar(2,4)=z0/rho0

        dxsph_dxcar(3,1)=0
        dxsph_dxcar(3,2)=-(Sqrt(rho0**2 - x0**2)/rho0**2)
        dxsph_dxcar(3,3)=(x0*y0)/(rho0**2*Sqrt(rho0**2 - x0**2))
        dxsph_dxcar(3,4)=(x0*z0)/(rho0**2*Sqrt(rho0**2 - x0**2))

        dxsph_dxcar(4,1)=0
        dxsph_dxcar(4,2)=0
        dxsph_dxcar(4,3)=-(z0/(y0**2 + z0**2))
        dxsph_dxcar(4,4)=y0/(y0**2 + z0**2)

        d2xsph_dxcardxcar(1,1,1)=
     -   0
        d2xsph_dxcardxcar(1,1,2)=
     -   0
        d2xsph_dxcardxcar(1,1,3)=
     -   0
        d2xsph_dxcardxcar(1,1,4)=
     -   0
        d2xsph_dxcardxcar(1,2,2)=
     -   0
        d2xsph_dxcardxcar(1,2,3)=
     -   0
        d2xsph_dxcardxcar(1,2,4)=
     -   0
        d2xsph_dxcardxcar(1,3,3)=
     -   0
        d2xsph_dxcardxcar(1,3,4)=
     -   0
        d2xsph_dxcardxcar(1,4,4)=
     -   0

        d2xsph_dxcardxcar(2,1,1)=
     -   0
        d2xsph_dxcardxcar(2,1,2)=
     -   0
        d2xsph_dxcardxcar(2,1,3)=
     -   0
        d2xsph_dxcardxcar(2,1,4)=
     -   0
        d2xsph_dxcardxcar(2,2,2)=
     -   (rho0**2 - x0**2)/rho0**3
        d2xsph_dxcardxcar(2,2,3)=
     -   -((x0*y0)/rho0**3)
        d2xsph_dxcardxcar(2,2,4)=
     -   -((x0*z0)/rho0**3)
        d2xsph_dxcardxcar(2,3,3)=
     -   (rho0**2 - y0**2)/rho0**3
        d2xsph_dxcardxcar(2,3,4)=
     -   -((y0*z0)/rho0**3)
        d2xsph_dxcardxcar(2,4,4)=
     -   (rho0**2 - z0**2)/rho0**3

        d2xsph_dxcardxcar(3,1,1)=
     -   0
        d2xsph_dxcardxcar(3,1,2)=
     -   0
        d2xsph_dxcardxcar(3,1,3)=
     -   0
        d2xsph_dxcardxcar(3,1,4)=
     -   0
        d2xsph_dxcardxcar(3,2,2)=
     -   (2*x0*Sqrt(rho0**2 - x0**2))/rho0**4
        d2xsph_dxcardxcar(3,2,3)=
     -   ((rho0**2 - 2*x0**2)*y0)/
     -  (rho0**4*Sqrt(rho0**2 - x0**2))
        d2xsph_dxcardxcar(3,2,4)=
     -   ((rho0**2 - 2*x0**2)*z0)/
     -  (rho0**4*Sqrt(rho0**2 - x0**2))
        d2xsph_dxcardxcar(3,3,3)=
     -   (x0*(rho0**4 + 2*x0**2*y0**2 - 
     -      rho0**2*(x0**2 + 3*y0**2)))/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)
        d2xsph_dxcardxcar(3,3,4)=
     -   ((-3*rho0**2*x0 + 2*x0**3)*y0*z0)/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)
        d2xsph_dxcardxcar(3,4,4)=
     -   (x0*(rho0**4 + 2*x0**2*z0**2 - 
     -      rho0**2*(x0**2 + 3*z0**2)))/
     -  (rho0**4*(rho0**2 - x0**2)**1.5)

        d2xsph_dxcardxcar(4,1,1)=
     -   0
        d2xsph_dxcardxcar(4,1,2)=
     -   0
        d2xsph_dxcardxcar(4,1,3)=
     -   0
        d2xsph_dxcardxcar(4,1,4)=
     -   0
        d2xsph_dxcardxcar(4,2,2)=
     -   0
        d2xsph_dxcardxcar(4,2,3)=
     -   0
        d2xsph_dxcardxcar(4,2,4)=
     -   0
        d2xsph_dxcardxcar(4,3,3)=
     -   (2*y0*z0)/(y0**2 + z0**2)**2
        d2xsph_dxcardxcar(4,3,4)=
     -   (-y0**2 + z0**2)/(y0**2 + z0**2)**2
        d2xsph_dxcardxcar(4,4,4)=
     -   (-2*y0*z0)/(y0**2 + z0**2)**2

        d3xsph_dxcardxcardxcar(1,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,2,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,3,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(1,4,4,4)=
     -   0



        d3xsph_dxcardxcardxcar(2,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(2,2,2,2)=
     -   (3*x0*(-rho0 + x0)*(rho0 + x0))/rho0**5
        d3xsph_dxcardxcardxcar(2,2,2,3)=
     -   -(((rho0**2 - 3*x0**2)*y0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,2,4)=
     -   -(((rho0**2 - 3*x0**2)*z0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,3,3)=
     -   -((x0*(rho0**2 - 3*y0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,2,3,4)=
     -   (3*x0*y0*z0)/rho0**5
        d3xsph_dxcardxcardxcar(2,2,4,4)=
     -   -((x0*(rho0**2 - 3*z0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,3,3,3)=
     -   (3*y0*(-rho0 + y0)*(rho0 + y0))/rho0**5
        d3xsph_dxcardxcardxcar(2,3,3,4)=
     -   -(((rho0**2 - 3*y0**2)*z0)/rho0**5)
        d3xsph_dxcardxcardxcar(2,3,4,4)=
     -   -((y0*(rho0**2 - 3*z0**2))/rho0**5)
        d3xsph_dxcardxcardxcar(2,4,4,4)=
     -   (3*z0*(-rho0 + z0)*(rho0 + z0))/rho0**5

        d3xsph_dxcardxcardxcar(3,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(3,2,2,2)=
     -   (2*(rho0**2 - 4*x0**2)*
     -   Sqrt(rho0**2 - x0**2))/rho0**6
        d3xsph_dxcardxcardxcar(3,2,2,3)=
     -   (-6*rho0**2*x0*y0 + 8*x0**3*y0)/
     -  (rho0**6*Sqrt(rho0**2 - x0**2))
        d3xsph_dxcardxcardxcar(3,2,2,4)=
     -   (-6*rho0**2*x0*z0 + 8*x0**3*z0)/
     -  (rho0**6*Sqrt(rho0**2 - x0**2))
        d3xsph_dxcardxcardxcar(3,2,3,3)=
     -   (rho0**6 - 8*x0**4*y0**2 - 
     -    3*rho0**4*(x0**2 + y0**2) + 
     -    2*rho0**2*x0**2*(x0**2 + 6*y0**2))/
     -  (rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,2,3,4)=
     -   ((-3*rho0**4 + 12*rho0**2*x0**2 - 8*x0**4)*y0*
     -    z0)/(rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,2,4,4)=
     -   (rho0**6 - 8*x0**4*z0**2 - 
     -    3*rho0**4*(x0**2 + z0**2) + 
     -    2*rho0**2*x0**2*(x0**2 + 6*z0**2))/
     -  (rho0**6*(rho0**2 - x0**2)**1.5)
        d3xsph_dxcardxcardxcar(3,3,3,3)=
     -   (x0*y0*(-9*rho0**6 + 8*x0**4*y0**2 + 
     -      15*rho0**4*(x0**2 + y0**2) - 
     -      2*rho0**2*x0**2*(3*x0**2 + 10*y0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,3,3,4)=
     -   (x0*(-3*rho0**6 + 8*x0**4*y0**2 + 
     -      5*rho0**4*(x0**2 + 3*y0**2) - 
     -      2*rho0**2*x0**2*(x0**2 + 10*y0**2))*z0)/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,3,4,4)=
     -   (x0*y0*(-3*rho0**6 + 8*x0**4*z0**2 + 
     -      5*rho0**4*(x0**2 + 3*z0**2) - 
     -      2*rho0**2*x0**2*(x0**2 + 10*z0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)
        d3xsph_dxcardxcardxcar(3,4,4,4)=
     -   (x0*z0*(-9*rho0**6 + 8*x0**4*z0**2 + 
     -      15*rho0**4*(x0**2 + z0**2) - 
     -      2*rho0**2*x0**2*(3*x0**2 + 10*z0**2)))/
     -  (rho0**6*(rho0**2 - x0**2)**2.5)

        d3xsph_dxcardxcardxcar(4,1,1,1)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,1,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,1,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,2)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,2,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,3,3)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,3,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,2,4,4)=
     -   0
        d3xsph_dxcardxcardxcar(4,3,3,3)=
     -   (2*z0*(-3*y0**2 + z0**2))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,3,3,4)=
     -   (2*y0*(y0**2 - 3*z0**2))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,3,4,4)=
     -   (-2*(-3*y0**2*z0 + z0**3))/(y0**2 + z0**2)**3
        d3xsph_dxcardxcardxcar(4,4,4,4)=
     -   (-2*y0*(y0**2 - 3*z0**2))/(y0**2 + z0**2)**3


        do a=1,4
         do b=1,3
          do c=b+1,4
           d2xsph_dxcardxcar(a,c,b)=d2xsph_dxcardxcar(a,b,c)
          end do
         end do
        end do


        do a=1,4
         do b=1,4
          do c=1,4
           do d=1,4
            if ((max(b,c,d).eq.d).and.(max(b,c).eq.c)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,b,c,d)
            else if ((max(b,c,d).eq.d).and.(max(b,c).eq.b)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,c,b,d)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,c,b,d)

            else if ((max(b,c,d).eq.c).and.(max(b,d).eq.d)) then
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,b,d,c)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,b,d,c)
            else if ((max(b,c,d).eq.c).and.(max(b,d).eq.b)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,d,b,c)

            else if ((max(b,c,d).eq.b).and.(max(c,d).eq.d)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
             d3xsph_dxcardxcardxcar(a,d,c,b)=
     &        d3xsph_dxcardxcardxcar(a,c,d,b)
            else if ((max(b,c,d).eq.b).and.(max(c,d).eq.c)) then
             d3xsph_dxcardxcardxcar(a,b,d,c)=
     &        d3xsph_dxcardxcardxcar(a,d,c,b)
             d3xsph_dxcardxcardxcar(a,c,b,d)=
     &        d3xsph_dxcardxcardxcar(a,d,c,b)
            end if
           end do
          end do
         end do
        end do

        !compute Cartesian quantities in terms of spherical ones
        do a=1,4
         do b=1,4
          gschwads_ll(a,b)=0
          do c=1,4
            gschwads_ll_x(a,b,c)=0
           do d=1,4
            gschwads_ll_xx(a,b,c,d)=0
            gschwads_ll(a,b)=gschwads_ll(a,b)+
     &         dxsph_dxcar(c,a)*dxsph_dxcar(d,b)*gschwads_ll_sph(c,d)
            do e=1,4
             gschwads_ll_x(a,b,c)=gschwads_ll_x(a,b,c)+
     &     d2xsph_dxcardxcar(d,c,a)*dxsph_dxcar(e,b)
     &     *gschwads_ll_sph(d,e)+
     &     dxsph_dxcar(d,a)*d2xsph_dxcardxcar(e,c,b)
     &     *gschwads_ll_sph(d,e)
             do f=1,4
              gschwads_ll_x(a,b,c)=gschwads_ll_x(a,b,c)+
     &     dxsph_dxcar(d,a)*dxsph_dxcar(e,b)*dxsph_dxcar(f,c)
     &     *gschwads_ll_sph_x(d,e,f)
              gschwads_ll_xx(a,b,c,d)=gschwads_ll_xx(a,b,c,d)+
     &     d3xsph_dxcardxcardxcar(e,d,c,a)*dxsph_dxcar(f,b)
     &     *gschwads_ll_sph(e,f)+
     &     d2xsph_dxcardxcar(e,c,a)*d2xsph_dxcardxcar(f,d,b)
     &     *gschwads_ll_sph(e,f)+
     &     d2xsph_dxcardxcar(e,d,a)*d2xsph_dxcardxcar(f,c,b)
     &     *gschwads_ll_sph(e,f)+
     &     dxsph_dxcar(e,a)*d3xsph_dxcardxcardxcar(f,d,c,b)
     &     *gschwads_ll_sph(e,f)
              do g=1,4
           gschwads_ll_xx(a,b,c,d)=gschwads_ll_xx(a,b,c,d)+
     &     d2xsph_dxcardxcar(e,c,a)*dxsph_dxcar(f,b)*dxsph_dxcar(g,d)
     &     *gschwads_ll_sph_x(e,f,g)+
     &     dxsph_dxcar(e,a)*d2xsph_dxcardxcar(f,c,b)*dxsph_dxcar(g,d)
     &     *gschwads_ll_sph_x(e,f,g)+
     &     d2xsph_dxcardxcar(e,d,a)*dxsph_dxcar(f,b)*dxsph_dxcar(g,c)
     &     *gschwads_ll_sph_x(e,f,g)+
     &     dxsph_dxcar(e,a)*d2xsph_dxcardxcar(f,d,b)*dxsph_dxcar(g,c)
     &     *gschwads_ll_sph_x(e,f,g)+
     &     dxsph_dxcar(e,a)*dxsph_dxcar(f,b)*d2xsph_dxcardxcar(g,d,c)
     &     *gschwads_ll_sph_x(e,f,g)
               do h=1,4
                gschwads_ll_xx(a,b,c,d)=gschwads_ll_xx(a,b,c,d)+
     &     dxsph_dxcar(e,a)*dxsph_dxcar(f,b)*
     &     dxsph_dxcar(g,c)*dxsph_dxcar(h,d)*gschwads_ll_sph_xx(e,f,g,h)
               end do
              end do
             end do
            end do
           end do
          end do
         end do
        end do

       !some of the dxsph_dxcar diverge at y=z=0 so we need to consider this case separately
       ! we only need to consider the quantities with a y or z index
       if ((abs(y0).lt.10.0d0**(-10)).and.
     &     (abs(z0).lt.10.0d0**(-10))) then
        gschwads_ll(1,3)=
     -   0
        gschwads_ll(1,4)=
     -   0
        gschwads_ll(2,3)=
     -   0
        gschwads_ll(2,4)=
     -   0
        gschwads_ll(3,3)=
     -   4/(-1 + x0**2)**2
        gschwads_ll(3,4)=
     -   0
        gschwads_ll(4,4)=
     -   4/(-1 + x0**2)**2


        gschwads_ll_x(1,1,3)=
     -   0
        gschwads_ll_x(1,1,4)=
     -   0
        gschwads_ll_xx(1,1,1,3)=
     -   0
        gschwads_ll_xx(1,1,1,4)=
     -   0
        gschwads_ll_xx(1,1,2,3)=
     -   0
        gschwads_ll_xx(1,1,2,4)=
     -   0
        gschwads_ll_xx(1,1,3,3)=
     -   (-4*(1 + x0**2))/(-1 + x0**2)**2 + 
     -  (4*(1 + 2*x0**2 + x0**4))/
     -   (-1 + x0**2)**3 - M0/Abs(x0)**3 - 
     -  M0/Abs(x0)
        gschwads_ll_xx(1,1,3,4)=
     -   0
        gschwads_ll_xx(1,1,4,4)=
     -   (-4*(1 + x0**2))/(-1 + x0**2)**2 + 
     -  (4*(1 + 2*x0**2 + x0**4))/
     -   (-1 + x0**2)**3 - M0/Abs(x0)**3 - 
     -  M0/Abs(x0)

        gschwads_ll_x(1,2,3)=
     -   0
        gschwads_ll_x(1,2,4)=
     -   0
        gschwads_ll_xx(1,2,1,3)=
     -   0
        gschwads_ll_xx(1,2,1,4)=
     -   0
        gschwads_ll_xx(1,2,2,3)=
     -   0
        gschwads_ll_xx(1,2,2,4)=
     -   0
        gschwads_ll_xx(1,2,3,3)=
     -    (-4*M0*(1 + x0**2 - 3*x0**4 + x0**6))/
     -  (x0**3*Sqrt((-1 + x0**2)**2)*
     -    (1 + x0**2)**2)
        gschwads_ll_xx(1,2,3,4)=
     -   0
        gschwads_ll_xx(1,2,4,4)=
     -   (-4*M0*(1 + x0**2 - 3*x0**4 + x0**6))/
     -  (x0**3*Sqrt((-1 + x0**2)**2)*
     -    (1 + x0**2)**2)

        gschwads_ll_x(1,3,1)=
     -   0
        gschwads_ll_x(1,3,2)=
     -   0
        gschwads_ll_x(1,3,3)=
     -   (2*M0*Sqrt((-1 + x0**2)**2))/
     -      (x0**2*(1 + x0**2))
        gschwads_ll_x(1,3,4)=
     -   0
        gschwads_ll_xx(1,3,1,1)=
     -   0
        gschwads_ll_xx(1,3,1,2)=
     -   0
        gschwads_ll_xx(1,3,1,3)=
     -   0
        gschwads_ll_xx(1,3,1,4)=
     -   0
        gschwads_ll_xx(1,3,2,2)=
     -   0
        gschwads_ll_xx(1,3,2,3)=
     -    (-4*M0*(1 + x0**2 - 3*x0**4 + x0**6))/
     -  (x0**3*Sqrt((-1 + x0**2)**2)*
     -    (1 + x0**2)**2)
        gschwads_ll_xx(1,3,2,4)=
     -   0
        gschwads_ll_xx(1,3,3,3)=
     -   0
        gschwads_ll_xx(1,3,3,4)=
     -   0
        gschwads_ll_xx(1,3,4,4)=
     -   0

        gschwads_ll_x(1,4,1)=
     -   0
        gschwads_ll_x(1,4,2)=
     -   0
        gschwads_ll_x(1,4,3)=
     -   0
        gschwads_ll_x(1,4,4)=
     -   (2*M0*Sqrt((-1 + x0**2)**2))/
     -      (x0**2*(1 + x0**2))
        gschwads_ll_xx(1,4,1,1)=
     -   0
        gschwads_ll_xx(1,4,1,2)=
     -   0
        gschwads_ll_xx(1,4,1,3)=
     -   0
        gschwads_ll_xx(1,4,1,4)=
     -   0
        gschwads_ll_xx(1,4,2,2)=
     -   0
        gschwads_ll_xx(1,4,2,3)=
     -   0
        gschwads_ll_xx(1,4,2,4)=
     -   (-4*M0*(1 + x0**2 - 3*x0**4 + x0**6))/
     -  (x0**3*Sqrt((-1 + x0**2)**2)*
     -    (1 + x0**2)**2)
        gschwads_ll_xx(1,4,3,3)=
     -   0
        gschwads_ll_xx(1,4,3,4)=
     -   0
        gschwads_ll_xx(1,4,4,4)=
     -   0

        gschwads_ll_x(2,2,3)=
     -   0
        gschwads_ll_x(2,2,4)=
     -   0
        gschwads_ll_xx(2,2,1,3)=
     -   0
        gschwads_ll_xx(2,2,1,4)=
     -   0
        gschwads_ll_xx(2,2,2,3)=
     -   0
        gschwads_ll_xx(2,2,2,4)=
     -   0
        gschwads_ll_xx(2,2,3,3)=
     -   -16/(-1 + x0**2)**3 - 
     -  (12*M0*x0**2)/
     -   ((1 + x0**2)**2*Abs(x0)**5) + 
     -  (12*M0*x0**4)/
     -   ((1 + x0**2)**2*Abs(x0)**5) - 
     -  (16*M0*x0**2)/
     -   ((1 + x0**2)**3*Abs(x0)**3) + 
     -  (16*M0*x0**4)/
     -   ((1 + x0**2)**3*Abs(x0)**3) - 
     -  (8*M0*x0**2)/((1 + x0**2)**2*Abs(x0)**3)
        gschwads_ll_xx(2,2,3,4)=
     -   0
        gschwads_ll_xx(2,2,4,4)=
     -   -16/(-1 + x0**2)**3 - 
     -  (12*M0*x0**2)/
     -   ((1 + x0**2)**2*Abs(x0)**5) + 
     -  (12*M0*x0**4)/
     -   ((1 + x0**2)**2*Abs(x0)**5) - 
     -  (16*M0*x0**2)/
     -   ((1 + x0**2)**3*Abs(x0)**3) + 
     -  (16*M0*x0**4)/
     -   ((1 + x0**2)**3*Abs(x0)**3) - 
     -  (8*M0*x0**2)/((1 + x0**2)**2*Abs(x0)**3)

        gschwads_ll_x(2,3,1)=
     -   0
        gschwads_ll_x(2,3,2)=
     -   0
        gschwads_ll_x(2,3,3)=
     -   (-4*M0*x0*(-1 + x0**2))/
     -  ((1 + x0**2)**2*Abs(x0)**3)
        gschwads_ll_x(2,3,4)=
     -   0
        gschwads_ll_xx(2,3,1,1)=
     -   0
        gschwads_ll_xx(2,3,1,2)=
     -   0
        gschwads_ll_xx(2,3,1,3)=
     -   0
        gschwads_ll_xx(2,3,1,4)=
     -   0
        gschwads_ll_xx(2,3,2,2)=
     -   0
        gschwads_ll_xx(2,3,2,3)=
     -   (4*(-3*M0*x0**2 + 3*M0*x0**6 + 
     -      M0*Abs(x0)**2 - 
     -      6*M0*x0**2*Abs(x0)**2 + 
     -      M0*x0**4*Abs(x0)**2))/
     -  ((1 + x0**2)**3*Abs(x0)**5)
        gschwads_ll_xx(2,3,2,4)=
     -   0
        gschwads_ll_xx(2,3,3,3)=
     -   0
        gschwads_ll_xx(2,3,3,4)=
     -   0
        gschwads_ll_xx(2,3,4,4)=
     -   0

        gschwads_ll_x(2,4,1)=
     -   0
        gschwads_ll_x(2,4,2)=
     -   0
        gschwads_ll_x(2,4,3)=
     -   0
        gschwads_ll_x(2,4,4)=
     -   (-4*M0*x0*(-1 + x0**2))/
     -  ((1 + x0**2)**2*Abs(x0)**3)
        gschwads_ll_xx(2,4,1,1)=
     -   0
        gschwads_ll_xx(2,4,1,2)=
     -   0
        gschwads_ll_xx(2,4,1,3)=
     -   0
        gschwads_ll_xx(2,4,1,4)=
     -   0
        gschwads_ll_xx(2,4,2,2)=
     -   0
        gschwads_ll_xx(2,4,2,3)=
     -   0
        gschwads_ll_xx(2,4,2,4)=
     -   (4*(-3*M0*x0**2 + 3*M0*x0**6 + 
     -      M0*Abs(x0)**2 - 
     -      6*M0*x0**2*Abs(x0)**2 + 
     -      M0*x0**4*Abs(x0)**2))/
     -  ((1 + x0**2)**3*Abs(x0)**5)
        gschwads_ll_xx(2,4,3,3)=
     -   0
        gschwads_ll_xx(2,4,3,4)=
     -   0
        gschwads_ll_xx(2,4,4,4)=
     -   0

        gschwads_ll_x(3,3,1)=
     -   0
        gschwads_ll_x(3,3,2)=
     -   (-16*x0)/(-1 + x0**2)**3
        gschwads_ll_x(3,3,3)=
     -   0
        gschwads_ll_x(3,3,4)=
     -   0
        gschwads_ll_xx(3,3,1,1)=
     -   0
        gschwads_ll_xx(3,3,1,2)=
     -   0
        gschwads_ll_xx(3,3,1,3)=
     -   0
        gschwads_ll_xx(3,3,1,4)=
     -   0
        gschwads_ll_xx(3,3,2,2)=
     -   (16*(1 + 5*x0**2))/(-1 + x0**2)**4
        gschwads_ll_xx(3,3,2,3)=
     -   0
        gschwads_ll_xx(3,3,2,4)=
     -   0
        gschwads_ll_xx(3,3,3,3)=
     -   -16/(-1 + x0**2)**3 + 
     -  (8*M0)/((1 + x0**2)**2*Abs(x0)**3) - 
     -  (8*M0*x0**2)/((1 + x0**2)**2*Abs(x0)**3)
        gschwads_ll_xx(3,3,3,4)=
     -   0
        gschwads_ll_xx(3,3,4,4)=
     -   -16/(-1 + x0**2)**3

        gschwads_ll_x(3,4,1)=
     -   0
        gschwads_ll_x(3,4,2)=
     -   0
        gschwads_ll_x(3,4,3)=
     -   0
        gschwads_ll_x(3,4,4)=
     -   0
        gschwads_ll_xx(3,4,1,1)=
     -   0
        gschwads_ll_xx(3,4,1,2)=
     -   0
        gschwads_ll_xx(3,4,1,3)=
     -   0
        gschwads_ll_xx(3,4,1,4)=
     -   0
        gschwads_ll_xx(3,4,2,2)=
     -   0
        gschwads_ll_xx(3,4,2,3)=
     -   0
        gschwads_ll_xx(3,4,2,4)=
     -   0
        gschwads_ll_xx(3,4,3,3)=
     -   0
        gschwads_ll_xx(3,4,3,4)=
     -   (-4*M0*(-1 + x0**2))/
     -      ((1 + x0**2)**2*Abs(x0)**3)
        gschwads_ll_xx(3,4,4,4)=
     -   0

        gschwads_ll_x(4,4,1)=
     -   0
        gschwads_ll_x(4,4,2)=
     -   (-16*x0)/(-1 + x0**2)**3
        gschwads_ll_x(4,4,3)=
     -   0
        gschwads_ll_x(4,4,4)=
     -   0
        gschwads_ll_xx(4,4,1,1)=
     -   0
        gschwads_ll_xx(4,4,1,2)=
     -   0
        gschwads_ll_xx(4,4,1,3)=
     -   0
        gschwads_ll_xx(4,4,1,4)=
     -   0
        gschwads_ll_xx(4,4,2,2)=
     -   (16*(1 + 5*x0**2))/(-1 + x0**2)**4
        gschwads_ll_xx(4,4,2,3)=
     -   0
        gschwads_ll_xx(4,4,2,4)=
     -   0
        gschwads_ll_xx(4,4,3,3)=
     -   -16/(-1 + x0**2)**3
        gschwads_ll_xx(4,4,3,4)=
     -   0
        gschwads_ll_xx(4,4,4,4)=
     -    -16/(-1 + x0**2)**3 + 
     -  (8*M0)/((1 + x0**2)**2*Abs(x0)**3) - 
     -  (8*M0*x0**2)/((1 + x0**2)**2*Abs(x0)**3)


        do a=1,3
          do b=a+1,4
            gschwads_ll(b,a)=gschwads_ll(a,b)
            do c=1,4
              gschwads_ll_x(b,a,c)=gschwads_ll_x(a,b,c)
            end do
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                gschwads_ll_xx(a,b,c,d)=
     &             gschwads_ll_xx(min(a,b),max(a,b),min(c,d),max(c,d))
              end do
            end do
          end do
        end do

       end if !closes condition on y=z=0

        call calc_g0uu(gschwads_ll(1,1),gschwads_ll(1,2),
     &         gschwads_ll(1,3),gschwads_ll(1,4),
     &         gschwads_ll(2,2),gschwads_ll(2,3),
     &         gschwads_ll(2,4),
     &         gschwads_ll(3,3),gschwads_ll(3,4),
     &         gschwads_ll(4,4),
     &         gschwads_uu(1,1),gschwads_uu(1,2),
     &         gschwads_uu(1,3),gschwads_uu(1,4),
     &         gschwads_uu(2,2),gschwads_uu(2,3),
     &         gschwads_uu(2,4),
     &         gschwads_uu(3,3),gschwads_uu(3,4),
     &         gschwads_uu(4,4),detg0_schwads0)

        do a=1,3
          do b=a+1,4
            gschwads_uu(b,a)=gschwads_uu(a,b) 
          end do
        end do

        !set values for the FULL bulk scalar field value of the analytic solution

        phi1schwads=0
        phi1schwads_x(1)=0
        phi1schwads_x(2)=0
        phi1schwads_x(3)=0
        phi1schwads_x(4)=0
        phi1schwads_xx(1,1)=0
        phi1schwads_xx(1,2)=0
        phi1schwads_xx(1,3)=0
        phi1schwads_xx(1,4)=0
        phi1schwads_xx(2,2)=0
        phi1schwads_xx(2,3)=0
        phi1schwads_xx(2,4)=0
        phi1schwads_xx(3,3)=0
        phi1schwads_xx(3,4)=0

        do a=1,3
          do b=a+1,4
            phi1schwads_xx(b,a)=phi1schwads_xx(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            do c=1,4
              gschwads_uu_x(a,b,c)=
     &              -gschwads_ll_x(1,1,c)*gschwads_uu(a,1)
     &                  *gschwads_uu(b,1)
     &              -gschwads_ll_x(1,2,c)*(gschwads_uu(a,1)
     &                  *gschwads_uu(b,2)
     &                     +gschwads_uu(a,2)*gschwads_uu(b,1))
     &              -gschwads_ll_x(1,3,c)*(gschwads_uu(a,1)
     &                  *gschwads_uu(b,3)
     &                     +gschwads_uu(a,3)*gschwads_uu(b,1))
     &              -gschwads_ll_x(1,4,c)*(gschwads_uu(a,1)
     &                  *gschwads_uu(b,4)
     &                     +gschwads_uu(a,4)*gschwads_uu(b,1))
     &              -gschwads_ll_x(2,2,c)*gschwads_uu(a,2)
     &                  *gschwads_uu(b,2)
     &              -gschwads_ll_x(2,3,c)*(gschwads_uu(a,2)
     &                  *gschwads_uu(b,3)
     &                     +gschwads_uu(a,3)*gschwads_uu(b,2))
     &              -gschwads_ll_x(2,4,c)*(gschwads_uu(a,2)
     &                  *gschwads_uu(b,4)
     &                     +gschwads_uu(a,4)*gschwads_uu(b,2))
     &              -gschwads_ll_x(3,3,c)*gschwads_uu(a,3)
     &                  *gschwads_uu(b,3)
     &              -gschwads_ll_x(3,4,c)*(gschwads_uu(a,3)
     &                  *gschwads_uu(b,4)
     &                     +gschwads_uu(a,4)*gschwads_uu(b,3))
     &              -gschwads_ll_x(4,4,c)*gschwads_uu(a,4)
     &                  *gschwads_uu(b,4)
            end do
          end do
        end do

        ! give values to the Christoffel symbols
        do a=1,4
          do b=1,4
            do c=1,4
              gammaschwads_ull(a,b,c)=0
              do d=1,4
                gammaschwads_ull(a,b,c)=gammaschwads_ull(a,b,c)
     &                          +0.5d0*gschwads_uu(a,d)
     &                                *(gschwads_ll_x(c,d,b)
     &                                 -gschwads_ll_x(b,c,d)
     &                                 +gschwads_ll_x(d,b,c))
              end do
            end do
          end do
        end do

              ! calculate boxx^c at point i,j
              ! (boxschwadsx^c = -gschwads^ab gammaschwads^c_ab)
              do c=1,4
                boxschwadsx_u(c)=
     &           -( gammaschwads_ull(c,1,1)*gschwads_uu(1,1)+
     &              gammaschwads_ull(c,2,2)*gschwads_uu(2,2)+
     &              gammaschwads_ull(c,3,3)*gschwads_uu(3,3)+
     &              gammaschwads_ull(c,4,4)*gschwads_uu(4,4)+
     &            2*(gammaschwads_ull(c,1,2)*gschwads_uu(1,2)+
     &              gammaschwads_ull(c,1,3)*gschwads_uu(1,3)+
     &              gammaschwads_ull(c,1,4)*gschwads_uu(1,4)+
     &              gammaschwads_ull(c,2,3)*gschwads_uu(2,3)+
     &              gammaschwads_ull(c,2,4)*gschwads_uu(2,4)+
     &              gammaschwads_ull(c,3,4)*gschwads_uu(3,4)) )
              end do

              !compute Hschwads_l(a) in Cartesian coordinates
              ! (Hschwads_a = gschwads_ab boxschwadsx^b)
              do a=1,4
                Hschwads_l(a)=boxschwadsx_u(1)*gschwads_ll(a,1)+
     &                    boxschwadsx_u(2)*gschwads_ll(a,2)+
     &                    boxschwadsx_u(3)*gschwads_ll(a,3)+
     &                    boxschwadsx_u(4)*gschwads_ll(a,4)
              end do

        ! calculate Christoffel symbol derivatives at point i,j
        !(gamma^a_bc,e = 1/2 g^ad_,e(g_bd,c  + g_cd,b  - g_bc,d)
        !              +   1/2 g^ad(g_bd,ce + g_cd,be - g_bc,de))
        do a=1,4
          do b=1,4
            do c=1,4
              do e=1,4
                gammaschwads_ull_x(a,b,c,e)=0
                do d=1,4
                  gammaschwads_ull_x(a,b,c,e)=
     &             gammaschwads_ull_x(a,b,c,e)
     &              +0.5d0*gschwads_uu_x(a,d,e)*
     &                  (gschwads_ll_x(b,d,c)+
     &                    gschwads_ll_x(c,d,b)
     &                  -gschwads_ll_x(b,c,d))
     &              +0.5d0*gschwads_uu(a,d)*
     &                  (gschwads_ll_xx(b,d,c,e)+
     &                    gschwads_ll_xx(c,d,b,e)
     &                  -gschwads_ll_xx(b,c,d,e))
                end do
              end do
            end do
          end do
        end do

        ! calculate riemann tensor at point i,j
        !(R^a_bcd =gamma^a_bd,c - gamma^a_bc,d
        !          +gamma^a_ce gamma^e_bd - gamma^a_de gamma^e_bc)
        do a=1,4
          do b=1,4
            do c=1,4
              do d=1,4
                riemannschwads_ulll(a,b,c,d)=
     &                gammaschwads_ull_x(a,b,d,c)
     &               -gammaschwads_ull_x(a,b,c,d)
                do e=1,4
                   riemannschwads_ulll(a,b,c,d)=
     &                riemannschwads_ulll(a,b,c,d)
     &               +gammaschwads_ull(a,c,e)*
     &                      gammaschwads_ull(e,b,d)
     &               -gammaschwads_ull(a,d,e)*
     &                      gammaschwads_ull(e,b,c)
                end do
              end do
            end do
          end do
        end do

        ! calculate Ricci tensor at point i,j
        !(R_bd = R^a_bad)
        do b=1,4
          do d=1,4
            riccischwads_ll(b,d)=0
            do a=1,4
              riccischwads_ll(b,d)=riccischwads_ll(b,d)
     &            +riemannschwads_ulll(a,b,a,d)
            end do
          end do
        end do

        ! calculate raised Ricci tensor at point i,j
        !(R_a^b = R_ad g^db)
        do a=1,4
          do b=1,4
            riccischwads_lu(a,b)=0
            do d=1,4
              riccischwads_lu(a,b)=riccischwads_lu(a,b)
     &         +riccischwads_ll(a,d)*gschwads_uu(d,b)
            end do
          end do
        end do

        ! calculate Ricci scalar
        !(R = R_a^a)
        riccischwads=0
        do a=1,4
          riccischwads=riccischwads+riccischwads_lu(a,a)
        end do
  
        ! calculates Einstein tensor at point i,j
        !(G_ab = R_ab - 1/2 R g_ab)
        do a=1,4
          do b=1,4
            einsteinschwads_ll(a,b)=riccischwads_ll(a,b)
     &       -0.5d0*riccischwads*gschwads_ll(a,b)
          end do
        end do

        ! calculates stress-energy tensor at point i,j 
        !(T_ab = 2*phi1,a phi1,b - (phi1,c phi1,d) g^cd g_ab + ...)
        grad_phi1schwads_sq=0
        do a=1,4
          do b=1,4
            grad_phi1schwads_sq=grad_phi1schwads_sq
     &        +phi1schwads_x(a)*phi1schwads_x(b)*gschwads_uu(a,b)
          end do
        end do

        do a=1,4
          do b=1,4
            setschwads_ll(a,b)=
     &            phi1schwads_x(a)*phi1schwads_x(b)
     &           -gschwads_ll(a,b)*(grad_phi1schwads_sq/2)
          end do
        end do

!!!!!!!!!DEBUG!!!!
!        if ((abs(x0-(-1.0+5*dx)).lt.10.0d0**(-10))
!     -  .and.(abs(y0-(0.0d0)).lt.10.0d0**(-10))
!     -  .and.(abs(z0-(0.0d0)).lt.10.0d0**(-10))) then
!        write (*,*) "x0,y0,z0=",x0,y0,z0
!        write (*,*) "rho0,theta0,phi0=",rho0,theta0,phi0
!        do a=1,4
!          write (*,*) "a,Hschwads_l(a)="
!     -                   ,a,Hschwads_l(a)
!         do b=1,4
!             write (*,*) "a,b,gschwads_ll(a,b)="
!     -                   ,a,b,gschwads_ll(a,b)
!             write (*,*) "a,b,gschwads_uu(a,b)="
!     -                   ,a,b,gschwads_uu(a,b)
!          do c=1,4
!             write (*,*) "a,b,c,gschwads_ll_x(a,b,c)="
!     -                   ,a,b,c,gschwads_ll_x(a,b,c)
!             write (*,*) "a,b,c,gschwads_ll_sph_x(a,b,c)="
!     -                   ,a,b,c,gschwads_ll_sph_x(a,b,c)
!
!           do d=1,4
!             write (*,*) "a,b,c,d,gschwads_ll_xx(a,b,c,d)="
!     -                   ,a,b,c,d,gschwads_ll_xx(a,b,c,d)
!           end do
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!            write(*,*) "a,b,c,gammaschwads_ull(a,b,c)=",
!     -           a,b,c,gammaschwads_ull(a,b,c)
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!           do d=1,4
!            write(*,*) "a,b,c,d,gammaschwads_ull_x(a,b,c,d)=",
!     -           a,b,c,d,gammaschwads_ull_x(a,b,c,d)
!           end do
!          end do
!         end do
!        end do
!        do a=1,4
!         do b=1,4
!          do c=1,4
!           do d=1,4
!            write(*,*) "a,b,c,d,riemannschwads_ulll(a,b,c,d)=",
!     -           a,b,c,d,riemannschwads_ulll(a,b,c,d)
!           end do
!          end do
!         end do
!        end do
!
!         do a=1,4
!          do b=1,4
!            write (*,*) "a,b,
!     -       einsteinschwads_ll(a,b)+ Lambda* gschwads_ll(a,b)
!     -       setschwads_ll(a,b)="
!     -       ,a,b,
!     -        einsteinschwads_ll(a,b)-3*gschwads_ll(a,b),
!     -        setschwads_ll(a,b)
!          end do
!         end do
!         do b=1,4
!          do c=1,4
!            write (*,*) "b,c,riccischwads_ll(b,c)="
!     -       ,b,c,riccischwads_ll(b,c)
!          end do
!         end do
!         write (*,*) "riccischwads=",riccischwads
!          stop
!        end if
!!!!!!!!!!!

        
        return
        end
