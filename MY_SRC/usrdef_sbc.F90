MODULE usrdef_sbc
   !!======================================================================
   !!                     ***  MODULE  usrdef_sbc  ***
   !!
   !!                     ===  TCMIX configuration  ===
   !!
   !! User defined :   surface forcing of a user configuration
   !!======================================================================
   !! History :  4.0   ! 2016-03  (S. Flavoni, G. Madec)  user defined interface
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usrdef_sbc    : user defined surface bounday conditions in TCMIX case
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! Surface boundary condition: ocean fields
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distribued memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE lib_fortran    ! for glob_2Dsum
   !
   USE usrdef_nam, ONLY: lon0, lat0, rmw, msw, pc, poci, ln_move, tc_speed

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usrdef_sbc_oce       ! routine called in sbcmod module
   PUBLIC   usrdef_sbc_ice_tau   ! routine called by icestp.F90 for ice dynamics
   PUBLIC   usrdef_sbc_ice_flx   ! routine called by icestp.F90 for ice thermo

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usrdef_sbc_oce( kt, Kbb )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE usrdef_sbc  ***
      !!
      !! ** Purpose :   provide at each time-step the TCMIX surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :  The storm surface wind field is computed using the parametric cyclone model of Holland (1980).
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, sfx
      !!
      !! Reference : Holland, Greg J. 1980. “An Analytic Model of the Wind and Pressure Profiles in Hurricanes.”
      !!             Monthly Weather Review 108 (8): 1212–18. https://doi.org/10.1175/1520-0493(1980)108<1212:AAMOTW>2.0.CO;2.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      INTEGER, INTENT(in) ::   Kbb  ! ocean time index
      !!
      INTEGER  ::   ji, jj                ! dummy loop indices
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   zcoef, zmod           ! temporary variables
      REAL(wp) ::   ztime, ztimemax1, ztimemin1  ! 21th June, and 21th decem. if date0 = 1st january
      REAL(wp) ::   zcos_sais1
      REAL(wp) ::   zyydd                 ! number of days in one year

      REAL(wp) :: dlon_tc, lon_min, lon_max, lat_min, lat_max
      REAL(wp), DIMENSION(jpi,jpj) :: rdist, vr
      !!---------------------------------------------------------------------

      zyydd = REAL(nyear_len(1),wp)

      IF( kt == nit000 .AND. lwp ) THEN
        WRITE(numout,*) 'Time step:', kt, ' TC center (Lon Lat): ', lon0, lat0
      ENDIF

      ! Compute TC position at time step kt
      dlon_tc = 0.0_wp
      IF( ln_move ) THEN
          dlon_tc = ( tc_speed * rn_Dt ) / ( COS(lat0 * rad ) * 111000._wp)
          lon0 = lon0 - kt * dlon_tc     ! westward movement; constant latitude at its initial value lat0
          lon0 = MIN( MAX (lon0, -180.0_wp) , +180.0_wp)
          IF( kt == nit000 ) THEN
            WRITE(numout,*) 'TC moving zonally at lat =', lat0
            WRITE(numout,*) 'Degrees per timestep =', dlon_tc
          ELSE
            WRITE(numout,*) 'Time step:', kt, ' TC center (Lon Lat): ', lon0, lat0
          ENDIF
      ENDIF

      ! Compute bounds from glamt and gphit
      lon_min = MINVAL(glamt(:,:))
      lon_max = MAXVAL(glamt(:,:))
      lat_min = MINVAL(gphit(:,:))
      lat_max = MAXVAL(gphit(:,:))

      ! Check if point is within bounds
      IF ( .NOT. (lon0 >= lon_min .AND. lon0 <= lon_max) .AND. (lat0 >= lat_min .AND. lat0 <= lat_max) ) THEN
        CALL ctl_warn('userdef_sbc: Cyclone OUTSIDE the model domain, so STOP ! you are wasting your time !!')
      ENDIF

      CALL haversine_distance( kt, lat0, lon0, gphit, glamt, rdist )
      CALL holland_wnd( kt, rdist, lon0, lat0, rmw, msw, pc, poci, vr, utau, vtau )

      ! ---------------------------- !
      !       momentum fluxes        !
      ! ---------------------------- !
      ! module of wind stress and wind speed at T-point
      zcoef = 1. / ( zrhoa * zcdrag )
      DO_2D( 0, 0, 0, 0 )
         zmod = SQRT( utau(ji,jj) * utau(ji,jj) + vtau(ji,jj) * vtau(ji,jj) )
         taum(ji,jj) = zmod
         wndm(ji,jj) = SQRT( zmod * zcoef )
      END_2D

      ! ---------------------------- !
      !  heat and freshwater fluxes  !
      ! ---------------------------- !
      ! current day (in hours) since january the 1st of the current year
      ztime = REAL( kt ) * rn_Dt / (rmmss * rhhmm)   &       !  total incrementation (in hours)
         &      - (nyear  - 1) * rjjhh * zyydd               !  minus years since beginning of experiment (in hours)
      ztimemax1 = ((5.*30.)+21.)* 24.                      ! 21th june     at 24h in hours
      ztimemin1 = ztimemax1 + rjjhh * zyydd / 2            ! 21th december        in hours

      ! 1/2 period between 21th June and 21th December and between 21th July and 21th January
      zcos_sais1 = COS( (ztime - ztimemax1) / (ztimemin1 - ztimemax1) * rpi )

      ! set ZEROS to isolate TC mechanical forcing; except from 'qsr' required by the PISCES model for light availability.
      DO_2D( 0, 0, 0, 0 )
        sfx (ji,jj) = 0.0_wp
        qns(ji,jj)  = 0.0_wp
        emp(ji,jj)  = 0.0_wp
!         qsr(ji,jj)  = 0.0_wp
        qsr (ji,jj) =  230 * COS( 3.1415 * ( gphit(ji,jj) - 23.5 * zcos_sais1 ) / ( 0.9 * 180 ) )
      END_2D

      ! ---------------------------------- !
      !  control print at first time-step  !
      ! ---------------------------------- !
      IF( kt == nit000 .AND. lwp ) THEN
         WRITE(numout,*)
         WRITE(numout,*)'usrdef_sbc_oce : analytical surface fluxes for TCMIX configuration'
         WRITE(numout,*)'utau, vtau stress is computed from Holland (1980) parametric model which provide a symmetrical wind field around the cyclone centre'
         WRITE(numout,*)' ~~~~~~~~~~~ Zero heat fluxes except qsr required by the PISCES model for light availability (i.e.,  qns = emp = sfx = 0)'
      ENDIF

   END SUBROUTINE usrdef_sbc_oce

   SUBROUTINE holland_wnd(kt, rdist, lon, lat, rmw, msw, pc, poci, vr, utau, vtau)
   !!----------------------------------------------------------------------
   !!                   ***  ROUTINE holland_wnd  ***
   !!
   !! ** Purpose :    Compute tangential wind speed according to Holland (1980) model
   !!
   !! ** Method  : -  Reconstruct winds based on Holland (1980) parametric model
   !!              -  A common empirical formula for Cd as a function of wind speed based on Large and Pond (1981)
   !!
   !! References : Holland, Greg J. 1980. “An Analytic Model of the Wind and Pressure Profiles in Hurricanes.”
   !!             Monthly Weather Review 108 (8): 1212–18. https://doi.org/10.1175/1520-0493(1980)108<1212:AAMOTW>2.0.CO;2.
   !!             Large, W. G., and S. Pond, 1981: Open ocean momentum flux measurements in moderate to strong winds.
   !!             J. Phys. Oceanogr., 11, 324–336, https://doi.org/10.1175/1520-0485(1981)011<0324:OOMFMI>2.0.CO;2.
   !!----------------------------------------------------------------------

     ! Inputs
     INTEGER, INTENT(IN) :: kt
     REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)  :: rdist
     REAL(wp), INTENT(IN) :: lon, lat, rmw, msw, pc, poci

     ! Outputs
     REAL(wp), DIMENSION(jpi,jpj), INTENT(OUT) :: vr, utau, vtau

     ! Locals
     INTEGER :: ji, jj
     REAL(wp) :: f, b, term1, term2
     REAL(wp) :: dx, dy, theta, uwind, vwind, magv, Cd
     REAL(wp), PARAMETER :: rho_air = 1.15_wp ! [kg/m3]

     ! compute Holland model's shape parameter (B)
     b = rho_air * EXP(1._wp) * msw**2 / (poci - pc)

     DO_2D(0,0,0,0)

       IF (rdist(ji,jj) <= 0._wp) THEN
         vr(ji,jj) = 0._wp
         utau(ji,jj) = 0._wp
         vtau(ji,jj) = 0._wp
       ELSE
         f = 2._wp * 7.2921e-5_wp * SIN( lat * rad )
         term1 = ( b / rho_air ) * ( rmw / rdist(ji,jj) )**b * ( poci - pc ) * EXP( -( rmw / rdist(ji,jj))**b )
         term2 = ( rdist(ji,jj) * f / 2._wp )**2
         vr(ji,jj) = SQRT( term1 + term2 ) - rdist(ji,jj) * f / 2._wp

         ! Displacement from TC center
         dx = glamt(ji,jj) - lon
         dy = gphit(ji,jj) - lat
         theta = ATAN2(dy, dx)

         ! Wind components
         uwind = -vr(ji,jj) * SIN( theta )
         vwind =  vr(ji,jj) * COS( theta )
         magv = SQRT( uwind**2 + vwind**2 )

         ! Variable drag coefficient Cd as a function of wind speed (magv) based on Large and Pond (1981) empirical formula
         IF (magv < 11._wp) THEN
           Cd = 1.2e-3_wp
         ELSEIF (magv <= 25._wp) THEN
           Cd = (0.49_wp + 0.065_wp * magv) * 1.0e-3_wp
         ELSE
           Cd = 2.6e-3_wp
         END IF

         ! Wind stress
         utau(ji,jj) = rho_air * Cd * magv * uwind
         vtau(ji,jj) = rho_air * Cd * magv * vwind
         END IF

     END_2D

   END SUBROUTINE holland_wnd

   SUBROUTINE haversine_distance(kt, lat0, lon0, gphit, glamt, rdist)
   !!----------------------------------------------------------------------
   !!                   ***  ROUTINE haversine_distance  ***
   !!
   !! ** Purpose : Calculation of the great-circle distance from the TC centre (in km).
   !!
   !! ** Method  : Harvesine formula
   !!----------------------------------------------------------------------

     !IMPLICIT NONE
     INTEGER, INTENT(IN) :: kt
     REAL(wp), INTENT(IN) :: lat0, lon0
     REAL(wp), DIMENSION(jpi,jpj), INTENT(IN)  :: gphit, glamt
     REAL(wp), DIMENSION(jpi,jpj), INTENT(OUT) :: rdist

     INTEGER  :: ji, jj
     REAL(wp) :: dlat, dlon, a, c

     DO_2D(0,0,0,0)
       dlat = (gphit(ji,jj) - lat0) * rad
       dlon = (glamt(ji,jj) - lon0) * rad

       a = SIN( dlat / 2._wp )**2 + COS( lat0 * rad ) * COS( gphit(ji,jj) * rad ) * SIN( dlon / 2._wp )**2
       c = 2._wp * ATAN2( SQRT(a), SQRT(1._wp - a) )

       rdist(ji,jj) = ( ra / 1000.0_wp ) * c
     END_2D

   END SUBROUTINE haversine_distance

   SUBROUTINE usrdef_sbc_ice_tau( kt )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
   END SUBROUTINE usrdef_sbc_ice_tau

   SUBROUTINE usrdef_sbc_ice_flx( kt, phs, phi )
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phs    ! snow thickness
      REAL(wp), DIMENSION(:,:,:), INTENT(in)  ::   phi    ! ice thickness
   END SUBROUTINE usrdef_sbc_ice_flx

   !!======================================================================
END MODULE usrdef_sbc
