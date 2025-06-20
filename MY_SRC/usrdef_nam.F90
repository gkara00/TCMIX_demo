MODULE usrdef_nam
   !!======================================================================
   !!                     ***  MODULE usrdef_nam   ***
   !!
   !!                     ===  TCMIX configuration  ===
   !!
   !! User defined : set the domain characteristics of a user configuration
   !!======================================================================
   !! History :  4.0  ! 2016-03  (S. Flavoni, G. Madec)  Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_nam   : read user defined namelist and set global domain size
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_nam   ! called in nemogcm.F90 module

   !                              !!* namusr_def namelist *!!
   INTEGER, PUBLIC ::   nn_GYRE    ! 1/nn_GYRE = the resolution chosen in degrees and thus defining the horizontal domain size

   !!* Storm parameters for TCMIX config *!!
   REAL(wp), PUBLIC ::   lon0  = -68.1_wp   ! initial cyclone center's Longitude coordinate [Eastern degree]
   REAL(wp), PUBLIC ::   lat0  = 33.0_wp    ! initial cyclone center's Latitude coordinate [Northern degree]
   REAL(wp), PUBLIC ::   rmw   = 50.0_wp   ! radius of maximum sustained wind speed [km] ! Use at least 5–7 points across the RMW, so set Rm = 50–75 km as a realistic value.
   REAL(wp), PUBLIC ::   msw   = 50.0_wp         ! maximum sustained wind speed [m/s]

   REAL(wp), PUBLIC ::   pc    = 95000.0_wp      ! pressure at the centre of the storm [Pa]
   REAL(wp), PUBLIC ::   poci  = 101000.0_wp     ! pressure at outermost closed isobar of the storm [Pa]

   LOGICAL,  PUBLIC ::   ln_move = .false.       ! =T: Moving TC along the initial Latitude: only change of Longitude | =F stationary TC
   REAL(wp), PUBLIC ::   tc_speed    = 5.0_wp    ! Tropical cyclone moving speed [m/s]

   !!* Parameters for initial Temperature and Salinity profiles *!!
   ! Temperature
   REAL(wp), PUBLIC :: t_surf = 28.0_wp       ! Surface temp (degC)
   REAL(wp), PUBLIC :: t_deep = 4.0_wp        ! Deep ocean temp (degC)
   REAL(wp), PUBLIC :: z_mid_t  = 100.0_wp    ! Midpoint thermocline depth (m)
   REAL(wp), PUBLIC :: dz_tcl = 40.0_wp       ! Thermocline thickness scale (m)

   ! Salinity
   REAL(wp), PUBLIC :: s_surf  = 36.0_wp      ! Surface salinity (PSU)
   REAL(wp), PUBLIC :: s_deep  = 36.5_wp      ! Deep ocean salinity (PSU)
   REAL(wp), PUBLIC :: z_mid_s = 180.0_wp     ! Salinity gradient midpoint (m)
   REAL(wp), PUBLIC :: dz_scl  = 80.0_wp      ! Salinity thickness scale (m)

   !! * Substitutions
#  include "read_nml_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_nam( cd_cfg, kk_cfg, kpi, kpj, kpk, ldIperio, ldJperio, ldNFold, cdNFtype )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dom_nam  ***
      !!
      !! ** Purpose :   read user defined namelist and define the domain size
      !!
      !! ** Method  :   read in namusr_def containing all the user specific namelist parameter
      !!
      !!                Here TCMIX configuration based on GYRE reference config
      !!
      !! ** input   : - namusr_def namelist found in namelist_cfg
      !!----------------------------------------------------------------------
      CHARACTER(len=*), INTENT(out) ::   cd_cfg               ! configuration name
      INTEGER         , INTENT(out) ::   kk_cfg               ! configuration resolution
      INTEGER         , INTENT(out) ::   kpi, kpj, kpk        ! global domain sizes
      LOGICAL         , INTENT(out) ::   ldIperio, ldJperio   ! i- and j- periodicity
      LOGICAL         , INTENT(out) ::   ldNFold              ! North pole folding
      CHARACTER(len=1), INTENT(out) ::   cdNFtype             ! Folding type: T or F
#if defined key_agrif
      INTEGER :: ighost_n, ighost_s, ighost_w, ighost_e
#endif
      !
      INTEGER ::   ios   ! Local integer
      !!
      NAMELIST/namusr_def/ nn_GYRE, jpkglo  &
         &                 , lon0, lat0, rmw, msw, pc, poci &
         &                 , ln_move, tc_speed &
         &                 , t_surf, t_deep, z_mid_t, dz_tcl  &
         &                 , s_surf, s_deep, z_mid_s, dz_scl
      !!----------------------------------------------------------------------
      !
      READ_NML_(numnam_cfg,cfg,namusr_def,.TRUE.)
      IF(lwm)   WRITE( numond, namusr_def )
      !
      cd_cfg = 'GYRE'               ! name & resolution (not used)
#if defined key_agrif
      IF (.NOT.Agrif_root()) nn_GYRE = Agrif_parent(nn_GYRE) * Agrif_irhox()
#endif
      kk_cfg = nn_GYRE
      !
#if defined key_agrif
      IF( Agrif_Root() ) THEN
#endif
         kpi = 30 * nn_GYRE + 2       !
         kpj = 20 * nn_GYRE + 2
#if defined key_agrif
      ELSE                          ! Global Domain size: add nbghostcells + 1 "land" point on each side
         ! At this stage, child ghosts have not been set
         ighost_w = nbghostcells
         ighost_e = nbghostcells
         ighost_s = nbghostcells
         ighost_n = nbghostcells

         ! In case one sets zoom boundaries over domain edges:
         IF  ( Agrif_Ix() == 2 - Agrif_Parent(nbghostcells_x_w) ) ighost_w = 1
         IF  ( Agrif_Ix() + nbcellsx/AGRIF_Irhox() == Agrif_Parent(Ni0glo) - Agrif_Parent(nbghostcells_x_w) ) ighost_e = 1
         IF  ( Agrif_Iy() == 2 - Agrif_Parent(nbghostcells_y_s) ) ighost_s = 1
         IF  ( Agrif_Iy() + nbcellsy/AGRIF_Irhoy() == Agrif_Parent(Nj0glo) - Agrif_Parent(nbghostcells_y_s) ) ighost_n = 1
         kpi  = nbcellsx + ighost_w + ighost_e
         kpj  = nbcellsy + ighost_s + ighost_n
      ENDIF
#endif
      kpk = jpkglo
      !                             ! Set the lateral boundary condition of the global domain
      ldIperio = .FALSE.   ;   ldJperio = .FALSE.   ! GYRE configuration : closed domain
      ldNFold  = .FALSE.   ;   cdNFtype = '-'
      !
      !                             ! control print
      IF(lwp) THEN
         WRITE(numout,*) '   '
         WRITE(numout,*) 'usr_def_nam  : read the user defined namelist (namusr_def) in namelist_cfg'
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namusr_def : TCMIX case for idealized Tropical Cyclone (TC) wind forcing derived from Holland (1980) model'
         WRITE(numout,*) '      inverse resolution & implied domain size         nn_GYRE   = ', nn_GYRE
         WRITE(numout,*) '   '
         WRITE(numout,*) '~~~~~~~~~~~ Storm parameters'
         WRITE(numout,*) '      Initial TCs Longitude                            lon0     = ', lon0
         WRITE(numout,*) '      Initial TCs Latitude                             lat0     = ', lat0
         WRITE(numout,*) '      Radius of maximum sustained wind speed [km]      rmw      = ', rmw
         WRITE(numout,*) '      Maximum sustained wind speed [km]                msw      = ', msw
         WRITE(numout,*) '      Pressure at the TC centre [Pa]                   pc       = ', pc
         WRITE(numout,*) '      Pressure at outermost closed isobar of TC [Pa]   poci     = ', poci
         WRITE(numout,*) '      Moving [true] or stationary TC [false]           ln_move  = ', ln_move
         WRITE(numout,*) '      TC moving speed [m/s]                            tc_speed = ', tc_speed
         WRITE(numout,*) '   '
         WRITE(numout,*) '~~~~~~~~~~~ Initial T/S profile parameters'
         WRITE(numout,*) '      Surface temp (degC)                              t_surf  = ', t_surf
         WRITE(numout,*) '      Deep ocean temp (degC)                           t_deep  = ', t_deep
         WRITE(numout,*) '      Midpoint thermocline depth (m)                   z_mid_t = ', z_mid_t
         WRITE(numout,*) '      Thermocline thickness scale (m)                  dz_tcl  = ', dz_tcl
         WRITE(numout,*) '      Surface salinity (psu)                           s_surf  = ', s_surf
         WRITE(numout,*) '      Deep ocean salinity (psu)                        s_deep  = ', s_deep
         WRITE(numout,*) '      Salinity gradient midpoint (m)                   z_mid_s = ', z_mid_s
         WRITE(numout,*) '      Salinity thickness scale (m)                     dz_scl  = ', dz_scl

#if defined key_agrif
         IF( Agrif_Root() ) THEN
#endif
         WRITE(numout,*) '      Ni0glo = 30*nn_GYRE                              Ni0glo = ', kpi
         WRITE(numout,*) '      Nj0glo = 20*nn_GYRE                              Nj0glo = ', kpj
#if defined key_agrif
         ENDIF
#endif
         WRITE(numout,*) '      number of model levels                           jpkglo = ', kpk
         WRITE(numout,*) '   '
      ENDIF
      !
   END SUBROUTINE usr_def_nam

   !!======================================================================
END MODULE usrdef_nam
