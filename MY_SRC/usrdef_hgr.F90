MODULE usrdef_hgr
   !!======================================================================
   !!                     ***  MODULE usrdef_hgr   ***
   !!
   !!                     ===  TCMIX configuration  ===
   !!
   !! User defined :   mesh and Coriolis parameter of a user configuration
   !!======================================================================
   !! History :  4.0 ! 2016-03  (S. Flavoni)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   usr_def_hgr   : initialize the horizontal mesh
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE par_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE usrdef_nam     !
   !
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC   usr_def_hgr   ! called in domhgr.F90

   !! * Substitutions
#  include "do_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OCE 5.0, NEMO Consortium (2024)
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE usr_def_hgr( plamt , plamu , plamv  , plamf  ,   &   ! geographic position (required)
      &                    pphit , pphiu , pphiv  , pphif  ,   &   !
      &                    kff   , pff_f , pff_t  ,            &   ! Coriolis parameter  (if domain not on the sphere)
      &                    pe1t  , pe1u  , pe1v   , pe1f   ,   &   ! scale factors       (required)
      &                    pe2t  , pe2u  , pe2v   , pe2f   ,   &   !
      &                    ke1e2u_v      , pe1e2u , pe1e2v     )   ! u- & v-surfaces (if gridsize reduction is used in strait(s))
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE usr_def_hgr  ***
      !!
      !! ** Purpose :   user defined mesh and Coriolis parameter
      !!
      !! ** Method  :   set all intent(out) argument to a proper value
      !!
      !!                Here TCMIX based on GYRE configuration :
      !!          Rectangular mid-latitude domain
      !!          - with axes rotated by 45 degrees
      !!          - a constant horizontal resolution of 106 km
      !!          - on a beta-plane
      !!
      !! ** Action  : - define longitude & latitude of t-, u-, v- and f-points (in degrees)
      !!              - define coriolis parameter at f-point if the domain in not on the sphere (on beta-plane)
      !!              - define i- & j-scale factors at t-, u-, v- and f-points (in meters)
      !!              - define u- & v-surfaces (if gridsize reduction is used in some straits) (in m2)
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   plamt, plamu, plamv, plamf   ! longitude outputs                     [degrees]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pphit, pphiu, pphiv, pphif   ! latitude outputs                      [degrees]
      INTEGER                 , INTENT(out) ::   kff                          ! =1 Coriolis parameter computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pff_f, pff_t                 ! Coriolis factor at f-point                [1/s]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1t, pe1u, pe1v, pe1f       ! i-scale factors                             [m]
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe2t, pe2u, pe2v, pe2f       ! j-scale factors                             [m]
      INTEGER                 , INTENT(out) ::   ke1e2u_v                     ! =1 u- & v-surfaces computed here, =0 otherwise
      REAL(wp), DIMENSION(:,:), INTENT(out) ::   pe1e2u, pe1e2v               ! u- & v-surfaces (if reduction in strait)   [m2]
      !
      INTEGER  ::   ji, jj               ! dummy loop indices
      INTEGER  ::   ioffset, joffset     ! agrif indexes offset
      REAL(wp) ::   zlam1, zlam0, zcos_alpha, zim1 , zjm1 , ze1  , ze1deg, zf0 ! local scalars
      REAL(wp) ::   zphi1, zphi0, zsin_alpha, zim05, zjm05, zbeta, znorme      !   -      -
      !!-------------------------------------------------------------------------------
      !
      !     !==  beta-plane with regular grid-spacing and rotated domain ==!  (GYRE configuration)
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'usr_def_hgr : GYRE configuration (beta-plane with rotated regular grid-spacing)'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      !
      !
      !                       !==  grid point position  ==!
      !
      zlam1 = -85._wp                           ! position of F-point at (i,j) = (1,Njglo0-1)
      zphi1 =  29._wp
      !
      ze1 = 106000._wp / REAL( nn_GYRE , wp )   ! gridspacing in meters
      !
      zsin_alpha = - SQRT( 2._wp ) * 0.5_wp     ! angle: 45 degrees
      zcos_alpha =   SQRT( 2._wp ) * 0.5_wp
      ze1deg = ze1 / (ra * rad)
      ! longitude/latitude of the reference F-point at (i,j) = (1,1)
      zlam0 = zlam1 + zcos_alpha * ze1deg * REAL( Nj0glo - 2, wp )
      zphi0 = zphi1 + zsin_alpha * ze1deg * REAL( Nj0glo - 2, wp )

#if defined key_agrif
      IF (.NOT.Agrif_root()) THEN
         ! Retrieve Parent longitude/latitude south west corner:
         zlam0 = zlam1 + zcos_alpha * ze1deg * (REAL(Agrif_Parent(Nj0glo)-2, wp) * Agrif_Rhoy())
         zphi0 = zphi1 + zsin_alpha * ze1deg * (REAL(Agrif_Parent(Nj0glo)-2, wp) * Agrif_Rhoy())
         !
         ! Offset this point according to the child grid position:
         ioffset = (Agrif_Ix() + Agrif_Parent(nbghostcells_x_w) - 2) * Agrif_Irhox() - (nbghostcells_x_w - 1)
         joffset = (Agrif_Iy() + Agrif_Parent(nbghostcells_y_s) - 2) * Agrif_Irhoy() - (nbghostcells_y_s - 1)
         zlam0   = zlam0 + REAL(ioffset, wp)  * ze1deg * zcos_alpha + REAL(joffset, wp)  * ze1deg * zsin_alpha
         zphi0   = zphi0 - REAL(ioffset, wp)  * ze1deg * zsin_alpha + REAL(joffset, wp)  * ze1deg * zcos_alpha
      ENDIF
#endif
      !
      IF( lwp )   THEN
         WRITE(numout,*) 'ze1', ze1, 'cosalpha', zcos_alpha, 'sinalpha', zsin_alpha
         WRITE(numout,*) 'ze1deg', ze1deg, 'zlam0', zlam0, 'zphi0', zphi0
      ENDIF
      !
      DO_2D( nn_hls, nn_hls, nn_hls, nn_hls )
         zim1 = REAL( mig(ji,0), wp ) - 1.   ;   zim05 = REAL( mig(ji,0), wp ) - 1.5
         zjm1 = REAL( mjg(jj,0), wp ) - 1.   ;   zjm05 = REAL( mjg(jj,0), wp ) - 1.5
         !
         !glamt(i,j) longitude at T-point
         !gphit(i,j) latitude at T-point
         plamt(ji,jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
         pphit(ji,jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
         !
         !glamu(i,j) longitude at U-point
         !gphiu(i,j) latitude at U-point
         plamu(ji,jj) = zlam0 + zim1  * ze1deg * zcos_alpha + zjm05 * ze1deg * zsin_alpha
         pphiu(ji,jj) = zphi0 - zim1  * ze1deg * zsin_alpha + zjm05 * ze1deg * zcos_alpha
         !
         !glamv(i,j) longitude at V-point
         !gphiv(i,j) latitude at V-point
         plamv(ji,jj) = zlam0 + zim05 * ze1deg * zcos_alpha + zjm1  * ze1deg * zsin_alpha
         pphiv(ji,jj) = zphi0 - zim05 * ze1deg * zsin_alpha + zjm1  * ze1deg * zcos_alpha
         !
         !glamf(i,j) longitude at F-point
         !gphif(i,j) latitude at F-point
         plamf(ji,jj) = zlam0 + zim1  * ze1deg * zcos_alpha + zjm1  * ze1deg * zsin_alpha
         pphif(ji,jj) = zphi0 - zim1  * ze1deg * zsin_alpha + zjm1  * ze1deg * zcos_alpha
      END_2D
      !
      !                       !== Horizontal scale factors ==! (in meters)
      !
      !                                         ! constant grid spacing
      pe1t(:,:) =  ze1     ;      pe2t(:,:) = ze1
      pe1u(:,:) =  ze1     ;      pe2u(:,:) = ze1
      pe1v(:,:) =  ze1     ;      pe2v(:,:) = ze1
      pe1f(:,:) =  ze1     ;      pe2f(:,:) = ze1
      !
      !                                         ! NO reduction of grid size in some straits
      ke1e2u_v = 0                              !    ==>> u_ & v_surfaces will be computed in dom_ghr routine
      pe1e2u(:,:) = 0._wp                       !    CAUTION: set to zero to avoid error with some compilers that
      pe1e2v(:,:) = 0._wp                       !             require an initialization of INTENT(out) arguments
      !
      !
      !                       !==  Coriolis parameter  ==!
      kff = 1                                            !  indicate not to compute ff afterward
      !
      zbeta = 2. * omega * COS( rad * zphi1 ) / ra       ! beta at latitude zphi1
      !SF we overwrite zphi0 (south point in latitude) used just above to define pphif (value of zphi0=15.5190567531966)
      !SF for computation of Coriolis we keep the parameter of Hazeleger, W., and S. S. Drijfhout, JPO 1998.
      zphi0 = 15._wp                                     !  latitude of the most southern grid point
      zf0   = 2. * omega * SIN( rad * zphi0 )            !  compute f0 1st point south
      !
      pff_f(:,:) = ( zf0 + zbeta * ABS( pphif(:,:) - zphi0 ) * rad * ra ) ! f = f0 +beta* y ( y=0 at south)
      pff_t(:,:) = ( zf0 + zbeta * ABS( pphit(:,:) - zphi0 ) * rad * ra ) ! f = f0 +beta* y ( y=0 at south)
      !
      IF(lwp) WRITE(numout,*) '                           beta-plane used. beta = ', zbeta, ' 1/(s.m)'
      !
   END SUBROUTINE usr_def_hgr

   !!======================================================================
END MODULE usrdef_hgr
