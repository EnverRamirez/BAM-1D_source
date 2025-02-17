MODULE Micro_Ferrier
  IMPLICIT NONE
  PRIVATE
  INTEGER, PARAMETER :: r8=8
  REAL(kind=r8),PARAMETER:: con_pi     =3.1415926535897931_r8 ! pi
  REAL(kind=r8),PARAMETER:: con_cvap   =1.8460e+3_r8      ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cliq   =4.1855e+3_r8      ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap   =2.5000e+6_r8      ! lat heat H2O cond   (J/kg)
  REAL(kind=r8),PARAMETER:: con_rv     =4.6150e+2_r8      ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_ttp    =2.7316e+2_r8      ! temp at H2O 3pt     (K)
  REAL(kind=r8),PARAMETER:: con_psat   =6.1078e+2_r8      ! pres at H2O 3pt     (Pa)  
  real(kind=r8),parameter:: con_hfus   =3.3358e+5_r8      ! lat heat H2O fusion (J/kg)
  real(kind=r8),parameter:: con_cp     =1.0046e+3_r8      ! spec heat air @p    (J/kg/K)
  real(kind=r8),parameter:: con_rd     =2.8705e+2_r8      ! gas constant air    (J/kg/K)
  real(kind=r8),parameter:: con_eps    =con_rd/con_rv
  real(kind=r8),parameter:: con_csol   =2.1060e+3_r8      ! spec heat H2O ice   (J/kg/K)
  real(kind=r8),parameter:: con_epsm1  =con_rd/con_rv-1.0_r8
  real(kind=r8),parameter:: con_g      =9.80665e+0_r8     ! gravity           (m/s2)
  real(kind=r8),parameter:: con_dldt   =con_cvap-con_cliq
  real(kind=r8),parameter:: con_hsub   = con_hvap+con_hfus

  real(kind=r8),parameter:: con_rocp   =con_rd/con_cp

  !-----------------------------------------------------------------------
  !-------------- Local arrays & parameters in GSMDRIVE -----------------
  !-----------------------------------------------------------------------
  !
  !---- Comments on 14 March 2002
  !    * EPSQ=1.E-12 is the universal lower limit for specific humidity and 
  !      total condensate, and is consistent throughout the Eta code.
  !
  !     real(kind=r8), PARAMETER :: EPSQ=1.E-12,  RHOL=1000., T0C=273.15, 
  !REAL(kind=r8), PARAMETER :: EPSQ=1.0E-20_r8
  REAL(kind=r8), PARAMETER :: EPSQ=1.0E-12_r8
  REAL(KIND=r8),PARAMETER ::  RGRAV=1.0_r8/con_g
  REAL(kind=r8), PARAMETER :: RHOL=1000.0_r8
  REAL(kind=r8), PARAMETER :: T0C=273.15_r8
  REAL(kind=r8), PARAMETER :: T_ICE=-40.0_r8
  REAL(kind=r8), PARAMETER :: T_ICEK=T0C+T_ICE
  REAL(kind=r8), PARAMETER :: RRHOL=1.0_r8/RHOL
  REAL(kind=r8), PARAMETER :: EPSQ1=1.001_r8*EPSQ
  !    & T_ICE=-10._r8, T_ICEK=T0C+T_ICE, RRHOL=1._r8/RHOL, EPSQ1=1.001_r8*EPSQ
  !

  !
  !--- Common block of constants used in column microphysics
  !
  REAL(kind=r8), PRIVATE :: ABFR
  REAL(kind=r8), PRIVATE :: CBFR
  REAL(kind=r8), PRIVATE :: CIACW
  REAL(kind=r8), PRIVATE :: CIACR
  REAL(kind=r8), PRIVATE :: C_N0r0
  REAL(kind=r8), PRIVATE :: CN0r0
  REAL(kind=r8), PRIVATE :: CN0r_DMRmin
  REAL(kind=r8), PRIVATE :: CN0r_DMRmax
  REAL(kind=r8), PRIVATE :: CRACW
  REAL(kind=r8), PRIVATE :: CRAUT
  REAL(kind=r8), PRIVATE :: ESW0
  REAL(kind=r8), PRIVATE :: QAUTx
  REAL(kind=r8), PRIVATE :: RFmax
  REAL(kind=r8), PRIVATE :: RQR_DR1
  REAL(kind=r8), PRIVATE :: RQR_DR2
  REAL(kind=r8), PRIVATE :: RQR_DR3
  REAL(kind=r8), PRIVATE :: RQR_DRmin
  REAL(kind=r8), PRIVATE :: RQR_DRmax
  REAL(kind=r8), PRIVATE :: RR_DRmin
  REAL(kind=r8), PRIVATE :: RR_DR1
  REAL(kind=r8), PRIVATE :: RR_DR2
  REAL(kind=r8), PRIVATE :: RR_DR3
  REAL(kind=r8), PRIVATE :: RR_DRmax
  INTEGER      , PRIVATE :: mic_step
  !
  !--- Common block for lookup table used in calculating growth rates of
  !    nucleated ice crystals growing in water saturated conditions
  !--- Discretized growth rates of small ice crystals after their nucleation
  !     at 1 C intervals from -1 C to -35 C, based on calculations by Miller
  !     and Young (1979, JAS) after 600 s of growth.  Resultant growth rates
  !     are multiplied by physics time step in GSMCONST.
  !
  INTEGER, PRIVATE,PARAMETER :: MY_T1=1
  INTEGER, PRIVATE,PARAMETER :: MY_T2=35
  REAL(kind=r8),PRIVATE       :: MY_GROWTH(MY_T1:MY_T2)
  !
  !--- Parameters for ice lookup tables, which establish the range of mean ice
  !    particle diameters; from a minimum mean diameter of 0.05 mm (DMImin) to a
  !    maximum mean diameter of 1.00 mm (DMImax).  The tables store solutions
  !    at 1 micron intervals (DelDMI) of mean ice particle diameter.
  !
  REAL(kind=r8), PRIVATE,PARAMETER :: DMImin=.05e-3_r8
  REAL(kind=r8), PRIVATE,PARAMETER :: DMImax=1.e-3_r8
  REAL(kind=r8), PRIVATE,PARAMETER :: XMImin=1.e6_r8*DMImin
  REAL(kind=r8), PRIVATE,PARAMETER :: XMImax=1.e6_r8*DMImax
  REAL(kind=r8), PRIVATE,PARAMETER :: DelDMI=1.e-6_r8
  INTEGER      , PRIVATE,PARAMETER :: MDImin=INT(XMImin)
  INTEGER      , PRIVATE,PARAMETER :: MDImax=INT(XMImax)
  !
  !--- Various ice lookup tables
  !
  REAL(kind=r8), PRIVATE ::    ACCRI (MDImin:MDImax)
  REAL(kind=r8), PRIVATE ::    MASSI (MDImin:MDImax)
  REAL(kind=r8), PRIVATE ::    SDENS (MDImin:MDImax)
  REAL(kind=r8), PRIVATE ::    VSNOWI(MDImin:MDImax)
  REAL(kind=r8), PRIVATE ::    VENTI1(MDImin:MDImax)
  REAL(kind=r8), PRIVATE ::    VENTI2(MDImin:MDImax)
  !
  !--- Mean rain drop diameters varying from 50 microns (0.05 mm) to 450 microns
  !      (0.45 mm), assuming an exponential size distribution.
  !
  REAL(kind=r8), PRIVATE,PARAMETER :: DMRmin=.05e-3_r8
  REAL(kind=r8), PRIVATE,PARAMETER :: DMRmax=.45e-3_r8
  REAL(kind=r8), PRIVATE,PARAMETER :: XMRmin=1.e6_r8*DMRmin
  REAL(kind=r8), PRIVATE,PARAMETER :: XMRmax=1.e6_r8*DMRmax
  REAL(kind=r8), PRIVATE,PARAMETER :: DelDMR=1.e-6_r8
  REAL(kind=r8), PRIVATE,PARAMETER :: NLImin=100.0_r8
  !    &,                          NLImin=100., NLImax=20.E3
  INTEGER      , PRIVATE,PARAMETER :: MDRmin=INT(XMRmin)
  INTEGER      , PRIVATE,PARAMETER :: MDRmax=INT(XMRmax)
  !
  !--- Factor of 1.5 for RECImin, RESNOWmin, & RERAINmin accounts for
  !    integrating exponential distributions for effective radius
  !    (i.e., the r**3/r**2 moments).
  !
  !     INTEGER, PRIVATE, PARAMETER :: INDEXSmin=300
  !!    INTEGER, PRIVATE, PARAMETER :: INDEXSmin=200
  INTEGER      , PRIVATE, PARAMETER :: INDEXSmin=100
  REAL(kind=r8), PRIVATE, PARAMETER :: RERAINmin=1.5_r8*XMRmin
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=8.0
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=7.5
  REAL(kind=r8), PRIVATE, PARAMETER :: RECImin=1.5_r8*XMImin
  REAL(kind=r8), PRIVATE, PARAMETER :: RESNOWmin=1.5_r8*INDEXSmin
  REAL(kind=r8), PRIVATE, PARAMETER :: RECWmin=10.0_r8
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=15.
  !    &, RECImin=1.5*XMImin, RESNOWmin=1.5*INDEXSmin, RECWmin=5.

  !
  !--- Various rain lookup tables
  !--- Rain lookup tables for mean rain drop diameters from DMRmin to DMRmax,
  !      assuming exponential size distributions for the rain drops
  !
  REAL(kind=r8), PRIVATE :: ACCRR (MDRmin:MDRmax)
  REAL(kind=r8), PRIVATE :: MASSR (MDRmin:MDRmax)
  REAL(kind=r8), PRIVATE :: RRATE (MDRmin:MDRmax)
  REAL(kind=r8), PRIVATE :: VRAIN (MDRmin:MDRmax)
  REAL(kind=r8), PRIVATE :: VENTR1(MDRmin:MDRmax)
  REAL(kind=r8), PRIVATE :: VENTR2(MDRmin:MDRmax)
  !
  !--- Common block for riming tables
  !--- VEL_RF - velocity increase of rimed particles as functions of crude
  !      particle size categories (at 0.1 mm intervals of mean ice particle
  !      sizes) and rime factor (different values of Rime Factor of 1.1**N,
  !      where N=0 to Nrime).
  !
  INTEGER      , PRIVATE,PARAMETER :: Nrime=40
  REAL(kind=r8), PRIVATE           :: VEL_RF(2:9,0:Nrime)
  !
  !--- The following variables are for microphysical statistics
  !
  INTEGER, PARAMETER :: ITLO=-60
  INTEGER, PARAMETER :: ITHI=40
  INTEGER            :: NSTATS(ITLO:ITHI,4)
  REAL(kind=r8)      :: QMAX(ITLO:ITHI,5)
  REAL(kind=r8)      :: QTOT(ITLO:ITHI,22)
  !
  REAL(kind=r8), PRIVATE,  PARAMETER ::   T_ICE_init=-15.0_r8
  !    &  T_ICE=-10., T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.                      !- Ver2
  !     &  T_ICE=-40., T_ICE_init=-15.     !- Ver2
  !    &  T_ICE=-30., T_ICE_init=-5.      !- Ver2
  !
  !     Some other miscellaneous parameters
  !
  REAL(kind=r8), PRIVATE, PARAMETER :: Thom=T_ICE
  REAL(kind=r8), PRIVATE, PARAMETER :: TNW=50.0_r8
  REAL(kind=r8), PRIVATE, PARAMETER :: TOLER=1.0E-20_r8
  !     real(kind=r8), PRIVATE, PARAMETER :: Thom=T_ICE, TNW=50._r8, TOLER=5.E-7_r8
  !     real(kind=r8), PRIVATE, PARAMETER :: Thom=-35._r8, TNW=50._r8, TOLER=5.E-7_r8

  ! Assume fixed cloud ice effective radius
  REAL(kind=r8), PRIVATE, PARAMETER :: RECICE=RECImin
  !     &, EPSQ=1.0E-20                                                    &
  !    &, EPSQ=1.E-12                                                     &
  REAL(kind=r8), PRIVATE, PARAMETER :: FLG0P1=0.1_r8
  REAL(kind=r8), PRIVATE, PARAMETER :: FLG0P2=0.2_r8
  REAL(kind=r8), PRIVATE, PARAMETER :: FLG1P0=1.0_r8
  !
  !
  INTEGER       , PARAMETER   :: ntrac=3
  INTEGER,PARAMETER :: nxpvsl=7501
  REAL(r8)          :: tbpvsl(nxpvsl)
  REAL(r8)          :: c1xpvsl
  REAL(r8)          :: c2xpvsl
  integer,parameter:: nxpvsi=7501
  REAL(r8)          :: tbpvsi(nxpvsi)
  REAL(r8)          :: c1xpvsi
  REAL(r8)          :: c2xpvsi


  INTEGER           :: ncwp(2)!     ncwp      - integer, range of droplet number concentrations for    !
!                         Ferrier microphysics                     2    !

  REAL(r8)          :: flgmin(2)!     flgmin   - real, range of  minimum large ice fraction for         !
!                         Ferrier microphys                        2    !

  REAL(r8)          :: crtrh(3)!     crtrh    - real, critical relative humidity at the surface, PBL   !
                               !                      top and at the top of the atmosphere        3    !
  integer :: nlons         !     nlons(im)    - integer, number of total grid points in a latitude     !
                               !                         circle through a point                   im   !
!     lonf,
  integer ::  latg  !latg- integer, number of lon/lat points                 1    !

  !    & dxmax=ln(1.0_r8/(5000.0_r8*2500.0_r8)),  dxmin=ln(1.0_r8/(192.0_r8*94.0_r8)
  !    & dxmax=-16.118095651_r8,dxmin=-9.800790154_r8,dxinv=1.0_r8/(dxmax-dxmin)
  real (kind=r8) :: dxmaxs!=-16.118095651_r8
  real (kind=r8) :: dxmins!= -9.800790154_r8
  real (kind=r8) :: dxinvs!=  1.0/(dxmaxs-dxmins)
  real (kind=r8) :: dxmax != dxmaxs
  real (kind=r8) :: dxmin != dxmins
  real (kind=r8) :: dxinv != dxinvs
  real(kind=r8), parameter :: rhc_max = 0.9999_r8
  real (kind=r8), ALLOCATABLE  :: sik  (:)
  real (kind=r8), ALLOCATABLE  :: sikp1(:)
  real (kind=r8), ALLOCATABLE  :: si   (:)

  real (kind=r8), ALLOCATABLE  ::  slk(:)
  real (kind=r8), ALLOCATABLE  ::  sl(:)
  PUBLIC :: Init_Micro_Ferrier
  PUBLIC :: RunMicro_FERRIER
CONTAINS
  SUBROUTINE Init_Micro_Ferrier (dt,iMax,jMax,ibmax,kMax,jbMax,F_ICE_PHY ,&
           F_RAIN_PHY ,F_RIMEF_PHY,si_in,sl_in)
    IMPLICIT NONE
    REAL (kind=r8), INTENT(IN ) :: dt
    INTEGER       , INTENT(IN ) :: iMax
    INTEGER       , INTENT(IN ) :: jMax
    INTEGER       , INTENT(IN ) :: ibMax
    INTEGER       , INTENT(IN ) :: kMax
    INTEGER       , INTENT(IN ) :: jbMax
    REAL (kind=r8), INTENT(OUT) :: F_ICE_PHY   (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RAIN_PHY  (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RIMEF_PHY (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(IN ) :: si_in (kMax+1)
    REAL (kind=r8), INTENT(IN ) :: sl_in (kMax)
    INTEGER  :: k
    ALLOCATE( sik  (kMax+1))
    ALLOCATE( sikp1(kMax+1))
    ALLOCATE( si   (kMax+1))
    ALLOCATE( slk  (kMax))
    ALLOCATE( sl   (kMax))
    
    crtrh(1)         = 0.95_r8 !critical relative humidity at the surface
    crtrh(2)         = 0.95_r8 !critical relative humidity at the PBL top
    crtrh(3)         = 0.95_r8 !critical relative humidity at the top of the atmosphere   

    flgmin(1)        = 0.20_r8 !range of  minimum large ice fraction for Ferrier microphys  maximo no polo
    flgmin(2)        = 0.20_r8 !range of  minimum large ice fraction for Ferrier microphys  maximo no equador

!     ncwp(1)           = 75
    ncwp(1)          = 50 !range of droplet number concentrations for Ferrier microphysics maximo no polo
    ncwp(2)          = 150!range of droplet number concentrations for Ferrier microphysics maximo no equador

    si=si_in
    sl=sl_in
    nlons=iMax
    latg =jMax
    dxmaxs=-16.118095651_r8!log(1.0_r8/(REAL(iMax*20)*REAL(jMax*20)))
    dxmins= -9.800790154_r8!log(1.0_r8/(REAL(iMax/(20))*REAL(jMax/(20))))
    dxinvs=  1.0_r8/(dxmaxs-dxmins)
    dxmax = dxmaxs
    dxmin = dxmins
    dxinv = dxinvs

!  move water from vapor to liquid should the liquid amount be negative
      DO k=1,kMax+1
          sik  (k)   = si(k) ** con_rocp
          sikp1(k)   = si(k) ** (1.0_r8 + con_rocp)
      END DO
      DO k=1,kMax
         slk(k)   = (sikp1(k)-sikp1(k+1))/((1.0_r8 + con_rocp) * (si(k) - si(k+1)))
      END DO


    CALL gpvsl()
    CALL gpvsi()
    CALL INIT_MICRO(DT,ibmax,kMax,jbMax,F_ICE_PHY , & 
         F_RAIN_PHY ,F_RIMEF_PHY)
  END SUBROUTINE Init_Micro_Ferrier
  !

  !-----------------------------------------------------------------------
  !-----------------------------------------------------------------------
  SUBROUTINE RunMicro_FERRIER (&
       nCols       , &!INTEGER      , INTENT(IN)     :: nCols
       kMax        , &!INTEGER      , INTENT(IN)     :: kMax
       DT          , &!REAL(KIND=r8), INTENT(IN)     :: DT
       RHO         , &!REAL(KIND=r8), INTENT(IN),    :: RHO       (ims:ime, kms:kme, jms:jme)
       gt          , &!REAL(KIND=r8), INTENT(INOUT), :: gt         (ims:ime, kms:kme, jms:jme)
       qv          , &!REAL(KIND=r8), INTENT(INOUT), :: qv         (ims:ime, kms:kme, jms:jme)
       QC          , &!REAL(KIND=r8), INTENT(INOUT), :: qc         (ims:ime, kms:kme, jms:jme)
       QI          , &!REAL(KIND=r8), INTENT(INOUT), :: qi         (ims:ime, kms:kme, jms:jme)
       QR          , &!REAL(KIND=r8), INTENT(INOUT), :: qr         (ims:ime, kms:kme, jms:jme)
       F_ICE_PHY   , &!REAL(KIND=r8), INTENT(INOUT), :: F_ICE_PHY  (ims:ime, kms:kme, jms:jme)
       F_RAIN_PHY  , &!REAL(KIND=r8), INTENT(INOUT), :: F_RAIN_PHY (ims:ime, kms:kme, jms:jme)
       F_RIMEF_PHY , &!REAL(KIND=r8), INTENT(INOUT), :: F_RIMEF_PHY(ims:ime, kms:kme, jms:jme)
       RAINNCV     , &!REAL(KIND=r8), INTENT(INOUT), :: RAINNCV    (ims:ime,          jms:jme) 
       SNOWNCV     , &!REAL(KIND=r8), INTENT(INOUT), :: SNOWNCV    (ims:ime,          jms:jme) 
       pgr          , &
       colrad        )
    !-----------------------------------------------------------------------
    IMPLICIT NONE
    !-----------------------------------------------------------------------

    INTEGER         ,INTENT(IN   ) :: nCols
    INTEGER         ,INTENT(IN   ) :: kMax
    REAL(KIND=r8)   ,INTENT(IN   ) :: DT
    REAL(KIND=r8)   ,INTENT(IN   ) :: RHO        (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: gt         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qv         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qc         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qi         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: qr         (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_ICE_PHY  (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RAIN_PHY (1:nCols,1:kMax)
    REAL(KIND=r8)   ,INTENT(INOUT) :: F_RIMEF_PHY(1:nCols,1:kMax)
!    REAL(KIND=r8)   ,INTENT(INOUT) :: RAINNC     (1:nCols)
    REAL(KIND=r8)   ,INTENT(OUT  ) :: RAINNCV    (1:nCols)
    REAL(KIND=r8)   ,INTENT(OUT  ) :: SNOWNCV    (1:nCols)    
    REAL(KIND=r8)   ,INTENT(IN   ) :: pgr(nCols)

    REAL(KIND=r8)   ,INTENT(IN   ) :: colrad (nCols)

    !-----------------------------------------------------------------------
    !     LOCAL VARS
    !-----------------------------------------------------------------------

    !     NSTATS,QMAX,QTOT are diagnostic vars

    !     SOME VARS WILL BE USED FOR DATA ASSIMILATION (DON'T NEED THEM NOW). 
    !     THEY ARE TREATED AS LOCAL VARS, BUT WILL BECOME STATE VARS IN THE 
    !     FUTURE. SO, WE DECLARED THEM AS MEMORY SIZES FOR THE FUTURE USE

    !     TLATGS_PHY,TRAIN_PHY,APREC,PREC,ACPREC,SR are not directly related 
    !     the microphysics scheme. Instead, they will be used by Eta precip 
    !     assimilation.
    REAL(KIND=r8)  :: SR         (1:nCols)
    REAL(KIND=r8)  :: TLATGS_PHY(1:nCols,1:kMax)
    REAL(KIND=r8)  :: TRAIN_PHY (1:nCols,1:kMax)
    REAL(KIND=r8)  :: APREC     (1:nCols)
    REAL(KIND=r8)  :: ASNOW     (1:nCols) 

    REAL(KIND=r8)  :: PREC      (1:nCols)
    REAL(KIND=r8)  :: ACPREC    (1:nCols)
    REAL(KIND=r8)  :: t_phy     (1:nCols,1:kMax)
    REAL(KIND=r8)  :: CWM_PHY   (1:nCols,1:kMax)
    REAL(KIND=r8)  :: WC
    INTEGER :: I,K

    !
    !-----------------------------------------------------------------------
    !**********************************************************************
    !-----------------------------------------------------------------------
    !
    !
    DO k = 1,kMax
       DO i = 1,nCols
          t_phy(i,k) = gt(i,k)
          !qr(i,k) = 0.0_r8!PK
       ENDDO
    ENDDO

    !     initial diagnostic variables and data assimilation vars
    !     (will need to delete this part in the future)


    ! initial data assimilation vars (will need to delete this part in the future)

    DO k = 1,kMax
       DO i = 1,nCols
          TLATGS_PHY (i,k)=0.0_r8
          TRAIN_PHY  (i,k)=0.0_r8
       ENDDO
    ENDDO
    DO i = 1,nCols
       ACPREC(i)=0.0_r8
       APREC (i)=0.0_r8 
       ASNOW (i)=0.0_r8
       PREC  (i)=0.0_r8
       SR    (i)=0.0_r8
    ENDDO

    !-- 6/11/2010: Update CWM_PHY, F_ice, F_rain arrays

    DO k = 1,kMax
       DO i = 1,nCols
          CWM_PHY(I,K)=QC(I,K)+QR(I,K)+QI(I,K)
          IF (QI(I,K) <= EPSQ) THEN
             F_ICE_PHY(I,K)=0.0_r8
             IF (T_PHY(I,K) < T_ICEK) F_ICE_PHY(I,K)=1.0_r8
          ELSE
             F_ICE_PHY(I,K)=MAX( 0.0_r8, MIN(1.0_r8, QI(I,K)/CWM_PHY(I,K) ) )
          ENDIF
          IF (QR(I,K) <= EPSQ) THEN
             F_RAIN_PHY(I,K)=0.0_r8
          ELSE
             F_RAIN_PHY(I,K)=QR(I,K)/(QR(I,K)+QC(I,K))
          ENDIF
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    CALL GSMDRIVE( &
      DT         , &!real (kind=r8) :: DT
      APREC      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: APREC       (1:nCols)
      ASNOW      , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ASNOW       (1:nCols)
      PREC       , &!REAL(KIND=r8)   ,INTENT(INOUT) :: PREC        (1:nCols)
      ACPREC     , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC      (1:nCols)
      SR         , &!real (kind=r8) :: SR     (nCols)
      RHO        , &!real (kind=r8) :: RHO    (nCols,kMax)      
      CWM_PHY    , &!real (kind=r8) :: CWM_PHY  (nCols,kMax)
      t_phy      , &!real (kind=r8) :: t_phy    (nCols,kMax)
      qv         , &!real (kind=r8) :: qv    (nCols,kMax)
      qc         , &!(1:nCols,1:kMax)
      qi         , &!(1:nCols,1:kMax)
      qr         , &!(1:nCols,1:kMax)
      F_ICE_PHY  , &!real (kind=r8) :: F_ICE_PHY  (nCols,kMax)
      F_RAIN_PHY , &!real (kind=r8) :: F_RAIN_PHY (nCols,kMax)
      F_RIMEF_PHY, &!real (kind=r8) :: F_RimeF(nCols,kMax)
      TLATGS_PHY , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
      TRAIN_PHY  , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
      colrad     , &!real (kind=r8) :: colrad (nCols)
      pgr        , &!real (kind=r8) :: pgr    (nCols)
      nCols      , &!integer, INTENT(IN   ) ::  nCols
      kMax         )!integer, INTENT(IN   ) ::  kMax

!
    !-----------------------------------------------------------------------

    DO k = 1,kMax
       DO i = 1,nCols
          gt(i,k) = t_phy    (i,k)
          WC      = CWM_PHY(I,K)
          !QI(I,K) = 0.0_r8
          !QR(I,K) = 0.0_r8
          !QC(I,K) = 0.0_r8
          !IF(F_ICE_PHY(I,K)>=1.0_r8)THEN
          !   QI(I,K)=WC
          !ELSEIF(F_ICE_PHY(I,K)<=0.0_r8)THEN
          !   QC(I,K)=WC
          !ELSE
          !   QI(I,K)=F_ICE_PHY(I,K)*WC
          !   QC(I,K)=WC-QI(I,K)
          !ENDIF
          !
          !IF(QC(I,K)>0.0_r8.AND.F_RAIN_PHY(I,K)>0.0_r8)THEN
          !   IF(F_RAIN_PHY(I,K).GE.1.0_r8)THEN
          !      QR(I,K)=QC(I,K)
          !      QC(I,K)=0.0_r8
          !   ELSE
          !      QR(I,K)=F_RAIN_PHY(I,K)*QC(I,K)
          !      QC(I,K)=QC(I,K)-QR(I,K)
          !   ENDIF
          !ENDIF
       ENDDO
    ENDDO
    ! 
    ! update rain (from m to mm)

    DO i=1,nCols
       RAINNCV(i)=APREC(i)*0.5_r8!*1000.0_r8
       SNOWNCV(i)=ASNOW(i)*0.5_r8!*1000.0_r8
    ENDDO
    !
    !-----------------------------------------------------------------------


  END SUBROUTINE RunMicro_FERRIER
!@PROCESS NOEXTCHK
!
!--- The 1st line is an inlined compiler directive that turns off -qextchk
!    during compilation, even if it's specified as a compiler option in the
!    makefile (Tuccillo, personal communication;  Ferrier, Feb '02).
!
!###############################################################################
!---------------------- Driver of the new microphysics -------------------------
!###############################################################################
!
      SUBROUTINE GSMDRIVE( &
      DT            , &!real (kind=r8) :: DT
      APREC         , &!real (kind=r8) :: APREC  (nCols)
      ASNOW         , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ASNOW       (1:nCols)
      PREC          , &!REAL(KIND=r8)   ,INTENT(INOUT) :: PREC        (1:nCols)
      ACPREC        , &!REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC      (1:nCols)
      SR            , &!real (kind=r8) :: SR     (nCols)
      RHO           , &!real (kind=r8) :: RHO    (nCols,kMax)
      CWM_PHY       , &!real (kind=r8) :: CWM_PHY   (nCols,kMax)
      T_PHY         , &!real (kind=r8) :: T_PHY    (nCols,kMax)
      Q_PHY         , &!real (kind=r8) :: Q_PHY   (nCols,kMax)
      qc            , &!(1:nCols,1:kMax)
      qi            , &!(1:nCols,1:kMax)
      qr            , &!(1:nCols,1:kMax)
      F_ICE_PHY     , &!real (kind=r8) :: F_ICE_PHY  (nCols,kMax)
      F_RAIN_PHY    , &!real (kind=r8) :: F_RAIN_PHY (nCols,kMax)
      F_RimeF_PHY   , &!real (kind=r8) :: F_RimeF_PHY(nCols,kMax)
      TLATGS_PHY    , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
      TRAIN_PHY     , &!REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
      colrad        , &!real (kind=r8) :: colrad (nCols)
      pgr           , &!real (kind=r8) :: pgr    (nCols)
      nCols         , &!integer, INTENT(IN   ) ::  nCols
      kMax            )!integer, INTENT(IN   ) ::  kMax

!
!-------------------------------------------------------------------------------
!----- NOTE:  Code is currently set up w/o threading!  
!-------------------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:  Grid-scale microphysical processes - condensation & precipitation
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: February 2001
!  2001-04-xx   Ferrier     - Beta-tested version
!  2001-05-21   Ferrier     - Added gradual latent heating to remove external waves
!  2001-05-30   Ferrier     - Changed default to uniform maritime conditions for testing
!  2001-11-09   Moorthi     - Modified for Global Spectral Model
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Merges original GSCOND & PRECPD subroutines.   
!   * Code has been substantially streamlined and restructured.
!   * Exchange between water vapor & small cloud condensate is calculated using
!     the original Asai (1965, J. Japan) algorithm.  See also references to
!     Yau and Austin (1979, JAS), Rutledge and Hobbs (1983, JAS), and Tao et al.
!     (1989, MWR).  This algorithm replaces the Sundqvist et al. (1989, MWR)
!     parameterization.  
!-------------------------------------------------------------------------------
! Prior PROGRAM HISTORY LOG:
!
! *** Heritage as Subroutine GSCOND:
!   94-~??  ZHAO         - ORIGINATOR
!   95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   95-03-28  BLACK      - ADDED EXTERNAL EDGE
!   98-11-02  BLACK      - MODIFIED FOR DISTRIBUTED MEMORY
!
! *** Heritage as Subroutine PRECPD:
!   94-~??  ZHAO       - ORIGINATOR
!   95-03-25  BLACK      - CONVERSION FROM 1-D TO 2-D IN HORIZONTAL
!   95-11-20  ABELES     - PARALLEL OPTIMIZATION
!   96-03-29  BLACK      - REMOVED SCRCH COMMON
!   96-07-18  ZHAO       - NEW WMIN CALCULATION
!   96-09-25  BALDWIN    - NEW SR CALCULATION
!   98-11-02  BLACK      - MODIFICATION FOR DISTRIBUTED MEMORY
!-------------------------------------------------------------------------------
!     
! USAGE: CALL GSMDRIVE FROM gbphys
!
!   INPUT ARGUMENT LIST:
!       LM,DT,SL,DEL,PS,TIN,QIN,CCIN,
!       F_ice, F_rain,  F_RimeF, APREC, SR, 
!       ilat,   RHC, Xncwp,me
!  
!   OUTPUT ARGUMENT LIST: 
!     TIN, QIN, CCIN, F_ice, F_rain,  F_RimeF, APREC
!     
!   OUTPUT FILES:
!     NONE
!     
! Subprograms & Functions called:
!   GSMCONST  - initialize rain & ice lookup tables, read from external file;
!               initialize constants
!   GSMCOLUMN - cloud microphysics calculations over vertical columns
!
! UNIQUE: NONE
!  
! LIBRARY: NONE
!  
!?--- COMMON BLOCKS (input for microphysics):
!?       CTLBLK, LOOPS, MASKS, PHYS, VRBLS, CLDWTR, PVRBLS, ACMCLH, PPTASM, C_FRACN
!
!--- COMMON BLOCKS ("triggers" for microphysics & statistics):
!       CMICRO_START, CMICRO_STATS
!   
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
!------------------------------------------------------------------------
      implicit none
!
      integer, INTENT(IN   ) :: nCols
      integer, INTENT(IN   ) :: kMax
      real (kind=r8)  ,INTENT(IN   ) :: DT
      real (kind=r8)  ,INTENT(INOUT) :: APREC  (nCols)
      real (kind=r8)  ,INTENT(INOUT) :: ASNOW  (nCols)
      
      REAL(KIND=r8)   ,INTENT(INOUT) :: PREC  (1:nCols)
      REAL(KIND=r8)   ,INTENT(INOUT) :: ACPREC(1:nCols)
      REAL(KIND=r8)   ,INTENT(INOUT) :: SR    (1:nCols)
  
      real (kind=r8)  ,INTENT(INOUT) :: T_PHY  (nCols,kMax)
      real (kind=r8)  ,INTENT(IN   ) :: RHO    (nCols,kMax)
      real (kind=r8)  ,INTENT(IN   ) :: pgr    (nCols)
      real (kind=r8)  ,INTENT(INOUT) :: CWM_PHY   (nCols,kMax)
      real (kind=r8)  ,INTENT(INOUT) :: F_ICE_PHY (nCols,kMax)
      real (kind=r8)  ,INTENT(INOUT) :: F_RAIN_PHY (nCols,kMax)
      real (kind=r8)  ,INTENT(INOUT) :: F_RimeF_PHY(nCols,kMax)
      REAL(KIND=r8)   ,INTENT(INOUT) :: TLATGS_PHY  ( 1:nCols, 1:kMax )
      real (kind=r8)  ,INTENT(INOUT) :: Q_PHY    (nCols,kMax)
      REAL(KIND=r8)   ,INTENT(INOUT) :: qc         (1:nCols,1:kMax)
      REAL(KIND=r8)   ,INTENT(INOUT) :: qi         (1:nCols,1:kMax)
      REAL(KIND=r8)   ,INTENT(INOUT) :: qr         (1:nCols,1:kMax)
      REAL(KIND=r8)   ,INTENT(INOUT) :: TRAIN_PHY   ( 1:nCols, 1:kMax )
      real (kind=r8)  ,INTENT(IN   ) :: colrad (nCols)
     

    
!      INTEGER :: ntcw =1!     ntcw     - integer, cloud condensate location in the tracer  1    !
!                         array                                    1    !
!     nmtvr    - integer, number of topographic variables such as  1    !

!      real (kind=r8) :: gq0    (nCols,kMax,ntrac) !     gq0      - real, updated tracers x,kMax,ntrac!
        ! logical lprnt
!
!----------------------------------------------------------------------
!-----  Key parameters passed to column microphysics (COLUMN_MICRO) ------
!------------------------------------------------------------------------- 
!
!--- Flag from INIT.F at start of model run, used in initiating statistics
!
!     COMMON /CMICRO_START/ MICRO_START
!     LOGICAL :: MICRO_START
!
!--- This variable is for debugging purposes (if .true.)
!
!     LOGICAL, PARAMETER :: PRINT_diag=.TRUE.
!      LOGICAL PRINT_diag
!
!--- The following variables are for microphysical statistics (non-essential)
!
!     INTEGER, PARAMETER :: ITLO=-60, ITHI=40, ITHILO=ITHI-ITLO+1,
!    & ITHILO_N=ITHILO*4, ITHILO_QM=ITHILO*5, ITHILO_QT=ITHILO*22
!     COMMON /CMICRO_STATS/ NSTATS(ITLO:ITHI,4), QMAX(ITLO:ITHI,5),
!    & QTOT(ITLO:ITHI,22)
!     INTEGER :: NSTATS, NSTATS_0(ITLO:ITHI,4)
!     REAL :: QMAX, QTOT, QMAX_0(ITLO:ITHI,5),QTOT_0(ITLO:ITHI,22)
!     REAL, SAVE :: Thour_print, 
!    &  PRECmax(2),PRECtot(2),PRECmax_0(2),PRECtot_0(2)
!     REAL, PARAMETER :: DThour_print=3.     ! Print statistics every 3 h
!      REAL, PARAMETER :: DThour_print=0.     ! Print statistics every tnColse step
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~ BEGIN section on hydrometeor fractions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~ Saved values use REAL (REAL*4) arrays rather than INTEGER*2 
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!
!       real (kind=r8) ::   Fice
!       real (kind=r8) ::   Frain
       real (kind=r8) ::   DUM
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!~~~~~~ END section on hydrometeor fractions
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
!-----------------------------------------------------------------------
!-------------- Local arrays & parameters in GSMDRIVE -----------------
!-----------------------------------------------------------------------
!
!---- Comments on 14 March 2002
!    * EPSQ=1.E-12 is the universal lower limit for specific humidity and 
!      total condensate, and is consistent throughout the Eta code.
!
!     REAL, PARAMETER :: EPSQ=1.E-1_r82,  RHOL=1000._r8, T0C=273.15_r8, 
!      REAL(KIND=r8), PARAMETER :: EPSQ=1.0E-20_r8
      REAL(KIND=r8), PARAMETER :: RHOL=1000.0_r8
      REAL(KIND=r8), PARAMETER :: T0C=273.15_r8
      REAL(KIND=r8), PARAMETER :: T_ICE=-40.0_r8
      REAL(KIND=r8), PARAMETER :: T_ICEK=T0C+T_ICE
      REAL(KIND=r8), PARAMETER :: RRHOL=1.0_r8/RHOL
      REAL(KIND=r8), PARAMETER :: EPSQ1=1.001_r8*EPSQ
      
!    & T_ICE=-10._r8, T_ICEK=T0C+T_ICE, RRHOL=1._r8/RHOL, EPSQ1=1.001_r8*EPSQ
!
      real (kind=r8) ::  ARAIN
      real (kind=r8) ::  ASNOW2
      REAL(KIND=r8)    :: TRAIN    (1:nCols,1:kMax)
      REAL(KIND=r8)    :: CWM      (1:nCols,1:kMax)
      REAL(KIND=r8)    :: TLATGS   (1:nCols,1:kMax)
      real (kind=r8) ::  F_ice  (nCols,kMax)
      REAL(KIND=r8)    :: T        (1:nCols,1:kMax)

      real (kind=r8) ::  F_rain (nCols,kMax)
      real (kind=r8)  :: F_RimeF(nCols,kMax)
      real (kind=r8) ::  P_col (kMax)
      real (kind=r8) ::  QI_col(kMax)
      real (kind=r8) ::  QR_col(kMax)
      real (kind=r8) ::  QV_col(kMax)
      real (kind=r8) ::  QW_col(kMax)
      real (kind=r8) ::  RimeF_col(kMax)
      real (kind=r8) ::  T_col(kMax)
      real (kind=r8) ::  THICK_col(kMax)
      real (kind=r8) ::  WC_col(kMax)
      real (kind=r8) ::  NCW(kMax)
!
!
!      real (kind=r8) :: tc, wc, qi, qr, qw
      real (kind=r8) :: tc, wc

      integer        :: L, LL, i,k,KFLIP
      real (kind=r8) :: rhbbot
      real (kind=r8) :: rhpbl 
      real (kind=r8) :: rhbtop
      real (kind=r8) :: work1(nCols),  work2(nCols),work3(nCols),coslat(nCols)
      real (kind=r8) :: flgmin_l(nCols)
      real (kind=r8) :: RHC(nCols,kMax)  
      real (kind=r8) :: XNCW(nCols)
      real (kind=r8) :: RHC_col(kMax)
      real (kind=r8), parameter :: pt01=0.01_r8

      real (kind=r8)  :: prslk(nCols,kMax)   !     prslk    - real, Exner function at layer  		nCols,kMax !
      real (kind=r8)  :: prsik(nCols,kMax+1) !     prsik    - real, Exner function at layer interface	nCols,kMax+1
      real (kind=r8)  :: prsi (nCols,kMax+1)
      real (kind=r8)  :: prsl(nCols,kMax)!   PRSL P_col      - vertical column of model pressure (Pa)
      real (kind=r8)  :: pgrk(nCols)
      REAL(KIND=r8)   :: DPCOL    (1:kMax)    !GFDL
      REAL(KIND=r8)   :: DEL      (nCols,1:kMax)    !GFDL

!
!------------------------------------------------------------------------
!
!#######################################################################
!########################## Begin Execution ############################
!#######################################################################
!
!------------------------------------------------------------------------
!---------------------- Microphysical constants -------------------------
!------------------------------------------------------------------------
!
!
!  move water from vapor to liquid should the liquid amount be negative

      do i=1,nCols
         prsi(i,kMax+1)  = si(kMax+1)*pgr(i)           ! prsi are now pressures(Pa)
         pgrk(i)       = (pgr(i)*pt01) ** con_rocp
         prsik(i,kMax+1) = sik(kMax+1) * pgrk(i)
      enddo

      do k=1,kmax
        do i=1,nCols
          prsi(i,k)  = si(k)*pgr(i)               ! prsi are now pressures(Pa)
          prsl(i,k)  = sl(k)*pgr(i)               !   P_col      - vertical column of model pressure (Pa)
          prsik(i,k) = sik(k) * pgrk(i)
          prslk(i,k) = slk(k) * pgrk(i)
        enddo
      enddo
      do k=1,kmax
         do i=1,nCols
           del(i,k) = PRSI(i,k) - PRSI(i,k+1)
         END DO
      END DO 
!      do L = 1, kMax
!        do i=1,nCols
!          if (CWM_PHY(i,L) .lt. 0.0_r8) then
!            Q_PHY(i,L)  = Q_PHY(i,L) + CWM_PHY(i,L)
!            if (T_PHY(i,l) .gt. t_icek) then
!              T_PHY(i,L)  = T_PHY(i,L) - CWM_PHY(i,L) * (con_hvap/con_cp)
!            else
!              T_PHY(i,L)  = T_PHY(i,L) - CWM_PHY(i,L) * (con_hsub/con_cp)
!            endif
!            CWM_PHY(i,L) = 0.0_r8
!          endif
!        enddo
!      enddo

      do i = 1, nCols
       ! colrad.....colatitude  colrad=0 - 3.14 (0-180)from np to sp in radians
       !IF((((colrad(i)*180.0_r8)/3.1415926e0_r8)-90.0_r8)  > 0.0_r8 ) THEN
!        coslat(i)   = (((colrad(i)*180.0_r8)/3.1415926e0_r8)-90.0_r8)
        coslat(i)   = cos(((colrad(i)))-(3.1415926e0_r8/2.0_r8))
         !log(1.0_r8/(REAL(iMax/(20))*REAL(jMax/(20))))

        work1(i)    = (log(coslat(i) / (nlons*latg)) - dxmin) * dxinv
        work1(i)    = max(0.0_r8, min(1.0_r8,work1(i)))
        work2(i)    = 1.0_r8 - work1(i)

        work3(i)    = prsik(i,1) / prslk(i,1)
!         PRINT*,work1(i),work2(i),coslat(i) ,colrad(i)

      enddo

      rhbbot = crtrh(1)
      rhpbl  = crtrh(2)
      rhbtop = crtrh(3)

      do k = 1, kMax
         do i = 1, nCols
            rhc(i,k) = rhbbot - (rhbbot-rhbtop) * (1.0_r8-prslk(i,k))
            rhc(i,k) = rhc_max * work1(i) + rhc(i,k) * work2(i)
            rhc(i,k) = max(0.0_r8, min(1.0_r8,rhc(i,k)))
         enddo
      enddo

      do i = 1, nCols
         flgmin_l(i) = flgmin(1)*work1(i) + flgmin(2)*work2(i)
      enddo
      
      do i = 1, nCols
         xncw(i) = ncwp(2) * work1(i) + ncwp(1) * work2(i)
      enddo    


      DO   I=1,nCols  
       DO L=1,kMax
          KFLIP=kMax+1-L
          CWM    (I,KFLIP)=CWM_PHY    (I,L)
          T      (I,KFLIP)=T_PHY      (I,L)
          TLATGS (I,KFLIP)=TLATGS_PHY (I,L)
          TRAIN  (I,KFLIP)=TRAIN_PHY  (I,L)
          F_ice  (I,KFLIP)=F_ice_PHY  (I,L)
          F_rain (I,KFLIP)=F_rain_PHY (I,L)
          F_RimeF(I,KFLIP)=F_RimeF_PHY(I,L)
       ENDDO
      END DO

!
!------------------------------------------------------------------------
!--------------- Initialize constants for statistics --------------------
!------------------------------------------------------------------------
!
!       Thour_print=-DTPH/3600.0_r8+FLOAT(NTSD-1)*DT/3600.0_r8
!       IF (PRINT_diag) THEN
!
!-------- Total and maximum quantities
!
!         DO I=ITLO,ITHI
!--- Microphysical statistics dealing w/ grid-point counts
!           DO J=1,4
!             NSTATS(I,J)=0
!           ENDDO
!--- Microphysical statistics dealing w/ maxima of hydrometeor mass
!           DO J=1,5
!             QMAX(I,J)=0.0_r8
!           ENDDO
!--- Microphysical statistics dealing w/ total hydrometeor mass
!           DO J=1,22
!             QTOT(I,J)=0.0_r8
!           ENDDO
!         ENDDO
!         DO I=1,2
!           PRECmax(I)=0.0_r8    ! Maximum precip rates (rain, snow) at surface (mm/h)
!           PRECtot(I)=0.0_r8    ! Total precipitation (rain, snow) accumulation at surface
!         ENDDO
!       ENDIF
!     ENDIF
!
      do i=1,nCols              ! Begining of the I loop!
!
!       if (lprnt .and. i .eq. ipr) then
!          PRINT_diag = .true.
!       else
!           PRINT_diag = .false.
!       endif
!       IF (PRINT_diag) THEN
!         print *,' printing for i=',i,' me=',me
!         print *,' CWM_PHY=',CWM_PHY(ipr,:)
!         print *,' Q_PHY=',Q_PHY(ipr,:)
!         print *,' F_rain=',F_rain(ipr,:)
!       endif
!
!--- Initialize column data (1D arrays)
!
      DO L=1,kMax
        LL = kMax + 1 - L
!        DPCOL    (L)=RHO(I,LL)*con_g*del(I,LL)
!              THICK  = THICK_col(L)   ! Layer thickness = RHO*DZ = -DP/G

!        THICK_col(L)  = DEL(I,LL) * (1.0_r8/con_g) !--- Layer thickness = RHO*DZ

        THICK_col(L) =del(I,LL)*RGRAV

        P_col    (L)  = PRSL(I,LL)
        T_col    (L)  = T_PHY(I,LL)
        QV_col   (L)  = max(EPSQ, Q_PHY(I,LL))
        RHC_col  (L)  = RHC(I,LL)
        WC_col   (L)  = CWM_PHY(I,LL)
        NCW      (L)  = XNCW(I)
!         PRINT*,NCW(L),XNCW(I)

      ENDDO
!     if (print_diag) print *,' wc_col=',wc_col
      DO L=1,kMax
        LL = kMax + 1 - L
        TC           = T_col(L)-T0C
        IF (WC_col(L) .LE. EPSQ1) THEN
          WC_col(L)  = 0.0_r8
          IF (TC .LT. T_ICE) THEN
            F_ice(I,LL) = 1.0_r8
          ELSE
            F_ice(I,LL) = 0.0_r8
          ENDIF
          F_rain (I,LL)  = 0.0_r8
          F_RimeF(I,LL)  = 1.0_r8
        ENDIF
!
!--- Determine composition of condensate in terms of
!      cloud water, ice, & rain
!
        WC    = WC_col(L)
        !QI    = 0.0_r8
        !QR    = 0.0_r8
        !QW    = 0.0_r8
!        Fice  = F_ice (I,LL)
!        Frain = F_rain(I,LL)
!
!--- REAL*4 array storage
!
        !IF (Fice .GE. 1.0_r8) THEN
        !  !QI = WC
        !ELSE IF (Fice .LE. 0.0_r8) THEN
        !  !QW = WC
        !ELSE
        !  QI = Fice*WC
        !  QW = WC-QI
        !ENDIF
        !IF (QW.GT.0.0_r8 .AND. Frain.GT.0.0_r8) THEN
        !  IF (Frain .GE. 1.0_r8) THEN
        !    QR = QW
        !    QW = 0.0_r8
        !  ELSE
        !    QR = Frain*QW
        !    QW = QW-QR
        !  ENDIF
        !ENDIF
        RimeF_col(L) = F_RimeF(I,LL)              ! (real)
!        QI_col(L) = QI
!        QR_col(L) = QR
!        QW_col(L) = QW
        QI_col(L) = qi         (i,LL)
        QR_col(L) = qr         (i,LL)
        QW_col(L) = qc         (i,LL)

      ENDDO
!
!#######################################################################
!
!--- Perform the microphysical calculations in this column
!
!      PRINT*,NCW

      CALL GSMCOLUMN ( &
      ARAIN     , & !real(kind=r8)   , INTENT(OUT  ) :: ARAING
      ASNOW2    , & !real(kind=r8)   , INTENT(OUT  ) :: ASNOWG
      DT        , & !real(kind=r8)   , INTENT(IN   ) :: dtpg
      kMax      , & !INTEGER     , INTENT(IN   ) :: LSFC
      P_col     , & !real(kind=r8)   , INTENT(IN   ) :: P_col(kMax)
      QI_col    , & !real(kind=r8)   , INTENT(INOUT) :: QI_col(kMax)
      QR_col    , & !real(kind=r8)   , INTENT(INOUT) :: QR_col(kMax)
      QV_col    , & !real(kind=r8)   , INTENT(INOUT) :: QV_col(kMax)
      QW_col    , & !real(kind=r8)   , INTENT(INOUT) :: QW_col(kMax)
      RimeF_col , & !real(kind=r8)   , INTENT(INOUT) :: RimeF_col(kMax)
      T_col     , & !real(kind=r8)   , INTENT(INOUT) :: T_col(kMax)
      THICK_col , & !real(kind=r8)   , INTENT(IN   ) :: THICK_col(kMax)
      WC_col    , & !real(kind=r8)   , INTENT(INOUT) :: WC_col(kMax)
      kMax        , & !integer	     , INTENT(IN   ) :: kMax
      RHC_col   , & !real(kind=r8)   , INTENT(IN   ) :: RHC_col(kMax)
      NCW       , & !real(kind=r8)   , INTENT(IN   ) :: XNCW(kMax)
      flgmin_l(i) ) !real(kind=r8)   , INTENT(IN   ) :: flgmin
!
!#######################################################################
!
!
!--- Update storage arrays
!
      DO L=1,kMax
        LL = kMax + 1 - L
        TRAIN (I,LL)  = (T_col(L)-T(I,L))/DT
        TLATGS(I,LL)  = T_col (L)-T(I,L)
        T_PHY (I,LL)  = T_col(L)
        IF (Q_PHY(I,LL) .LT. EPSQ) THEN
          Q_PHY(I,LL)  = Q_PHY(I,LL) + QV_col(L)
        else
          Q_PHY(I,LL)  = QV_col(L)
        endif
        IF (CWM_PHY(I,LL) .LT. EPSQ) THEN
          CWM_PHY(I,LL) = CWM_PHY(I,LL) + WC_col(L)
        else
          CWM_PHY(I,LL) = WC_col(L)
        endif
!
!--- REAL*4 array storage
!
        F_RimeF(I,LL)=MAX(1.0_r8, RimeF_col(L))
        IF (QI_col(L) .LE. EPSQ) THEN
          F_ice(I,LL)=0.0_r8
          IF (T_col(L) .LT. T_ICEK) F_ice(I,LL)=1.0_r8
        ELSE
          F_ice(I,LL)=MAX( 0.0_r8, MIN(1.0_r8, QI_col(L)/ max(EPSQ, WC_col(L))) )
        ENDIF
        IF (QR_col(L) .LE. EPSQ) THEN
          DUM=0
        ELSE
          DUM=QR_col(L)/MAX(QR_col(L)+QW_col(L),EPSQ)
        ENDIF
        F_rain(I,LL)=DUM
!
!
        qi   (i,LL)=  QI_col(L)!PK
        qr   (i,LL)=  QR_col(L)!PK
        qc   (i,LL)=  QW_col(L)!PK
      ENDDO
!
!
!--- Update accumulated precipitation statistics
!
!--- Surface precipitation statistics; SR is fraction of surface
!    precipitation (if >0) associated with snow
!
      APREC(I) = (ARAIN)*RRHOL    ! Accumulated surface precip (depth in m)
      ASNOW(I) = (ASNOW2)*RRHOL       ! Accumulated surface snow (depth in m)  !<--- Ying
      PREC  (I)= PREC (I)  + (ARAIN+ASNOW2)*RRHOL
      ACPREC(I)=ACPREC(I)  + (ARAIN+ASNOW2)*RRHOL
       IF((ARAIN+ASNOW2)*RRHOL .LT. 1.E-8_r8) THEN
          SR(I)  = 0.0_r8
      ELSE
          SR(I)=RRHOL*ASNOW2/(ARAIN+ASNOW2)*RRHOL
      ENDIF
!
!--- Debug statistics
!
!       IF (PRINT_diag) THEN
!         PRECtot(1)=PRECtot(1)+ARAIN
!         PRECtot(2)=PRECtot(2)+ASNOW2
!         PRECmax(1)=MAX(PRECmax(1), ARAIN)
!         PRECmax(2)=MAX(PRECmax(2), ASNOW2)
!       ENDIF
!#######################################################################
!#######################################################################
!
!-----------------------------------------------------------------------
!--------------------- END of main microphysics loop -------------------
!-----------------------------------------------------------------------
!
      ENDDO              ! End of the I loop
!
    DO  I=1,nCols  
       DO L=1,kMax
          KFLIP=kMax+1-L
          TLATGS_PHY (I,L)=TLATGS (I,KFLIP)
          TRAIN_PHY  (I,L)=TRAIN  (I,KFLIP)
          F_ice_PHY  (I,L)=F_ice  (I,KFLIP)
          F_rain_PHY (I,L)=F_rain (I,KFLIP)
          F_RimeF_PHY(I,L)=F_RimeF(I,KFLIP)
       ENDDO
    END DO

      RETURN
!-----------------------------------------------------------------------
!200   format(a2,i5,f6.2,4(1x,a10,g11.4))
!210   format(a2,i5,f6.2,4(1x,a10,i7))
!-----------------------------------------------------------------------
      END SUBROUTINE GSMDRIVE
!
!###############################################################################
! ***** VERSION OF MICROPHYSICS DESIGNED FOR HIGHER RESOLUTION MESO ETA MODEL
!       (1) Represents sedimentation by preserving a portion of the precipitation
!           through top-down integration from cloud-top.  Modified procedure to
!           Zhao and Carr (1997).
!       (2) Microphysical equations are modified to be less sensitive to time
!           steps by use of Clausius-Clapeyron equation to account for changes in
!           saturation mixing ratios in response to latent heating/cooling.  
!       (3) Prevent spurious temperature oscillations across 0C due to 
!           microphysics.
!       (4) Uses lookup tables for: calculating two different ventilation 
!           coefficients in condensation and deposition processes; accretion of
!           cloud water by precipitation; precipitation mass; precipitation rate
!           (and mass-weighted precipitation fall speeds).
!       (5) Assumes temperature-dependent variation in mean diameter of large ice
!           (Houze et al., 1979; Ryan et al., 1996).
!        -> 8/22/01: This relationship has been extended to colder temperatures
!           to parameterize smaller large-ice particles down to mean sizes of MDImin,
!           which is 50 microns reached at -55.9C.
!       (6) Attempts to differentiate growth of large and small ice, mainly for
!           improved transition from thin cirrus to thick, precipitating ice
!           anvils.
!        -> 8/22/01: This feature has been diminished by effectively adjusting to
!           ice saturation during depositional growth at temperatures colder than
!           -10C.  Ice sublimation is calculated more explicitly.  The logic is
!           that sources of are either poorly understood (e.g., nucleation for NWP) 
!           or are not represented in the Eta model (e.g., detrainment of ice from 
!           convection).  Otherwise the model is too wet compared to the radiosonde
!           observations based on 1 Feb - 18 March 2001 retrospective runs.  
!       (7) Top-down integration also attempts to treat mixed-phase processes,
!           allowing a mixture of ice and water.  Based on numerous observational
!           studies, ice growth is based on nucleation at cloud top &
!           subsequent growth by vapor deposition and riming as the ice particles 
!           fall through the cloud.  Effective nucleation rates are a function
!           of ice supersaturation following Meyers et al. (JAM, 1992).  
!        -> 8/22/01: The simulated relative humidities were far too moist compared 
!           to the rawinsonde observations.  This feature has been substantially 
!           diminished, limited to a much narrower temperature range of 0 to -10C.  
!       (8) Depositional growth of newly nucleated ice is calculated for large time
!           steps using Fig. 8 of Miller and Young (JAS, 1979), at 1 deg intervals
!           using their ice crystal masses calculated after 600 s of growth in water
!           saturated conditions.  The growth rates are normalized by time step
!           assuming 3D growth with time**1.5 following eq. (6.3) in Young (1993).
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!       (9) Ice precipitation rates can increase due to increase in response to
!           cloud water riming due to (a) increased density & mass of the rimed
!           ice, and (b) increased fall speeds of rimed ice.
!        -> 8/22/01: This feature has been effectively limited to 0 to -10C.  
!###############################################################################
!###############################################################################
!
      SUBROUTINE GSMCOLUMN (&
       ARAING    , &!real(kind=r8)   , INTENT(OUT  ) :: ARAING
       ASNOWG    , &!real(kind=r8)   , INTENT(OUT  ) :: ASNOWG
       DTPG      , &!real(kind=r8)   , INTENT(IN   ) :: dtpg
       LSFC      , &!INTEGER         , INTENT(IN   ) :: LSFC
       P_col     , &!real(kind=r8)   , INTENT(IN   ) :: P_col(LM)
       QI_col    , &!real(kind=r8)   , INTENT(INOUT) :: QI_col(LM)
       QR_col    , &!real(kind=r8)   , INTENT(INOUT) :: QR_col(LM)
       QV_col    , &!real(kind=r8)   , INTENT(INOUT) :: QV_col(LM)
       QW_col    , &!real(kind=r8)   , INTENT(INOUT) :: QW_col(LM)
       RimeF_col , &!real(kind=r8)   , INTENT(INOUT) :: RimeF_col(LM)
       T_col     , &!real(kind=r8)   , INTENT(INOUT) :: T_col(LM)
       THICK_col , &!real(kind=r8)   , INTENT(IN   ) :: THICK_col(LM)
       WC_col    , &!real(kind=r8)   , INTENT(INOUT) :: WC_col(LM)
       LM        , &!integer         , INTENT(IN   ) :: lm
       RHC_col   , &!real(kind=r8)   , INTENT(IN   ) :: RHC_col(LM)
       XNCW      , &!real(kind=r8)   , INTENT(IN   ) :: XNCW(LM)
       flgmin_l    )!real(kind=r8)   , INTENT(IN   ) :: flgmin
!
      implicit none
!
!###############################################################################
!###############################################################################
!
!-------------------------------------------------------------------------------
!----- NOTE:  In this version of the Code threading should be done outside!  
!-------------------------------------------------------------------------------
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:  Grid-scale microphysical processes - condensation & precipitation
!   PRGRMMR: Ferrier         ORG: W/NP22     DATE: 08-2001
!   Updated: Moorthi for GFS application
!-------------------------------------------------------------------------------
! ABSTRACT:
!   * Merges original GSCOND & PRECPD subroutines.   
!   * Code has been substantially streamlined and restructured.
!   * Exchange between water vapor & small cloud condensate is calculated using
!     the original Asai (1965, J. Japan) algorithm.  See also references to
!     Yau and Austin (1979, JAS), Rutledge and Hobbs (1983, JAS), and Tao et al.
!     (1989, MWR).  This algorithm replaces the Sundqvist et al. (1989, MWR)
!     parameterization.  
!-------------------------------------------------------------------------------
!     
! USAGE: 
!   * CALL GSMCOLUMN FROM SUBROUTINE GSMDRIVE
!   * SUBROUTINE GSMDRIVE CALLED FROM MAIN PROGRAM EBU
!
! INPUT ARGUMENT LIST:
!   DTPH       - physics time step (s)
!   I_index    - I index
!   J_index    - J index
!   LSFC       - Eta level of level above surface, ground
!   P_col      - vertical column of model pressure (Pa)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!   THICK_col  - vertical column of model mass thickness (density*height increment)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   
!
! OUTPUT ARGUMENT LIST: 
!   ARAIN      - accumulated rainfall at the surface (kg)
!   ASNOW      - accumulated snowfall at the surface (kg)
!   QV_col     - vertical column of model water vapor specific humidity (kg/kg)
!   WC_col     - vertical column of model mixing ratio of total condensate (kg/kg)
!   QW_col     - vertical column of model cloud water mixing ratio (kg/kg)
!   QI_col     - vertical column of model ice mixing ratio (kg/kg)
!   QR_col     - vertical column of model rain ratio (kg/kg)
!   RimeF_col  - vertical column of rime factor for ice in model (ratio, defined below)
!   T_col      - vertical column of model temperature (deg K)
!     
! OUTPUT FILES:
!     NONE
!     
! Subprograms & Functions called:
!   * real(kind=r8) Function CONDENSE  - cloud water condensation
!   * real(kind=r8) Function DEPOSIT   - ice deposition (not sublimation)
!
! UNIQUE: NONE
!  
! LIBRARY: NONE
!  
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!   MACHINE : IBM SP
!
!------------------------------------------------------------------------- 
!--------------- Arrays & constants in argument list --------------------- 
!------------------------------------------------------------------------- 
!
      integer, INTENT(IN   ) :: lm
      INTEGER, INTENT(IN   ) :: LSFC
      real(kind=r8)   , INTENT(OUT  ) :: ARAING
      real(kind=r8)   , INTENT(OUT  ) :: ASNOWG
      real(kind=r8)   , INTENT(IN   ) :: P_col (LM)
      real(kind=r8)   , INTENT(INOUT) :: QI_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: QR_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: QV_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: QW_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: RimeF_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: T_col(LM)
      real(kind=r8)   , INTENT(IN   ) :: THICK_col(LM)
      real(kind=r8)   , INTENT(INOUT) :: WC_col(LM)
      real(kind=r8)   , INTENT(IN   ) :: RHC_col(LM)
      real(kind=r8)   , INTENT(IN   ) :: XNCW(LM)
      real(kind=r8)   , INTENT(IN   ) :: dtpg
      real(kind=r8)   , INTENT(IN   ) :: flgmin_l
!
!
!
!------------------------------------------------------------------------- 
!
!--- Mean ice-particle diameters varying from 50 microns to 1000 microns
!      (1 mm), assuming an exponential size distribution.  
!
!---- Meaning of the following arrays: 
!        - mdiam - mean diameter (m)
!        - VENTI1 - integrated quantity associated w/ ventilation effects 
!                   (capacitance only) for calculating vapor deposition onto ice
!        - VENTI2 - integrated quantity associated w/ ventilation effects 
!                   (with fall speed) for calculating vapor deposition onto ice
!        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
!        - MASSI  - integrated quantity associated w/ ice mass 
!        - VSNOWI - mass-weighted fall speed of snow (large ice), used to calculate 
!                   precipitation rates
!
      real(kind=r8),    PARAMETER :: DMImin=.05e-3_r8
      real(kind=r8),    PARAMETER :: DMImax=1.e-3_r8
      real(kind=r8),    PARAMETER :: DelDMI=1.e-6_r8
      real(kind=r8),    PARAMETER :: XMImin=1.e6_r8*DMImin
      real(kind=r8),    PARAMETER :: XMImax=1.e6_r8*DMImax
      INTEGER, PARAMETER :: MDImin=INT(XMImin)
      INTEGER, PARAMETER :: MDImax=INT(XMImax)
!
!------------------------------------------------------------------------- 
!------- Key parameters, local variables, & important comments ---------
!-----------------------------------------------------------------------
!
!--- KEY Parameters:
!
!---- Comments on 14 March 2002
!    * Set EPSQ to the universal value of 1.e-12 throughout the code
!      condensate.  The value of EPSQ will need to be changed in the other 
!      subroutines in order to make it consistent throughout the Eta code.  
!    * Set CLIMIT=10.*EPSQ as the lower limit for the total mass of 
!      condensate in the current layer and the input flux of condensate
!      from above (TOT_ICE, TOT_ICEnew, TOT_RAIN, and TOT_RAINnew).
!
!-- NLImax - maximum number concentration of large ice crystals (20,000 /m**3, 20 per liter)
!-- NLImin - minimum number concentration of large ice crystals (100 /m**3, 0.1 per liter)
!
      real(kind=r8), PARAMETER ::   RHOL=1000.0_r8
      real(kind=r8), PARAMETER ::   XLS=con_hvap+con_hfus

!    &, T_ICE=-10.          !- Ver1
!    &, T_ICE_init=-5.      !- Ver1
!!!  &, T_ICE=-20.          !- Ver2
!    &, T_ICE=-40.          !- Ver2
!    &, T_ICE_init=-15.,    !- Ver2
!
!    & CLIMIT=10.*EPSQ, EPS1=RV/RD-1., RCP=1./con_cp,
      real(kind=r8), PARAMETER :: EPS1=con_rv/con_rd-1.0_r8
      real(kind=r8), PARAMETER :: CLIMIT=10.0_r8*EPSQ
      real(kind=r8), PARAMETER :: RCP=1.0_r8/con_cp
      real(kind=r8), PARAMETER :: RCPRV=RCP/con_rv
      real(kind=r8), PARAMETER :: RRHOL=1.0_r8/RHOL
      real(kind=r8), PARAMETER :: XLS1=XLS*RCP
      real(kind=r8), PARAMETER :: XLS2=XLS*XLS*RCPRV
      real(kind=r8), PARAMETER :: XLS3=XLS*XLS/con_rv
      real(kind=r8), PARAMETER :: C1=1.0_r8/3.0_r8
      real(kind=r8), PARAMETER :: C2=1.0_r8/6.0_r8
      real(kind=r8), PARAMETER :: C3=3.31_r8/6.0_r8
      real(kind=r8), PARAMETER :: DMR1=0.1E-3_r8
      real(kind=r8), PARAMETER :: DMR2=0.2E-3_r8
      real(kind=r8), PARAMETER :: DMR3=0.32E-3_r8
      real(kind=r8), PARAMETER :: N0r0=8.E6_r8
      real(kind=r8), PARAMETER :: N0rmin=1.e4_r8

      real(kind=r8), PARAMETER :: N0s0=4.E6_r8
      real(kind=r8), PARAMETER :: RHO0=1.194_r8
      real(kind=r8), PARAMETER :: XMR1=1.e6_r8*DMR1
      real(kind=r8), PARAMETER :: XMR2=1.e6_r8*DMR2
      real(kind=r8), PARAMETER :: XMR3=1.e6_r8*DMR3
      real(kind=r8)   , PARAMETER :: Xratio=.025_r8
      INTEGER, PARAMETER :: MDR1=INT(XMR1)
      INTEGER, PARAMETER :: MDR2=INT(XMR2)
      INTEGER, PARAMETER :: MDR3=INT(XMR3)
!
!--- If BLEND=1:
!      precipitation (large) ice amounts are estimated at each level as a 
!      blend of ice falling from the grid point above and the precip ice
!      present at the start of the time step (see TOT_ICE below).
!--- If BLEND=0:
!      precipitation (large) ice amounts are estimated to be the precip
!      ice present at the start of the time step.
!
!--- Extended to include sedimentation of rain on 2/5/01 
!
      real(kind=r8), PARAMETER :: BLEND=1.0_r8
!
!--- This variable is for debugging purposes (if .true.)
!
!      LOGICAL  PRINT_diag
!
!--- Local variables
!
      real(kind=r8)    :: EMAIRI, N0r,         NLICE,       NSmICE, NLImax, pfac
      LOGICAL :: CLEAR,  ICE_logical, DBG_logical, RAIN_logical
 
      integer :: lbef, ipass, ixrf, ixs,  idr                         &
     &,       index_my, indexr, indexr1, indexs                         &
     &,       l, ntimes, item
!    &,       i, j, k, my_600, i1, i2, l, ntimes

      real(kind=r8) flimass,  xlimass, vsnow,   qi_min, dum,    piloss           &
     &,    tot_ice,  xsimass, vel_inc, vrimef, rimef1, dum1             &
     &,    dum2,     fws,     denomi                              &
     &,    xrf,      qw0,     dli,     xli,    fsmall                   &
     &,    prevp,    tk2,     dtph                                      &
     &,    pievp,    picnd,   piacr,   pracw                            &
     &,    praut,    pimlt,   qtice,   qlice                            &
     &,    gammar,   flarge,  wvqw,    dynvis                           &
     &,    tfactor,  gammas,  diffus, therm_cond               &
     &,    wvnew,    delv,    tnew,    tot_icenew, rimef                &
     &,    deli,     fwr,     crevp,   ventr,      delt                 &
     &,    delw,     fir,     delr               &
     &,    budget,   vrain2,  tot_rainnew                      &
     &,    qtnew,    qt,      wcnew,   abw                              &
     &,    aievp,    tcc,     denomf,  abi                              &
     &,    sfactor,  pidep_max,        didep,      ventis, ventil       &
     &,    dievp,    rqr,     rfactor, dwvr,       rr,     tot_rain     &
     &,    dwv0,     qsw0,    prloss,  qtrain,     vrain1               &
     &,    qsw,      ws,      esi,     esw, wv, wc, rhgrd, rho          &
     &,    rrho,     dtrho,   wsgrd,   qsi, qswgrd, qsigrd              &
     &,    tk,       tc,      pp,      bldtrh                           &
     &,    xlv,      xlv1,    xlf,     xlf1,  xlv2, denomw, denomwi     &
     &,    qwnew,    pcond,   pidep,   qrnew, qi,   qr,     qw          &
     &,    piacw,    piacwi,  piacwr,  qv,    dwvi                      &
     &,    arainnew, thick,   asnownew                                  &
     &,    qinew,    qi_min_0c, QSW_l, QSI_l, QSW0_l, SCHMIT_FAC
      real(kind=r8)    :: ARAIN
      real(kind=r8)    :: ASNOW

!
!
!#######################################################################
!########################## Begin Execution ############################
!#######################################################################
!
      INDEXR = 0
      DTPH   = DTPG / mic_step
      ARAING = 0.0_r8    ! Total Accumulated rainfall at surface (kg/m**2)
      ASNOWG = 0.0_r8    ! Total Accumulated snowfall at surface (kg/m**2)
!
      do ntimes =1,mic_step
!
        QI_min_0C = 10.E3_r8*MASSI(MDImin)   !- Ver5
        ARAIN = 0.0_r8   ! Accumulated rainfall at surface for this step (kg/m**2)
        ASNOW = 0.0_r8   ! Accumulated snowfall at surface for this step (kg/m**2)
!
!-----------------------------------------------------------------------
!
        DO L=1,LSFC      !      Loop from top (L=1) to surface (L=LSFC)

!---      Skip this level and go to the next lower level if no condensate 
!         and very low specific humidities
!
          IF (QV_col(L) > EPSQ .OR. WC_col(L) > EPSQ) THEN
!
!-----------------------------------------------------------------------
!------------ Proceed with cloud microphysics calculations -------------
!-----------------------------------------------------------------------
!
            TK = T_col(L)         ! Temperature (deg K)
            TC = TK-T0C           ! Temperature (deg C)
            PP = P_col(L)         ! Pressure (Pa)
            QV = QV_col(L)        ! Specific humidity of water vapor (kg/kg)
!           WV = QV/(1.-QV)       ! Water vapor mixing ratio (kg/kg)
            WV = QV               ! Water vapor specific humidity (kg/kg)
            WC = WC_col(L)        ! Grid-scale mixing ratio of total condensate
                                ! (water or ice; kg/kg)
!           WC = WC/(1.-WC)
            RHgrd = RHC_col(L)

!
!   Pressure dependen scaling factor for flgmin_l (tunable)
!
!!!         pfac = max(0.5, (min(1.0, pp*0.00002))**2)   ! commented on 02182011
!go         pfac = max(0.5, (sqrt(min(1.0, pp*0.00004))))
            pfac = 1.0_r8
!
            CLEAR = .TRUE.
!    
!--- Check grid-scale saturation when no condensate is present
!    
            ESW = min(PP, FPVSL(TK))     ! Saturation vapor pressure w/r/t water
!           QSW = con_eps*ESW/(PP-ESW)       ! Saturation mixing ratio w/r/t water
            QSW = con_eps*ESW/(PP+con_epsm1*ESW) ! Saturation specific humidity  w/r/t water
            WS  = QSW                    ! General saturation mixing ratio (water/ice)
            QSI = QSW
            IF (TC < 0.0_r8) THEN
              ESI = min(PP,FPVSI(TK))      ! Saturation vapor pressure w/r/t ice
!             QSI = con_eps*ESI/(PP-ESI)       ! Saturation mixing ratio w/r/t water
              QSI = con_eps*ESI/(PP+con_epsm1*ESI) ! Saturation specific humidity w/r/t water
              WS  = QSI                    ! General saturation mixing ratio (water/ice)
              if (pp <= esi) ws = wv / rhgrd
            ENDIF
!
            dum  = min(PP, ESW0)
            QSW0 = con_eps*dum/(PP+con_epsm1*dum)  ! Saturation specific Humidity at 0C
!
!--- Effective grid-scale Saturation mixing ratios
!
            QSWgrd = RHgrd*QSW
            QSIgrd = RHgrd*QSI
            WSgrd  = RHgrd*WS
            QSW_l  = QSWgrd
            QSI_l  = QSIgrd
            QSW0_l = QSW0*RHgrd
!
!--- Check if air is subsaturated and w/o condensate
!
            IF (WV > WSgrd .OR. WC > EPSQ) CLEAR = .FALSE.  ! Cloudy case
            IF (ARAIN > CLIMIT) THEN ! If any rain is falling into layer from above
              CLEAR = .FALSE.
            ELSE
              ARAIN = 0.0_r8
            ENDIF
!
!--- Check if any ice is falling into layer from above
!
!--- NOTE that "SNOW" in variable names is synonomous with 
!    large, precipitation ice particles
!
            IF (ASNOW > CLIMIT) THEN
              CLEAR = .FALSE.
            ELSE
              ASNOW = 0.0_r8
            ENDIF
!
!-----------------------------------------------------------------------
!-- Loop to the end if in clear, subsaturated air free of condensate ---
!-----------------------------------------------------------------------
!
            IF (.not. CLEAR) THEN
!
!-----------------------------------------------------------------------
!--------- Initialize RHO, THICK & microphysical processes -------------
!-----------------------------------------------------------------------
!
!
!--- Virtual temperature, TV=T*(1./con_eps-1)*Q, Q is specific humidity;
!    (see pp. 63-65 in Fleagle & Businger, 1963)
!
              RHO    = PP/(con_rd*TK*(1.0_r8+EPS1*QV)) ! Air density (kg/m**3)
              RRHO   = 1.0_r8/RHO                  ! Reciprocal of air density
              DTRHO  = DTPH*RHO                ! Time step * air density
              BLDTRH = BLEND*DTRHO             ! Blend parameter * time step * air density
              THICK  = THICK_col(L)   ! Layer thickness = RHO*DZ = -DP/G
!
              ARAINnew = 0.0_r8           ! Updated accumulated rainfall at surface
              ASNOWnew = 0.0_r8           ! Updated accumulated snowfall at surface
              QI       = QI_col(L)    ! Ice mixing ratio
              QInew    = 0.0_r8           ! Updated ice mixing ratio
              QR       = QR_col(L)    ! Rain mixing ratio
              QRnew    = 0.0_r8           ! Updated rain ratio
              QW       = QW_col(L)    ! Cloud water mixing ratio
              QWnew    = 0.0_r8           ! Updated cloud water ratio
!
              PCOND    = 0.0_r8       ! Condensation (>0) or evaporation (<0)
                                  ! of cloud water (kg/kg)
              PIDEP    = 0.0_r8       ! Deposition (>0) or sublimation (<0)
                                  ! of ice crystals (kg/kg)
              PIACW   = 0.0_r8        ! Cloud water collection (riming)
                                  ! by precipitation ice (kg/kg; >0)
              PIACWI  = 0.0_r8        ! Growth of precip ice by riming (kg/kg; >0)
              PIACWR  = 0.0_r8        ! Shedding of accreted cloud water
                                  ! to form rain (kg/kg; >0)
              PIACR   = 0.0_r8        ! Freezing of rain onto large ice
                                  ! at supercooled temps (kg/kg; >0)
              PICND   = 0.0_r8        ! Condensation (>0) onto wet, melting
                                  ! ice (kg/kg)
              PIEVP   = 0.0_r8        ! Evaporation (<0) from wet, melting
                                  ! ice (kg/kg)
              PIMLT   = 0.0_r8        ! Melting ice (kg/kg; >0)
              PRAUT   = 0.0_r8        ! Cloud water autoconversion to rain (kg/kg; >0)
              PRACW   = 0.0_r8        ! Cloud water collection (accretion) by rain (kg/kg; >0)
              PREVP   = 0.0_r8        ! Rain evaporation (kg/kg; <0)
!
!---------------------------------------------------------------------------
!--- Double check input hydrometeor mixing ratios
!
!             DUM  = WC - (QI+QW+QR)
!             DUM1 = ABS(DUM)
!             DUM2 = TOLER * MIN(WC, QI+QW+QR)
!             IF (DUM1 >  DUM2) THEN
!               WRITE(6,"(/2(a,i4),a,i2)") '{@ i=',I_index,' j=',J_index,
!     &                                     ' L=',L
!               WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{@ TCold=',TC,'P=',.01*PP,'DIFF=',DUM,'WCold=',WC,
!     & '{@ QIold=',QI,'QWold=',QW,'QRold=',QR
!             ENDIF
!
!***********************************************************************
!*********** MAIN MICROPHYSICS CALCULATIONS NOW FOLLOW! ****************
!***********************************************************************
!
!--- Calculate a few variables, which are used more than once below
!
!--- Latent heat of vaporization as a function of temperature from
!      Bolton (1980, JAS)
!
              XLV    = 3.148E6_r8 - 2370*TK     ! Latent heat of vaporization (Lv)
              XLF    = XLS-XLV               ! Latent heat of fusion (Lf)
              XLV1   = XLV*RCP               ! Lv/Cp
              XLF1   = XLF*RCP               ! Lf/Cp
              TK2    = 1.0_r8/(TK*TK)            ! 1./TK**2
              XLV2   = XLV*XLV*QSW_l*TK2/con_rv  ! Lv**2*Qsw_l/(Rv*TK**2)
              DENOMW = 1.0_r8 + XLV2*RCP         ! Denominator term, Clausius-Clapeyron correction
!
!--- Basic thermodynamic quantities
!      *      DYNVIS     - dynamic viscosity           [ kg/(m*s) ]
!      *      THERM_COND - thermal conductivity        [ J/(m*s*K) ]
!      *      DIFFUS     - diffusivity of water vapor  [ m**2/s ]
!
!             TFACTOR    = TK**1.5_r8/(TK+120.0_r8)
              TFACTOR    = TK*sqrt(TK)/(TK+120.0_r8)
              DYNVIS     = 1.496E-6_r8*TFACTOR
              THERM_COND = 2.116E-3_r8*TFACTOR
              DIFFUS     = 8.794E-5_r8*TK**1.81_r8/PP
              SCHMIT_FAC = (RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
!
!--- Air resistance term for the fall speed of ice following the
!      basic research by Heymsfield, Kajikawa, others 
!
              GAMMAS = (1.E5_r8/PP)**C1
!
!--- Air resistance for rain fall speed (Beard, 1985, JAOT, p.470)
!
              GAMMAR = (RHO0/RHO)**0.4_r8
!
!----------------------------------------------------------------------
!-------------  IMPORTANT MICROPHYSICS DECISION TREE  -----------------
!----------------------------------------------------------------------
!
!--- Determine if conditions supporting ice are present
!
              IF (TC < 0.0_r8 .OR. QI > EPSQ .OR. ASNOW > CLIMIT) THEN
                ICE_logical = .TRUE.
              ELSE
                ICE_logical = .FALSE.
                QLICE = 0.0_r8
                QTICE = 0.0_r8
              ENDIF
!
!--- Determine if rain is present
!
              RAIN_logical = .FALSE.
              IF (ARAIN > CLIMIT .OR. QR > EPSQ) RAIN_logical = .TRUE.
!
              IF (ICE_logical) THEN
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!---
!  * FLARGE  - ratio of number of large ice to total (large & small) ice
!  * FSMALL  - ratio of number of small ice crystals to large ice particles
!  ->  Small ice particles are assumed to have a mean diameter of 50 microns.
!  * XSIMASS - used for calculating small ice mixing ratio
!---
!  * TOT_ICE - total mass (small & large) ice before microphysics,
!              which is the sum of the total mass of large ice in the 
!              current layer and the input flux of ice from above
!  * PILOSS  - greatest loss (<0) of total (small & large) ice by
!              sublimation, removing all of the ice falling from above
!              and the ice within the layer
!  * RimeF1  - Rime Factor, which is the mass ratio of total (unrimed & rimed) 
!              ice mass to the unrimed ice mass (>=1)
!  * VrimeF  - the velocity increase due to rime factor or melting (ratio, >=1)
!  * VSNOW   - Fall speed of rimed snow w/ air resistance correction
!  * EMAIRI  - equivalent mass of air associated layer and with fall of snow into layer
!  * XLIMASS - used for calculating large ice mixing ratio
!  * FLIMASS - mass fraction of large ice
!  * QTICE   - time-averaged mixing ratio of total ice
!  * QLICE   - time-averaged mixing ratio of large ice
!  * NLICE   - time-averaged number concentration of large ice
!  * NSmICE  - number concentration of small ice crystals at current level
!---
!--- Assumed number fraction of large ice particles to total (large & small) 
!    ice particles, which is based on a general impression of the literature.
!
                WVQW = WV + QW                ! Water vapor + cloud water
!
!--- 6/19/03 - Deleted some code here ....
!
!  *********************************************************

!               IF (TC >= 0. .OR. WVQW < QSIgrd) THEN
!  !
!  !--- Eliminate small ice particle contributions for melting & sublimation
!  !
!                 FLARGE = FLARGE1
!               ELSE
!  !
!  !--- Enhanced number of small ice particles during depositional growth
!  !    (effective only when 0C > T >= T_ice [-10C] )
!  !
!                 FLARGE = FLARGE2
!  !
!  !--- Larger number of small ice particles due to rime splintering
!  !
!                 IF (TC >= -8. .AND. TC <= -3.) FLARGE=.5*FLARGE
!
!               ENDIF            ! End IF (TC >= 0. .OR. WVQW < QSIgrd)
!               FSMALL=(1.-FLARGE)/FLARGE
!               XSIMASS=RRHO*MASSI(MDImin)*FSMALL
!  *********************************************************
!
                IF (QI <= EPSQ .AND. ASNOW <= CLIMIT) THEN
                  INDEXS  = MDImin
                  FLARGE  = 0.0_r8                   !--- Begin 6/19/03 changes
                  FSMALL  = 1.0_r8
                  XSIMASS = RRHO*MASSI(MDImin)   !--- End 6/19/03 changes
                  TOT_ICE = 0.0_r8
                  PILOSS  = 0.0_r8
                  RimeF1  = 1.0_r8
                  VrimeF  = 1.0_r8
                  VEL_INC = GAMMAS
                  VSNOW   = 0.0_r8
                  EMAIRI  = THICK
                  XLIMASS = RRHO*RimeF1*MASSI(INDEXS)
                  FLIMASS = XLIMASS/(XLIMASS+XSIMASS)
                  QLICE   = 0.0_r8
                  QTICE   = 0.0_r8
                  NLICE   = 0.0_r8
                  NSmICE  = 0.0_r8
                ELSE
   !
   !--- For T<0C mean particle size follows Houze et al. (JAS, 1979, p. 160), 
   !    converted from Fig. 5 plot of LAMDAs.  Similar set of relationships 
   !    also shown in Fig. 8 of Ryan (BAMS, 1996, p. 66).
   !
   !--- Begin 6/19/03 changes => allow NLImax to increase & FLARGE to 
   !    decrease at COLDER temperatures; set FLARGE to zero (i.e., only small
   !    ice) if the ice mixing ratio is below QI_min

!                 DUM    = MAX(0.05_r8, MIN(1.0_r8, EXP(.0536_r8*TC)) )
                  DUM    = MAX(0.05_r8, MIN(1.0_r8, EXP(.0564_r8*TC)) )
                  INDEXS = MIN(MDImax, MAX(MDImin, INT(XMImax*DUM) ) )
!   Pressure dependen scaling factor for flgmin_l (tunable)
!
!   Pressure dependen scaling factor for flgmin_l (tunable)
!
!!!         pfac = max(0.5, (     min(1.0, pp*0.00002))**2)   ! commented on 02182011
!go         pfac = max(0.5, (sqrt(min(1.0, pp*0.00004))))
!           pfac = 1.0_r8

                  DUM    = MAX(flgmin_l*pfac, DUM)

                  QI_min = QI_min_0C * dum  !- Ver5    ----Not used ----
!!                QI_min = QI_min_0C        !- Ver5
!!!               QI_min = QI_min_0C/DUM    !- Ver5

                  NLImax = 10.E3_r8/sqrt(DUM)  !- Ver3
                  IF (TC < 0.0_r8) THEN
                    FLARGE = DUM            !- Ver4
                  ELSE
                    FLARGE = 1.0_r8
                  ENDIF
                  FSMALL  = (1.0_r8-FLARGE)/FLARGE
                  XSIMASS = RRHO*MASSI(MDImin)*FSMALL
                  TOT_ICE = THICK*QI + BLEND*ASNOW
                  PILOSS  = -TOT_ICE/THICK
                  LBEF    = MAX(1,L-1)
                  RimeF1  = (RimeF_col(L)*THICK*QI                      &
     &                    + RimeF_col(LBEF)*BLEND*ASNOW)/TOT_ICE
                  RimeF1  = MIN(RimeF1, RFmax)
                  VSNOW   = 0.0_r8
                  DO IPASS=0,1
                    IF (RimeF1 .LE. 1.0_r8) THEN
                      RimeF1 = 1.0_r8
                      VrimeF = 1.0_r8
                    ELSE
                      IXS  = MAX(2, MIN(INDEXS/100, 9))
                      XRF  = 10.492_r8*LOG(RimeF1)
                      IXRF = MAX(0, MIN(INT(XRF), Nrime))
                      IF (IXRF .GE. Nrime) THEN
                        VrimeF = VEL_RF(IXS,Nrime)
                      ELSE
                        VrimeF = VEL_RF(IXS,IXRF)+(XRF-FLOAT(IXRF))*    &
     &                          (VEL_RF(IXS,IXRF+1)-VEL_RF(IXS,IXRF))
                      ENDIF
                    ENDIF            ! End IF (RimeF1 <= 1.)
                    VEL_INC = GAMMAS*VrimeF
                    VSNOW   = VEL_INC*VSNOWI(INDEXS)
                    EMAIRI  = THICK + BLDTRH*VSNOW
                    XLIMASS = RRHO*RimeF1*MASSI(INDEXS)
                    FLIMASS = XLIMASS/(XLIMASS+XSIMASS)
                    QTICE   = TOT_ICE/EMAIRI
                    QLICE   = FLIMASS*QTICE
                    NLICE   = QLICE/XLIMASS
                    NSmICE  = Fsmall*NLICE
   !
                    IF ( (NLICE >= NLImin .AND. NLICE <= NLImax)        & 
     &                    .OR. IPASS == 1) THEN
                      EXIT
                    ELSE
                      IF(TC < 0) THEN
                        XLI = RHO*(QTICE/DUM-XSIMASS)/RimeF1
                        IF (XLI <= MASSI(MDImin) ) THEN
                          INDEXS = MDImin
                        ELSE IF (XLI <= MASSI(450) ) THEN
                          DLI    = 9.5885E5_r8*XLI**0.42066_r8       ! DLI in microns
                          INDEXS = MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                        ELSE IF (XLI <= MASSI(MDImax) ) THEN
                          DLI    = 3.9751E6_r8*XLI**0.49870_r8       ! DLI in microns
                          INDEXS = MIN(MDImax, MAX(MDImin, INT(DLI) ) )
                        ELSE 
                          INDEXS = MDImax
                        ENDIF             ! End IF (XLI <= MASSI(MDImin) ) 
                      ENDIF               ! End IF (TC < 0)
!
!--- Reduce excessive accumulation of ice at upper levels
!    associated with strong grid-resolved ascent
!
!--- Force NLICE to be between NLImin and NLImax
!
!--- 8/22/01: Increase density of large ice if maximum limits
!    are reached for number concentration (NLImax) and mean size
!    (MDImax).  Done to increase fall out of ice.
!
!

                      DUM = MAX(NLImin, MIN(NLImax, NLICE) )
                      IF (DUM >= NLImax .AND. INDEXS >= MDImax)         &
     &                 RimeF1 = RHO*(QTICE/NLImax-XSIMASS)/MASSI(INDEXS)
!
!                WRITE(6,"(4(a12,g11.4,1x))") 
!     & '{$ TC=',TC,'P=',.01*PP,'NLICE=',NLICE,'DUM=',DUM,
!     & '{$ XLI=',XLI,'INDEXS=',FLOAT(INDEXS),'RHO=',RHO,'QTICE=',QTICE,
!     & '{$ XSIMASS=',XSIMASS,'RimeF1=',RimeF1

                    ENDIF    ! End IF ( (NLICE >=NLImin .AND. NLICE >= NLImax)
                  ENDDO      ! End DO IPASS=0,1
                ENDIF        ! End IF (QI <= EPSQ .AND. ASNOW <= CLIMIT)
              ENDIF          ! End IF (ICE_logical)
!
!----------------------------------------------------------------------
!--------------- Calculate individual processes -----------------------
!----------------------------------------------------------------------
!
!--- Cloud water autoconversion to rain and collection by rain
!
              IF (QW > EPSQ .AND. TC >= T_ICE) THEN
   !
   !--- QW0 could be modified based on land/sea properties, 
   !      presence of convection, etc.  This is why QAUT0 and CRAUT
   !      are passed into the subroutine as externally determined
   !      parameters.  Can be changed in the future if desired.
   !
!               QW0   = QAUT0*RRHO
               ! PRINT*,QAUTx,RRHO,XNCW(L)
                QW0   = QAUTx*RRHO*XNCW(L)
                PRAUT = MAX(0.0_r8, QW-QW0)*CRAUT
                IF (QLICE  >  EPSQ) THEN
      !
      !--- Collection of cloud water by large ice particles ("snow")
      !    PIACWI=PIACW for riming, PIACWI=0 for shedding
      !
!Moor              FWS   = MIN(1.0_r8, CIACW*VEL_INC*NLICE*ACCRI(INDEXS)/PP**C1) ! 20050422
                   FWS   = MIN(0.1_r8, CIACW*VEL_INC*NLICE*ACCRI(INDEXS)/PP**C1)
                   PIACW = FWS*QW
                   IF (TC  < 0.0_r8) PIACWI = PIACW    ! Large ice riming

                ENDIF             ! End IF (QLICE > EPSQ)
              ENDIF               ! End IF (QW > EPSQ .AND. TC >= T_ICE)
!
!----------------------------------------------------------------------
!--- Loop around some of the ice-phase processes if no ice should be present
!----------------------------------------------------------------------
!
              IF (ICE_logical) THEN
!
!--- Now the pretzel logic of calculating ice deposition
!
                IF (TC < T_ICE .AND. (WV > QSIgrd .OR. QW > EPSQ)) THEN
!
!--- Adjust to ice saturation at T<T_ICE (-10C) if supersaturated.
!    Sources of ice due to nucleation and convective detrainment are
!    either poorly understood, poorly resolved at typical NWP 
!    resolutions, or are not represented (e.g., no detrained 
!    condensate in BMJ Cu scheme).
!    
                  PCOND = -QW
                  DUM1  = TK + XLV1*PCOND              ! Updated (dummy) temperature (deg K)
                  DUM2  = WV+QW                        ! Updated (dummy) water vapor mixing ratio
                  DUM   = min(pp,FPVSI(DUM1))          ! Updated (dummy) saturation vapor pressure w/r/t ice
                  DUM   = RHgrd*con_eps*DUM/(pp+con_epsm1*dum) ! Updated (dummy) saturation specific humidity w/r/t ice
!                 DUM   = RHgrd*EPS*DUM/(PP-DUM)       ! Updated (dummy) saturation mixing ratio w/r/t ice

                  IF (DUM2 > DUM) PIDEP = DEPOSIT(PP, RHgrd, DUM1, DUM2)

                  DWVi = 0.0_r8                            ! Used only for debugging
!
                ELSE IF (TC < 0.0_r8) THEN
!
!--- These quantities are handy for ice deposition/sublimation
!    PIDEP_max - max deposition or minimum sublimation to ice saturation
!
                  DENOMI    = 1.0_r8 + XLS2*QSI_l*TK2
                  DWVi      = MIN(WVQW,QSW_l)-QSI_l
                  PIDEP_max = MAX(PILOSS, DWVi/DENOMI)
                  IF (QTICE > 0.0_r8) THEN
!
!--- Calculate ice deposition/sublimation
!      * SFACTOR - [VEL_INC**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
!        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
!      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
!               VENTIL, VENTIS - m**-2 ;  VENTI1 - m ;  
!               VENTI2 - m**2/s**.5 ; DIDEP - unitless
!
!                   SFACTOR = VEL_INC**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                    SFACTOR = sqrt(VEL_INC)*SCHMIT_FAC
                    ABI     = 1.0_r8/(RHO*XLS3*QSI*TK2/THERM_COND+1.0_r8/DIFFUS)
!
!--- VENTIL - Number concentration * ventilation factors for large ice
!--- VENTIS - Number concentration * ventilation factors for small ice
!
!--- Variation in the number concentration of ice with time is not
!      accounted for in these calculations (could be in the future).
!
                    VENTIL = (VENTI1(INDEXS) + SFACTOR*VENTI2(INDEXS))  &
     &                                       * NLICE
                    VENTIS = (VENTI1(MDImin) + SFACTOR*VENTI2(MDImin))  &
     &                                       * NSmICE
                    DIDEP  = ABI*(VENTIL+VENTIS)*DTPH
!
!--- Account for change in water vapor supply w/ time
!
                    IF (DIDEP >= Xratio)                                & 
     &                DIDEP = (1.0_r8-EXP(-DIDEP*DENOMI))/DENOMI
                    IF (DWVi > 0.0_r8) THEN
                      PIDEP = MIN(DWVi*DIDEP, PIDEP_max)
                    ELSE IF (DWVi < 0.0_r8) THEN
                      PIDEP = MAX(DWVi*DIDEP, PIDEP_max)
                    ENDIF
!
                  ELSE IF (WVQW > QSI_l .AND. TC <= T_ICE_init) THEN
!
!--- Ice nucleation in near water-saturated conditions.  Ice crystal
!    growth during time step calculated using Miller & Young (1979, JAS).
!--- These deposition rates could drive conditions below water saturation,
!    which is the basis of these calculations.  Intended to approximate
!    more complex & computationally intensive calculations.
!
                    INDEX_MY = MAX(MY_T1, MIN( INT(0.5_r8-TC), MY_T2 ) )
!
!--- DUM1 is the supersaturation w/r/t ice at water-saturated conditions
!
!--- DUM2 is the number of ice crystals nucleated at water-saturated 
!    conditions based on Meyers et al. (JAM, 1992).
!
!--- Prevent unrealistically large ice initiation (limited by PIDEP_max)
!      if DUM2 values are increased in future experiments
!
                    DUM1  = QSW/QSI - 1.0_r8      
                    DUM2  = 1.E3_r8*EXP(12.96_r8*DUM1-0.639_r8)
                    PIDEP = MIN(PIDEP_max,DUM2*MY_GROWTH(INDEX_MY)*RRHO)
!
                  ENDIF ! End IF (QTICE > 0.0_r8)
!
                ENDIF   ! End IF (TC < T_ICE .AND. (WV > QSIgrd .OR. QW > EPSQ))
!
!------------------------------------------------------------------------
!
              ENDIF     ! End of IF (ICE_logical)then loop
!
!------------------------------------------------------------------------
!
!--- Cloud water condensation
!
              IF (TC >= T_ICE .AND. (QW > EPSQ .OR. WV > QSWgrd)) THEN
                IF (PIACWI == 0.0_r8 .AND. PIDEP == 0.0_r8) THEN
                  PCOND = CONDENSE (PP, QW, RHgrd, TK, WV)
                ELSE  !-- Modify cloud condensation in response to ice processes
                  DUM     = XLV*QSWgrd*RCPRV*TK2
                  DENOMWI = 1.0_r8 + XLS*DUM
                  DENOMF  = XLF*DUM
                  DUM     = MAX(0.0_r8, PIDEP)
                  PCOND   = (WV-QSWgrd-DENOMWI*DUM-DENOMF*PIACWI)/DENOMW
                  DUM1    = -QW
                  DUM2    = PCOND - PIACW
                  IF (DUM2  <  DUM1) THEN    !--- Limit cloud water sinks
                    DUM    = DUM1/DUM2
                    PCOND  = DUM*PCOND
                    PIACW  = DUM*PIACW
                    PIACWI = DUM*PIACWI
                  ENDIF ! EndIF (DUM2 <  DUM1)
                ENDIF   ! EndIF (PIACWI == 0. .AND. PIDEP == 0.)
              ENDIF     ! EndIF (TC >= T_ICE .AND. (QW > EPSQ .OR. WV > QSWgrd))
!
!--- Limit freezing of accreted rime to prevent temperature oscillations,
!    a crude Schumann-Ludlam limit (p. 209 of Young, 1993). 
!
              TCC = TC + XLV1*PCOND + XLS1*PIDEP + XLF1*PIACWI
              IF (TCC  >  0.0_r8) THEN
                PIACWI = 0.0_r8
                TCC    = TC + XLV1*PCOND + XLS1*PIDEP
              ENDIF
!
              IF (TC > 0.0_r8 .AND. TCC > 0. .AND. ICE_logical) THEN
!
!--- Calculate melting and evaporation/condensation
!      * Units: SFACTOR - s**.5/m ;  ABI - m**2/s ;  NLICE - m**-3 ;
!               VENTIL - m**-2 ;  VENTI1 - m ;  
!               VENTI2 - m**2/s**.5 ; CIEVP - /s
!
!               SFACTOR = VEL_INC**0.5_r8*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                SFACTOR = sqrt(VEL_INC)*SCHMIT_FAC
                VENTIL  = NLICE*(VENTI1(INDEXS)+SFACTOR*VENTI2(INDEXS))
                AIEVP   = VENTIL*DIFFUS*DTPH
                IF (AIEVP  <  Xratio) THEN
                  DIEVP = AIEVP
                ELSE
                  DIEVP = 1.0_r8 - EXP(-AIEVP)
                ENDIF
!               QSW0 = con_eps*ESW0/(PP-ESW0)
!               QSW0 = con_eps*ESW0/(PP+con_epsm1*ESW0)
!!              dum  = min(PP, ESW0)
!!              QSW0 = con_eps*dum/(PP+con_epsm1*dum)
!               DWV0 = MIN(WV,QSW)-QSW0
                DWV0 = MIN(WV,QSW_l)-QSW0_l
                DUM  = QW + PCOND
                IF (WV < QSW_l .AND. DUM <= EPSQ) THEN
   !
   !--- Evaporation from melting snow (sink of snow) or shedding
   !    of water condensed onto melting snow (source of rain)
   !
                  DUM   = DWV0*DIEVP
                  PIEVP = MAX( MIN(0.0_r8, DUM), PILOSS)
                  PICND = MAX(0.0_r8, DUM)
                ENDIF            ! End IF (WV < QSW_l .AND. DUM <= EPSQ)
                PIMLT = THERM_COND*TCC*VENTIL*RRHO*DTPH/XLF
   !
   !--- Limit melting to prevent temperature oscillations across 0C
   !
                DUM1  = MAX( 0.0_r8, (TCC+XLV1*PIEVP)/XLF1 )
                PIMLT = MIN(PIMLT, DUM1)
   !
   !--- Limit loss of snow by melting (>0) and evaporation
   !
                DUM = PIEVP - PIMLT
                IF (DUM < PILOSS) THEN
                  DUM1  = PILOSS/DUM
                  PIMLT = PIMLT*DUM1
                  PIEVP = PIEVP*DUM1
                ENDIF       ! End IF (DUM  > QTICE)
              ENDIF         ! End IF (TC > 0. .AND. TCC > 0. .AND. ICE_logical)
!
!--- IMPORTANT:  Estimate time-averaged properties.
!
!  * TOT_RAIN - total mass of rain before microphysics, which is the sum of
!               the total mass of rain in the current layer and the input 
!               flux of rain from above
!  * VRAIN1   - fall speed of rain into grid from above (with air resistance correction)
!  * QTRAIN   - time-averaged mixing ratio of rain (kg/kg)
!  * PRLOSS   - greatest loss (<0) of rain, removing all rain falling from
!               above and the rain within the layer
!  * RQR      - rain content (kg/m**3)
!  * INDEXR   - mean size of rain drops to the nearest 1 micron in size
!  * N0r      - intercept of rain size distribution (typically 10**6 m**-4)
!
              TOT_RAIN = 0.0_r8
              VRAIN1   = 0.0_r8
              QTRAIN   = 0.0_r8
              PRLOSS   = 0.0_r8
              RQR      = 0.0_r8
              N0r      = 0.0_r8
              INDEXR1  = INDEXR    ! For debugging only
              INDEXR   = MDRmin
              IF (RAIN_logical) THEN
                IF (ARAIN <= 0.0_r8) THEN
                  INDEXR = MDRmin
                  VRAIN1 = 0.0_r8
                ELSE
   !
   !--- INDEXR (related to mean diameter) & N0r could be modified 
   !      by land/sea properties, presence of convection, etc.
   !
   !--- Rain rate normalized to a density of 1.194 kg/m**3
   !
                  RR = ARAIN / (DTPH*GAMMAR)
   !
                  IF (RR <= RR_DRmin) THEN
        !
        !--- Assume fixed mean diameter of rain (0.2 mm) for low rain rates, 
        !      instead vary N0r with rain rate
        !
                    INDEXR = MDRmin
                  ELSE IF (RR <= RR_DR1) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.05 and 0.10 mm:
        !      V(Dr)=5.6023e4*Dr**1.136, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*5.6023e4*Dr**(4+1.136) = 1.408e15*Dr**5.136,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.123e-3*RR**.1947 -> Dr (microns) = 1.123e3*RR**.1947
        !
                    INDEXR = INT( 1.123E3_r8*RR**0.1947_r8 + 0.5_r8 )
                    INDEXR = MAX( MDRmin, MIN(INDEXR, MDR1) )

                  ELSE IF (RR <= RR_DR2) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.10 and 0.20 mm:
        !      V(Dr)=1.0867e4*Dr**.958, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*1.0867e4*Dr**(4+.958) = 2.731e14*Dr**4.958,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.225e-3*RR**.2017 -> Dr (microns) = 1.225e3*RR**.2017
        !
                    INDEXR = INT( 1.225E3_r8*RR**0.2017_r8 + 0.5_r8 )
                    INDEXR = MAX( MDR1, MIN(INDEXR, MDR2) )

                  ELSE IF (RR <= RR_DR3) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.20 and 0.32 mm:
        !      V(Dr)=2831.*Dr**.80, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*2831.*Dr**(4+.80) = 7.115e13*Dr**4.80, 
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.3006e-3*RR**.2083 -> Dr (microns) = 1.3006e3*RR**.2083
        !
                    INDEXR = INT( 1.3006E3_r8*RR**0.2083_r8 + 0.5_r8 )
                    INDEXR = MAX( MDR2, MIN(INDEXR, MDR3) )

                  ELSE IF (RR <= RR_DRmax) THEN
        !
        !--- Best fit to mass-weighted fall speeds (V) from rain lookup tables 
        !      for mean diameters (Dr) between 0.32 and 0.45 mm:
        !      V(Dr)=944.8*Dr**.6636, V in m/s and Dr in m
        !      RR = PI*1000.*N0r0*944.8*Dr**(4+.6636) = 2.3745e13*Dr**4.6636,
        !        RR in kg/(m**2*s)
        !      Dr (m) = 1.355e-3*RR**.2144 -> Dr (microns) = 1.355e3*RR**.2144
        !
                    INDEXR = INT( 1.355E3_r8*RR**0.2144_r8 + 0.5_r8 )
                    INDEXR = MAX( MDR3, MIN(INDEXR, MDRmax) )
                  ELSE 
        !
        !--- Assume fixed mean diameter of rain (0.45 mm) for high rain rates, 
        !      instead vary N0r with rain rate
        !
                    INDEXR = MDRmax
                  ENDIF               ! End IF (RR <= RR_DRmin) etc. 
!
                  VRAIN1 = GAMMAR*VRAIN(INDEXR)
                ENDIF                 ! End IF (ARAIN <= 0.)
!
                INDEXR1  = INDEXR     ! For debugging only
                TOT_RAIN = THICK*QR+BLEND*ARAIN
                QTRAIN   = TOT_RAIN/(THICK+BLDTRH*VRAIN1)
                PRLOSS   = -TOT_RAIN/THICK
                RQR      = RHO*QTRAIN
   !
   !--- RQR - time-averaged rain content (kg/m**3)
   !
                IF (RQR <= RQR_DRmin) THEN
                  N0r    = MAX(N0rmin, CN0r_DMRmin*RQR)
                  INDEXR = MDRmin
                ELSE IF (RQR >= RQR_DRmax) THEN
                  N0r    = CN0r_DMRmax*RQR
                  INDEXR = MDRmax
                ELSE
                  N0r    = N0r0
!                 INDEXR = MAX( XMRmin, MIN(CN0r0*RQR**.25, XMRmax) )
                  item   = INT(CN0r0*sqrt(sqrt(RQR)) )              ! Moorthi 07/31/08
                  INDEXR = MAX( MDRmin, MIN(item, MDRmax) )         ! Moorthi 07/31/08
                ENDIF
   !
                IF (TC < T_ICE) THEN
                  PIACR = -PRLOSS
                ELSE
                  DWVr = WV - PCOND - QSW_l
                  DUM  = QW + PCOND
                  IF (DWVr < 0.0_r8 .AND. DUM <= EPSQ) THEN
!
!--- Rain evaporation
!
!    * RFACTOR - [GAMMAR**.5]*[Schmidt**(1./3.)]*[(RHO/DYNVIS)**.5],
!        where Schmidt (Schmidt Number) =DYNVIS/(RHO*DIFFUS)
!
!    * Units: RFACTOR - s**.5/m ;  ABW - m**2/s ;  VENTR - m**-2 ;  
!             N0r - m**-4 ;  VENTR1 - m**2 ;  VENTR2 - m**3/s**.5 ;
!             CREVP - unitless
!
!                   RFACTOR = GAMMAR**.5*(RHO/(DIFFUS*DIFFUS*DYNVIS))**C2
                    RFACTOR = sqrt(GAMMAR)*SCHMIT_FAC
                    ABW     = 1.0_r8/(RHO*XLV2/THERM_COND+1.0_r8/DIFFUS)
!
!--- Note that VENTR1, VENTR2 lookup tables do not include the 
!      1/Davg multiplier as in the ice tables
!
                    VENTR = N0r*(VENTR1(INDEXR)+RFACTOR*VENTR2(INDEXR))
                    CREVP = ABW*VENTR*DTPH
                    IF (CREVP < Xratio) THEN
                      DUM = DWVr*CREVP
                    ELSE
                      DUM = DWVr*(1.0_r8-EXP(-CREVP*DENOMW))/DENOMW
                    ENDIF
                    PREVP = MAX(DUM, PRLOSS)
                  ELSE IF (QW > EPSQ) THEN
                    FWR   = CRACW*GAMMAR*N0r*ACCRR(INDEXR)
!Moor               PRACW = MIN(1.0_r8,FWR)*QW               ! 20050422
                     PRACW = MIN(0.1_r8,FWR)*QW

                  ENDIF           ! End IF (DWVr < 0._r8 .AND. DUM <= EPSQ)
!
                  IF (TC < 0.0_r8 .AND. TCC < 0.0_r8) THEN
!
!--- Biggs (1953) heteorogeneous freezing (e.g., Lin et al., 1983)
!   - Rescaled mean drop diameter from microns (INDEXR) to mm (DUM) to prevent underflow
!
                    DUM   = 0.001_r8*FLOAT(INDEXR)
                    dum1  = dum * dum
                    DUM   = (EXP(ABFR*TC)-1.0_r8)*DUM1*DUM1*DUM1*DUM
                    PIACR = MIN(CBFR*N0r*RRHO*DUM, QTRAIN)
                    IF (QLICE > EPSQ) THEN
            !
            !--- Freezing of rain by collisions w/ large ice
            !
                      DUM  = GAMMAR*VRAIN(INDEXR)
                      DUM1 = DUM-VSNOW
            !
            !--- DUM2 - Difference in spectral fall speeds of rain and
            !      large ice, parameterized following eq. (48) on p. 112 of 
            !      Murakami (J. Meteor. Soc. Japan, 1990)
            !
                      DUM2 = (DUM1*DUM1+0.04_r8*DUM*VSNOW)**0.5_r8
                      DUM1 = 5.E-12_r8*INDEXR*INDEXR+2.E-12_r8*INDEXR*INDEXS      &
     &                     +0.5E-12_r8*INDEXS*INDEXS
                      FIR = MIN(1.0_r8, CIACR*NLICE*DUM1*DUM2)
            !
            !--- Future?  Should COLLECTION BY SMALL ICE SHOULD BE INCLUDED???
            !
                      PIACR = MIN(PIACR+FIR*QTRAIN, QTRAIN)
                    ENDIF        ! End IF (QLICE >  EPSQ)
                    DUM = PREVP - PIACR
                    If (DUM < PRLOSS) THEN
                      DUM1  = PRLOSS/DUM
                      PREVP = DUM1*PREVP
                      PIACR = DUM1*PIACR
                    ENDIF        ! End If (DUM < PRLOSS)
                  ENDIF          ! End IF (TC < 0.0_r8 .AND. TCC < 0.0_r8)
                ENDIF            ! End IF (TC < T_ICE)
              ENDIF              ! End IF (RAIN_logical) 
!
!----------------------------------------------------------------------
!---------------------- Main Budget Equations -------------------------
!----------------------------------------------------------------------
!
!
!-----------------------------------------------------------------------
!--- Update fields, determine characteristics for next lower layer ----
!-----------------------------------------------------------------------
!
!--- Carefully limit sinks of cloud water
!
              DUM1 = PIACW + PRAUT + PRACW - MIN(0.0_r8,PCOND)
              IF (DUM1 > QW) THEN
                DUM    = QW/DUM1
                PIACW  = DUM*PIACW
                PIACWI = DUM*PIACWI
                PRAUT  = DUM*PRAUT
                PRACW  = DUM*PRACW
                IF (PCOND < 0.0_r8) PCOND=DUM*PCOND
              ENDIF
              PIACWR = PIACW - PIACWI          ! TC >= 0C
!
!--- QWnew - updated cloud water mixing ratio
!
              DELW  = PCOND - PIACW - PRAUT - PRACW
              QWnew = QW+DELW
              IF (QWnew <=  EPSQ) QWnew = 0.0_r8
              IF (QW > 0.0_r8 .AND. QWnew /= 0.0_r8) THEN
                DUM = QWnew/QW
                IF (DUM  < TOLER) QWnew = 0.0_r8
              ENDIF
!
!--- Update temperature and water vapor mixing ratios
!
              DELT = XLV1 * (PCOND+PIEVP+PICND+PREVP)                   &
     &             + XLS1 * PIDEP + XLF1*(PIACWI+PIACR-PIMLT)
              Tnew = TK + DELT
!
              DELV  = -PCOND - PIDEP - PIEVP - PICND - PREVP
              WVnew = WV + DELV
!
!--- Update ice mixing ratios
!
!---
!  * TOT_ICEnew - total mass (small & large) ice after microphysics,
!                 which is the sum of the total mass of large ice in the 
!                 current layer and the flux of ice out of the grid box below
!  * RimeF      - Rime Factor, which is the mass ratio of total (unrimed & 
!                 rimed) ice mass to the unrimed ice mass (>=1)
!  * QInew      - updated mixing ratio of total (large & small) ice in layer
!      -> TOT_ICEnew=QInew*THICK+BLDTRH*QLICEnew*VSNOW
!        -> But QLICEnew=QInew*FLIMASS, so
!      -> TOT_ICEnew=QInew*(THICK+BLDTRH*FLIMASS*VSNOW)
!  * ASNOWnew   - updated accumulation of snow at bottom of grid cell
!---
!
              DELI  = 0.0_r8
              RimeF = 1.0_r8
              IF (ICE_logical) THEN
                DELI       = PIDEP + PIEVP + PIACWI + PIACR - PIMLT
                TOT_ICEnew = TOT_ICE + THICK*DELI
                IF (TOT_ICE > 0.0_r8 .AND. TOT_ICEnew /= 0.0_r8) THEN
                  DUM = TOT_ICEnew/TOT_ICE
                  IF (DUM  < TOLER) TOT_ICEnew = 0.0_r8
                ENDIF
                IF (TOT_ICEnew <= CLIMIT) THEN
                  TOT_ICEnew = 0.0_r8
                  RimeF      = 1.0_r8
                  QInew      = 0.0_r8
                  ASNOWnew   = 0.0_r8
                ELSE
      !
      !--- Update rime factor if appropriate
      !
                  DUM = PIACWI + PIACR
                  IF (DUM <= EPSQ .AND. PIDEP <= EPSQ) THEN
                    RimeF = RimeF1
                  ELSE
         !
         !--- Rime Factor, RimeF = (Total ice mass)/(Total unrimed ice mass)
         !      DUM1 - Total ice mass, rimed & unrimed
         !      DUM2 - Estimated mass of *unrimed* ice
         !
                    DUM1 = TOT_ICE + THICK*(PIDEP+DUM)
                    DUM2 = TOT_ICE/RimeF1 + THICK*PIDEP
                    IF (DUM2 <= 0.0_r8) THEN
                      RimeF = RFmax
                    ELSE
                      RimeF = MIN(RFmax, MAX(1.0_r8, DUM1/DUM2) )
                    ENDIF
                  ENDIF       ! End IF (DUM <= EPSQ .AND. PIDEP <= EPSQ)
                  QInew = TOT_ICEnew/(THICK+BLDTRH*FLIMASS*VSNOW)
                  IF (QInew  <=  EPSQ) QInew = 0.0_r8
                  IF (QI > 0.0_r8 .AND. QInew /= 0.0_r8) THEN
                    DUM = QInew/QI
                    IF (DUM < TOLER) QInew = 0.0_r8
                  ENDIF
                  ASNOWnew = BLDTRH*FLIMASS*VSNOW*QInew
                  IF (ASNOW > 0.0_r8 .AND. ASNOWnew /= 0.0_r8) THEN
                    DUM = ASNOWnew/ASNOW
                    IF (DUM < TOLER) ASNOWnew = 0.0_r8
                  ENDIF
                ENDIF         ! End IF (TOT_ICEnew <= CLIMIT)
              ENDIF           ! End IF (ICE_logical)
!
!--- Update rain mixing ratios
!
!---
! * TOT_RAINnew - total mass of rain after microphysics
!                 current layer and the input flux of ice from above
! * VRAIN2      - time-averaged fall speed of rain in grid and below 
!                 (with air resistance correction)
! * QRnew       - updated rain mixing ratio in layer
!      -> TOT_RAINnew=QRnew*(THICK+BLDTRH*VRAIN2)
!  * ARAINnew  - updated accumulation of rain at bottom of grid cell
!---
!
              DELR        = PRAUT + PRACW + PIACWR - PIACR + PIMLT      &
     &                    + PREVP + PICND
              TOT_RAINnew = TOT_RAIN+THICK*DELR
              IF (TOT_RAIN > 0.0_r8 .AND. TOT_RAINnew /= 0.0_r8) THEN
                DUM = TOT_RAINnew/TOT_RAIN
                IF (DUM < TOLER) TOT_RAINnew = 0.0_r8
              ENDIF
              IF (TOT_RAINnew <= CLIMIT) THEN
                TOT_RAINnew = 0.0_r8
                VRAIN2      = 0.0_r8
                QRnew       = 0.0_r8
                ARAINnew    = 0.0_r8
              ELSE
   !
   !--- 1st guess time-averaged rain rate at bottom of grid box
   !
                RR = TOT_RAINnew/(DTPH*GAMMAR)
   !
   !--- Use same algorithm as above for calculating mean drop diameter
   !      (IDR, in microns), which is used to estimate the time-averaged
   !      fall speed of rain drops at the bottom of the grid layer.  This
   !      isn't perfect, but the alternative is solving a transcendental 
   !      equation that is numerically inefficient and nasty to program
   !      (coded in earlier versions of GSMCOLUMN prior to 8-22-01).
   !
                IF (RR <= RR_DRmin) THEN
                  IDR = MDRmin
                ELSE IF (RR <= RR_DR1) THEN
                  IDR = INT( 1.123E3_r8*RR**0.1947_r8 + 0.5_r8 )
                  IDR = MAX( MDRmin, MIN(IDR, MDR1) )
                ELSE IF (RR <= RR_DR2) THEN
                  IDR = INT( 1.225E3_r8*RR**0.2017_r8 + 0.5_r8 )
                  IDR = MAX( MDR1, MIN(IDR, MDR2) )
                ELSE IF (RR <= RR_DR3) THEN
                  IDR = INT( 1.3006E3_r8*RR**0.2083_r8 + 0.5_r8 )
                  IDR = MAX( MDR2, MIN(IDR, MDR3) )
                ELSE IF (RR <= RR_DRmax) THEN
                  IDR = INT( 1.355E3_r8*RR**0.2144_r8 + 0.5_r8 )
                  IDR = MAX( MDR3, MIN(IDR, MDRmax) )
                ELSE 
                  IDR = MDRmax
                ENDIF              ! End IF (RR <= RR_DRmin)
                VRAIN2 = GAMMAR*VRAIN(IDR)
                QRnew  = TOT_RAINnew / (THICK+BLDTRH*VRAIN2)
                IF (QRnew <= EPSQ) QRnew = 0.0_r8
                IF (QR > 0.0_r8 .AND. QRnew /= 0.0_r8) THEN
                  DUM = QRnew / QR
                  IF (DUM < TOLER) QRnew = 0.0_r8
                ENDIF
                ARAINnew = BLDTRH*VRAIN2*QRnew
                IF (ARAIN > 0.0_r8 .AND. ARAINnew /= 0.0_r8) THEN
                  DUM = ARAINnew/ARAIN
                  IF (DUM < TOLER) ARAINnew = 0.0_r8
                ENDIF
              ENDIF                ! End IF (TOT_RAINnew < CLIMIT)
!
              WCnew = QWnew + QRnew + QInew
!
!----------------------------------------------------------------------
!-------------- Begin debugging & verification ------------------------
!----------------------------------------------------------------------
!
!--- QT, QTnew - total water (vapor & condensate) before & after microphysics, resp.
!
!             QT=THICK*(QV+WC_col(l))+ARAIN+ASNOW
!             QTnew  = THICK*(WVnew/(1.+WVnew)+WCnew/(1.+wcnew))
!    &               + ARAINnew + ASNOWnew

              QT     = THICK*(WV+WC)       + ARAIN    + ASNOW
              QTnew  = THICK*(WVnew+WCnew) + ARAINnew + ASNOWnew
              BUDGET = QT-QTnew
!
!--- Additional check on budget preservation, accounting for truncation effects
!
              DBG_logical=.FALSE.
!             DUM=ABS(BUDGET)
!             IF (DUM > TOLER) THEN
!               DUM=DUM/MIN(QT, QTnew)
!               IF (DUM > TOLER) DBG_logical=.TRUE.
!             ENDIF
!
!             DUM=(RHgrd+.001)*QSInew
!             IF ( (QWnew > EPSQ .OR. QRnew > EPSQ .OR. WVnew > DUM)
!     &           .AND. TC < T_ICE )  DBG_logical=.TRUE.
!
!             IF (TC > 5. .AND. QInewr > EPSQ) DBG_logical=.TRUE.
!
!              IF ((WVnew < EPSQ .OR. DBG_logical) .AND. PRINT_diag) THEN
!!
!                WRITE(6,"(/2(a,i4),2(a,i2))") '{} i=',I_index,' j=',    &
!     &                                J_index, ' L=',L,' LSFC=',LSFC
!!
!                ESW    = min(PP, FPVSL(Tnew))
!!               QSWnew = EPS*ESW/(PP-ESW)
!                QSWnew = EPS*ESW/(PP+con_epsm1*ESW)
!                IF (TC < 0. .OR. Tnew  <  0.) THEN
!                  ESI    = min(PP, FPVSI(Tnew))
!!                 QSInew = EPS*ESI/(PP-ESI)
!                  QSInew = EPS*ESI/(PP+con_epsm1*ESI)
!                ELSE
!                  QSI    = QSW
!                  QSInew = QSWnew
!                ENDIF
!                WSnew = QSInew
!                WRITE(6,"(4(a12,g11.4,1x))")                            &
!     & '{} TCold=',TC,'TCnew=',Tnew-T0C,'P=',.01*PP,'RHO=',RHO,         &
!     & '{} THICK=',THICK,'RHold=',WV/WS,'RHnew=',WVnew/WSnew,           &
!     &   'RHgrd=',RHgrd,                                                &
!     & '{} RHWold=',WV/QSW,'RHWnew=',WVnew/QSWnew,'RHIold=',WV/QSI,     &
!     &   'RHInew=',WVnew/QSInew,                                        &
!     & '{} QSWold=',QSW,'QSWnew=',QSWnew,'QSIold=',QSI,'QSInew=',QSInew,&
!     & '{} WSold=',WS,'WSnew=',WSnew,'WVold=',WV,'WVnew=',WVnew,        &
!     & '{} WCold=',WC,'WCnew=',WCnew,'QWold=',QW,'QWnew=',QWnew,        &
!     & '{} QIold=',QI,'QInew=',QInew,'QRold=',QR,'QRnew=',QRnew,        &
!     & '{} ARAINold=',ARAIN,'ARAINnew=',ARAINnew,'ASNOWold=',ASNOW,     &
!     &   'ASNOWnew=',ASNOWnew,                                          &
!     & '{} TOT_RAIN=',TOT_RAIN,'TOT_RAINnew=',TOT_RAINnew,              &
!     &   'TOT_ICE=',TOT_ICE,'TOT_ICEnew=',TOT_ICEnew,                   &
!     & '{} BUDGET=',BUDGET,'QTold=',QT,'QTnew=',QTnew
!!
!                WRITE(6,"(4(a12,g11.4,1x))")                            &
!     & '{} DELT=',DELT,'DELV=',DELV,'DELW=',DELW,'DELI=',DELI,          &
!     & '{} DELR=',DELR,'PCOND=',PCOND,'PIDEP=',PIDEP,'PIEVP=',PIEVP,    &
!     & '{} PICND=',PICND,'PREVP=',PREVP,'PRAUT=',PRAUT,'PRACW=',PRACW,  &
!     & '{} PIACW=',PIACW,'PIACWI=',PIACWI,'PIACWR=',PIACWR,'PIMLT=',    &
!     &    PIMLT,                                                        &
!     & '{} PIACR=',PIACR
!!
!                IF (ICE_logical) WRITE(6,"(4(a12,g11.4,1x))")           &
!     & '{} RimeF1=',RimeF1,'GAMMAS=',GAMMAS,'VrimeF=',VrimeF,           &
!     &   'VSNOW=',VSNOW,                                                &
!     & '{} INDEXS=',FLOAT(INDEXS),'FLARGE=',FLARGE,'FSMALL=',FSMALL,    &
!     &   'FLIMASS=',FLIMASS,                                            &
!     & '{} XSIMASS=',XSIMASS,'XLIMASS=',XLIMASS,'QLICE=',QLICE,         &
!     &   'QTICE=',QTICE,                                                &
!     & '{} NLICE=',NLICE,'NSmICE=',NSmICE,'PILOSS=',PILOSS,             &
!     &   'EMAIRI=',EMAIRI,                                              &
!     & '{} RimeF=',RimeF
!!
!                IF (TOT_RAIN > 0. .OR. TOT_RAINnew > 0.)                &
!     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
!     & '{} INDEXR1=',FLOAT(INDEXR1),'INDEXR=',FLOAT(INDEXR),            &
!     &   'GAMMAR=',GAMMAR,'N0r=',N0r,                                   &
!     & '{} VRAIN1=',VRAIN1,'VRAIN2=',VRAIN2,'QTRAIN=',QTRAIN,'RQR=',RQR,&
!     & '{} PRLOSS=',PRLOSS,'VOLR1=',THICK+BLDTRH*VRAIN1,                &
!     &   'VOLR2=',THICK+BLDTRH*VRAIN2
!!
!                IF (PRAUT > 0.) WRITE(6,"(a12,g11.4,1x)") '{} QW0=',QW0
!!
!                IF (PRACW > 0.) WRITE(6,"(a12,g11.4,1x)") '{} FWR=',FWR
!!
!                IF (PIACR > 0.) WRITE(6,"(a12,g11.4,1x)") '{} FIR=',FIR
!!
!                DUM = PIMLT + PICND - PREVP - PIEVP
!                IF (DUM > 0. .or. DWVi /= 0.)                           &
!     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
!     & '{} TFACTOR=',TFACTOR,'DYNVIS=',DYNVIS,                          &
!     &   'THERM_CON=',THERM_COND,'DIFFUS=',DIFFUS
!!
!                IF (PREVP < 0.) WRITE(6,"(4(a12,g11.4,1x))")            &
!     & '{} RFACTOR=',RFACTOR,'ABW=',ABW,'VENTR=',VENTR,'CREVP=',CREVP,  &
!     & '{} DWVr=',DWVr,'DENOMW=',DENOMW
!!
!                IF (PIDEP /= 0. .AND. DWVi /= 0.)                       &
!     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
!     & '{} DWVi=',DWVi,'DENOMI=',DENOMI,'PIDEP_max=',PIDEP_max,         &
!     &   'SFACTOR=',SFACTOR,                                            &
!     & '{} ABI=',ABI,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),        &
!     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                             &
!     & '{} VENTIS=',VENTIS,'DIDEP=',DIDEP
!!
!                IF (PIDEP > 0. .AND. PCOND /= 0.)                       &
!     &            WRITE(6,"(4(a12,g11.4,1x))")                          &
!     & '{} DENOMW=',DENOMW,'DENOMWI=',DENOMWI,'DENOMF=',DENOMF,         &
!     &    'DUM2=',PCOND-PIACW
!!
!                IF (FWS > 0.) WRITE(6,"(4(a12,g11.4,1x))") '{} FWS=',FWS
!!
!                DUM = PIMLT + PICND - PIEVP
!                IF (DUM >  0.) WRITE(6,"(4(a12,g11.4,1x))")             &
!     & '{} SFACTOR=',SFACTOR,'VENTIL=',VENTIL,'VENTIL1=',VENTI1(INDEXS),&
!     &   'VENTIL2=',SFACTOR*VENTI2(INDEXS),                             &
!     & '{} AIEVP=',AIEVP,'DIEVP=',DIEVP,'QSW0=',QSW0,'DWV0=',DWV0
!   !
!              ENDIF
!!
!!----------------------------------------------------------------------
!!-------------- Water budget statistics & maximum values --------------
!!----------------------------------------------------------------------
!!
!              IF (PRINT_diag) THEN
!                ITdx = MAX( ITLO, MIN( INT(Tnew-T0C), ITHI ) )
!                IF (QInew > EPSQ) NSTATS(ITdx,1) = NSTATS(ITdx,1)+1
!                IF (QInew > EPSQ .AND.  QRnew+QWnew > EPSQ)             &
!     &            NSTATS(ITdx,2) = NSTATS(ITdx,2)+1
!                IF (QWnew > EPSQ) NSTATS(ITdx,3) = NSTATS(ITdx,3)+1 
!                IF (QRnew > EPSQ) NSTATS(ITdx,4) = NSTATS(ITdx,4)+1
!  !
!                QMAX(ITdx,1)  = MAX(QMAX(ITdx,1), QInew)
!                QMAX(ITdx,2)  = MAX(QMAX(ITdx,2), QWnew)
!                QMAX(ITdx,3)  = MAX(QMAX(ITdx,3), QRnew)
!                QMAX(ITdx,4)  = MAX(QMAX(ITdx,4), ASNOWnew)
!                QMAX(ITdx,5)  = MAX(QMAX(ITdx,5), ARAINnew)
!                QTOT(ITdx,1)  = QTOT(ITdx,1)+QInew*THICK
!                QTOT(ITdx,2)  = QTOT(ITdx,2)+QWnew*THICK
!                QTOT(ITdx,3)  = QTOT(ITdx,3)+QRnew*THICK
!  !
!                QTOT(ITdx,4)  = QTOT(ITdx,4)+PCOND*THICK
!                QTOT(ITdx,5)  = QTOT(ITdx,5)+PICND*THICK
!                QTOT(ITdx,6)  = QTOT(ITdx,6)+PIEVP*THICK
!                QTOT(ITdx,7)  = QTOT(ITdx,7)+PIDEP*THICK
!                QTOT(ITdx,8)  = QTOT(ITdx,8)+PREVP*THICK
!                QTOT(ITdx,9)  = QTOT(ITdx,9)+PRAUT*THICK
!                QTOT(ITdx,10) = QTOT(ITdx,10)+PRACW*THICK
!                QTOT(ITdx,11) = QTOT(ITdx,11)+PIMLT*THICK
!                QTOT(ITdx,12) = QTOT(ITdx,12)+PIACW*THICK
!                QTOT(ITdx,13) = QTOT(ITdx,13)+PIACWI*THICK
!                QTOT(ITdx,14) = QTOT(ITdx,14)+PIACWR*THICK
!                QTOT(ITdx,15) = QTOT(ITdx,15)+PIACR*THICK
!  !
!                QTOT(ITdx,16) = QTOT(ITdx,16)+(WVnew-WV)*THICK
!                QTOT(ITdx,17) = QTOT(ITdx,17)+(QWnew-QW)*THICK
!                QTOT(ITdx,18) = QTOT(ITdx,18)+(QInew-QI)*THICK
!                QTOT(ITdx,19) = QTOT(ITdx,19)+(QRnew-QR)*THICK
!                QTOT(ITdx,20) = QTOT(ITdx,20)+(ARAINnew-ARAIN)
!                QTOT(ITdx,21) = QTOT(ITdx,21)+(ASNOWnew-ASNOW)
!                IF (QInew > 0.)                                         &
!     &            QTOT(ITdx,22) = QTOT(ITdx,22)+QInew*THICK/RimeF
!  !
!              ENDIF
!
!----------------------------------------------------------------------
!------------------------- Update arrays ------------------------------
!----------------------------------------------------------------------
!
              T_col(L)     = Tnew                        ! temperature
!
!             QV_col(L)    = max(EPSQ, WVnew/(1.+WVnew)) ! specific humidity
              QV_col(L)    = max(0.0_r8, WVnew)          ! specific humidity
              WC_col(L)    = max(0.0_r8, WCnew)          ! total condensate mixing ratio
              QI_col(L)    = max(0.0_r8, QInew)          ! ice mixing ratio
              QR_col(L)    = max(0.0_r8, QRnew)          ! rain mixing ratio
              QW_col(L)    = max(0.0_r8, QWnew)          ! cloud water mixing ratio
              RimeF_col(L) = RimeF                       ! rime factor
              ASNOW        = ASNOWnew                    ! accumulated snow
              ARAIN        = ARAINnew                    ! accumulated rain
!
!#######################################################################
!
            ENDIF   ! End of IF (.NOT. CLEAR) THEN
          ENDIF     ! End of IF (QV_col(L) <= EPSQ .AND. WC_col(L) <= EPSQ) THEN
!
        ENDDO       ! ##### End "L" loop through model levels #####
!
        ARAING = ARAING + ARAIN
        ASNOWG = ASNOWG + ASNOW
      enddo              ! do for ntimes=1,mic_step
!
!#######################################################################
!
!-----------------------------------------------------------------------
!--------------------------- Return to GSMDRIVE -----------------------
!-----------------------------------------------------------------------
!
      CONTAINS
!     END  SUBROUTINE GSMCOLUMN
!
!#######################################################################
!--------- Produces accurate calculation of cloud condensation ---------
!#######################################################################
!
      real(kind=r8) FUNCTION CONDENSE (PP, QW, RHgrd, TK, WV)
!
      implicit none
!
!---------------------------------------------------------------------------------
!------ The Asai (1965) algorithm takes into consideration the release of ------
!------ latent heat in increasing the temperature & in increasing the     ------
!------ saturation mixing ratio (following the Clausius-Clapeyron eqn.).  ------
!---------------------------------------------------------------------------------
!
      real(kind=r8) pp, qw, rhgrd, tk, wv
!     INTEGER, PARAMETER :: HIGH_PRES=Selected_Real_Kind(15)
      real(kind=r8), PARAMETER ::                               &
     & RHLIMIT=0.001_r8, RHLIMIT1=-RHLIMIT
      real(kind=r8), PARAMETER :: RCP=1.0_r8/con_cp, RCPRV=RCP/con_rv
      real(kind=r8) :: COND, SSAT, WCdum
      real(kind=r8) :: wvdum, tdum, xlv, xlv1, xlv2, ws, dwv, esw, rfac
!
!-----------------------------------------------------------------------
!
!--- LV (T) is from Bolton (JAS, 1980)
!
!     XLV=3.148E6_r8-2370.0_r8*TK
!     XLV1=XLV*RCP
!     XLV2=XLV*XLV*RCPRV
!
      Tdum     = TK
      WVdum    = WV
      WCdum    = QW
      ESW      = min(PP, FPVSL(Tdum))          ! Saturation vapor press w/r/t water
      WS       = RHgrd*con_eps*ESW/(PP+con_epsm1*ESW)  ! Saturation specific hum w/r/t water
      DWV      = WVdum - WS                    ! Deficit grid-scale specific humidity
      SSAT     = DWV / WS                      ! Supersaturation ratio
      CONDENSE = 0.0_r8
      rfac     = 0.5_r8                           ! converges faster with 0.5_r8
      DO WHILE ((SSAT < RHLIMIT1 .AND. WCdum > EPSQ)  .OR. SSAT > RHLIMIT)
!
        XLV   = 3.148E6_r8-2370.0_r8*Tdum
        XLV1  = XLV*RCP
        XLV2  = XLV*XLV*RCPRV
!
!       COND  = DWV/(1._r8+XLV2*WS/(Tdum*Tdum))  ! Asai (1965, J. Japan)
        COND  = rfac*DWV*(Tdum*Tdum)/((Tdum*Tdum)+XLV2*WS)    ! Asai (1965, J. Japan)
!       COND  =      DWV*tsq/(tsq+XLV2*WS)    ! Asai (1965, J. Japan)
        COND  = MAX(COND, -WCdum)             ! Limit cloud water evaporation
        Tdum  = Tdum+XLV1*COND                ! Updated temperature
        WVdum = WVdum-COND                    ! Updated water vapor mixing ratio
        WCdum = WCdum+COND                    ! Updated cloud water mixing ratio
        CONDENSE = CONDENSE + COND            ! Total cloud water condensation
        ESW  = min(PP, FPVSL(Tdum))           ! Updated saturation vapor press w/r/t water
!       WS   = RHgrd*con_eps*ESW/(PP-ESW)         ! Updated saturation mixing ratio w/r/t water
        WS   = RHgrd*con_eps*ESW/(PP+con_epsm1*ESW)   ! Updated saturation mixing ratio w/r/t water
        DWV  = WVdum-WS                       ! Deficit grid-scale water vapor mixing ratio
        SSAT = DWV / WS                       ! Grid-scale supersaturation ratio
        rfac = 1.0_r8
      ENDDO

      END FUNCTION CONDENSE
!
!#######################################################################
!---------------- Calculate ice deposition at T<T_ICE ------------------
!#######################################################################
!
      real(kind=r8) FUNCTION DEPOSIT (PP, RHgrd, Tdum, WVdum)
!
      implicit none
!
!--- Also uses the Asai (1965) algorithm, but uses a different target
!      vapor pressure for the adjustment
!
      real(kind=r8) PP, RHgrd, Tdum, WVdum
      INTEGER, PARAMETER :: HIGH_PRES=r8
!     INTEGER, PARAMETER :: HIGH_PRES=Selected_real(kind=r8)_Kind(15)
      REAL (KIND=r8), PARAMETER :: RHLIMIT=0.001_r8 
      REAL (KIND=r8), PARAMETER :: RHLIMIT1=-RHLIMIT
      real(kind=r8), PARAMETER :: RCP=1.0_r8/con_cp
      real(kind=r8), PARAMETER :: RCPRV=RCP/con_rv
      real(kind=r8), PARAMETER :: XLS=con_hvap+con_hfus
      real(kind=r8), PARAMETER :: XLS1=XLS*RCP
      real(kind=r8), PARAMETER :: XLS2=XLS*XLS*RCPRV
      REAL (KIND=r8)  :: DEP
      REAL (KIND=r8)  :: SSAT
      real(kind=r8) esi, ws, dwv
!
!-----------------------------------------------------------------------
!
      ESI=min(PP, FPVSI(Tdum))                  ! Saturation vapor press w/r/t ice
!     WS=RHgrd*con_eps*ESI/(PP-ESI)                 ! Saturation mixing ratio
      WS=RHgrd*con_eps*ESI/(PP+con_epsm1*ESI)           ! Saturation mixing ratio
      DWV=WVdum-WS                              ! Deficit grid-scale water vapor mixing ratio
      SSAT=DWV/WS                               ! Supersaturation ratio
      DEPOSIT=0.0_r8
      DO WHILE (SSAT > RHLIMIT .OR. SSAT < RHLIMIT1)
   !
   !--- Note that XLVS2=LS*LV/(CP*RV)=LV*WS/(RV*T*T)*(LS/CP*DEP1),
   !     where WS is the saturation mixing ratio following Clausius-
   !     Clapeyron (see Asai,1965; Young,1993,p.405)
   !
        DEP=DWV/(1.0_r8+XLS2*WS/(Tdum*Tdum))        ! Asai (1965, J. Japan)
        Tdum=Tdum+XLS1*DEP                      ! Updated temperature
        WVdum=WVdum-DEP                         ! Updated ice mixing ratio
        DEPOSIT=DEPOSIT+DEP                     ! Total ice deposition
        ESI=min(PP, FPVSI(Tdum))                ! Updated saturation vapor press w/r/t ice
!       WS=RHgrd*con_eps*ESI/(PP-ESI)               ! Updated saturation mixing ratio w/r/t ice
        WS=RHgrd*con_eps*ESI/(PP+con_epsm1*ESI)         ! Updated saturation mixing ratio w/r/t ice
        DWV=WVdum-WS                            ! Deficit grid-scale water vapor mixing ratio
        SSAT=DWV/WS                             ! Grid-scale supersaturation ratio
      ENDDO
      END FUNCTION DEPOSIT
!
      END  SUBROUTINE GSMCOLUMN



  SUBROUTINE INIT_MICRO(DT,ibmax,kMax,jbMax,F_ICE_PHY  ,F_RAIN_PHY ,F_RIMEF_PHY)
    !
    IMPLICIT NONE
    !

    REAL (kind=r8), INTENT(IN ) :: dt
    INTEGER       , INTENT(IN ) :: ibMax
    INTEGER       , INTENT(IN ) :: kMax
    INTEGER       , INTENT(IN ) :: jbMax
    REAL (kind=r8), INTENT(OUT) :: F_ICE_PHY   (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RAIN_PHY  (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RIMEF_PHY (ibMax,kMax,jbMax)
    REAL (kind=r8), PARAMETER   :: dtlast=0.0_r8
    REAL (kind=r8), PARAMETER   :: fhour=0.0_r8

    REAL (kind=r8) :: dtp
    REAL (kind=r8) :: dtphys
    INTEGER        :: nsphys

    !
    dtphys=3600.0_r8
    NSPHYS=MAX(INT(2.0_r8*dt/DTPHYS+0.9999_r8),1)
    DTP=(dt+dt) / NSPHYS
    !
    !IF (num_p3d .EQ. 3 .AND. dtp .NE. dtlast) THEN

       !      Initialization and/or constant evaluation for Ferrier's microphysics

       CALL MICRO_INIT(ibMax,kMax,jbMax, F_ICE_PHY  ,&
                       F_RAIN_PHY ,F_RIMEF_PHY, DTP, FHOUR)
    !ENDIF

    RETURN
  END SUBROUTINE INIT_MICRO



  SUBROUTINE MICRO_INIT(ibMax,kMax,jbmax,F_ICE_PHY,F_RAIN_PHY , &
                        F_RIMEF_PHY , DT ,FHOUR)
    !
    !     This subroutine initializes the necessary constants and
    !     tables for Brad Ferrier's cloud microphysics package
    !
    IMPLICIT NONE
    !
    INTEGER       , INTENT(IN ) :: ibMax
    INTEGER       , INTENT(IN ) :: kMax
    INTEGER       , INTENT(IN ) :: jbMax
    REAL (kind=r8), INTENT(OUT) :: F_ICE_PHY   (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RAIN_PHY  (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(OUT) :: F_RIMEF_PHY (ibMax,kMax,jbMax)
    REAL (kind=r8), INTENT(IN ) :: DT
    REAL (kind=r8), INTENT(IN ) :: FHOUR
    !
    IF (fhour .LT. 0.1_r8) THEN
       F_ICE_PHY = 0.0_r8   ! Initialize ice  fraction array (real)
       F_RAIN_PHY= 0.0_r8   ! Initialize rain fraction array (real)
       F_RIMEF_PHY= 1.0_r8   ! Initialize rime factor   array (real)
    ENDIF

    CALL GSMCONST (DT) ! Initialize lookup tables & constants
    !
    RETURN
  END SUBROUTINE MICRO_INIT
  !
  !#######################################################################
  !------- Initialize constants & lookup tables for microphysics ---------
  !#######################################################################
  !
  SUBROUTINE GSMCONST (dt)
    !
    IMPLICIT NONE
    !-------------------------------------------------------------------------------
    !---  SUBPROGRAM DOCUMENTATION BLOCK
    !   PRGRMMR: Ferrier         ORG: W/NP22     DATE: February 2001
    !-------------------------------------------------------------------------------
    ! ABSTRACT:
    !   * Reads various microphysical lookup tables used in COLUMN_MICRO
    !   * Lookup tables were created "offline" and are read in during execution
    !   * Creates lookup tables for saturation vapor pressure w/r/t water & ice
    !-------------------------------------------------------------------------------
    !
    ! USAGE: CALL GSMCONST FROM SUBROUTINE GSMDRIVE AT MODEL START TIME
    !
    !   INPUT ARGUMENT LIST:
    !       DTPH - physics time step (s)
    !
    !   OUTPUT ARGUMENT LIST:
    !     NONE
    !
    !   OUTPUT FILES:
    !     NONE
    !
    !
    !   SUBROUTINES:
    !     MY_GROWTH_RATES - lookup table for growth of nucleated ice
    !
    !   UNIQUE: NONE
    !
    !   LIBRARY: NONE
    !
    !
    ! ATTRIBUTES:
    !   LANGUAGE: FORTRAN 90
    !   MACHINE : IBM SP
    !
    REAL(kind=r8)    , INTENT(IN   ) :: dt
    !
    !--- Parameters & data statement for local calculations
    !
    REAL(kind=r8), PARAMETER :: C1=1.0_r8/3.0_r8
    REAL(kind=r8), PARAMETER :: DMR1=0.1E-3_r8
    REAL(kind=r8), PARAMETER :: DMR2=0.2E-3_r8
    REAL(kind=r8), PARAMETER :: DMR3=0.32E-3_r8
    REAL(kind=r8), PARAMETER :: N0r0=8.E6_r8
    REAL(kind=r8), PARAMETER :: N0s0=4.E6_r8
    REAL(kind=r8), PARAMETER :: RHOL=1000.0_r8
    REAL(kind=r8), PARAMETER :: RHOS=100.0_r8
    REAL(kind=r8), PARAMETER :: XMR1=1.e6_r8*DMR1
    REAL(kind=r8), PARAMETER :: XMR2=1.e6_r8*DMR2
    REAL(kind=r8), PARAMETER :: XMR3=1.e6_r8*DMR3
    INTEGER      , PARAMETER :: MDR1=INT(XMR1)
    INTEGER      , PARAMETER :: MDR2=INT(XMR2)
    INTEGER      , PARAMETER :: MDR3=INT(XMR3)
    REAL(KIND=r8)   ,PARAMETER   :: GSMDT=0.0_r8
    !
    REAL(kind=r8) :: dtph
    REAL(kind=r8) :: bbfr
    INTEGER i
    !
    !--- Added on 5/16/01 for Moorthi
    !
    !      logical, parameter :: read_lookup=.false., write_lookup=.false.
    !
    !------------------------------------------------------------------------
    !  *************  Parameters used in ETA model -- Not used in Global Model *****
    !
    !--- DPHD, DLMD are delta latitude and longitude at the model (NOT geodetic) equator
    !    => "DX" is the hypotenuse of the model zonal & meridional grid increments.
    !
    !     DX=111.*(DPHD**2+DLMD**2)**.5         ! Resolution at MODEL equator (km)
    !     DX=MIN(100., MAX(5., DX) )
    !
    !--- Assume the following functional relationship for key constants that
    !    depend on grid resolution from DXmin (5 km) to DXmax (100 km) resolution:
    !
    !     DXmin=5.
    !     DXmax=100.
    !     DX=MIN(DXmax, MAX(DXmin, DX) )
    !
    !--- EXtune determines the degree to which the coefficients change with resolution.
    !    The larger EXtune is, the more sensitive the parameter.
    !
    !     EXtune=1.

    !
    !--- FXtune ==> F(DX) is the grid-resolution tuning parameter (from 0 to 1)
    !
    !     FXtune=((DXmax-DX)/(DXmax-DXmin))**EXtune
    !     FXtune=MAX(0., MIN(1., FXtune))
    !
    !--- Calculate grid-averaged RH for the onset of condensation (RHgrd) based on
    !    simple ***ASSUMED*** (user-specified) values at DXmax and at DXmin.
    !
    !     RH_DXmax=.90              !-- 90% RH at DXmax=100 km
    !     RH_DXmin=.98              !-- 98% RH at DXmin=5 km
    !
    !--- Note that RHgrd is right now fixed throughout the domain!!
    !
    !     RHgrd=RH_DXmax+(RH_DXmin-RH_DXmax)*FXtune
    !   ********************************************************************************
    !
    !
    !      if (first) then
    !
    !--- Read in various lookup tables
    !
    !      if ( read_lookup ) then
    !        OPEN (UNIT=1,FILE="eta_micro_lookup.dat",FORM="UNFORMATTED")
    !        READ(1) VENTR1
    !        READ(1) VENTR2
    !        READ(1) ACCRR
    !        READ(1) MASSR
    !        READ(1) VRAIN
    !        READ(1) RRATE
    !        READ(1) VENTI1
    !        READ(1) VENTI2
    !        READ(1) ACCRI
    !        READ(1) MASSI
    !        READ(1) VSNOWI
    !        READ(1) VEL_RF
    !!       read(1) my_growth    ! Applicable only for DTPH=180 s for offline testing
    !        CLOSE (1)
    !      else
    !        etime1=timef()
    CALL ICE_LOOKUP()                   ! Lookup tables for ice
    !        etime2=timef()
    !        if (mype == 0)                                                  &
    !     &  print *,'CPU time (sec) in ICE_LOOKUP = ',(etime2-etime1)*0.001
    CALL RAIN_LOOKUP()                  ! Lookup tables for rain
    !        etime1=timef()
    !        if (mype == 0)                                                  &
    !     &  print *,'CPU time (sec) in RAIN_LOOKUP = ',(etime1-etime2)*0.001
    !        if (write_lookup) then
    !          open(unit=1,file='micro_lookup.dat',form='unformatted')
    !          write(1) ventr1
    !          write(1) ventr2
    !          write(1) accrr
    !          write(1) massr
    !          write(1) vrain
    !          write(1) rrate
    !          write(1) venti1
    !          write(1) venti2
    !          write(1) accri
    !          write(1) massi
    !          write(1) vsnowi
    !          write(1) vel_rf
    !!         write(1) my_growth    ! Applicable only for DTPH=180 s ????
    !          CLOSE (1)
    !        endif
    !      endif
    !!
    !--- Constants associated with Biggs (1953) freezing of rain, as parameterized
    !    following Lin et al. (JCAM, 1983) & Reisner et al. (1998, QJRMS).
    !
    ABFR=-0.66_r8
    BBFR=100.0_r8
    CBFR=20.0_r8*con_pi*con_pi*BBFR*RHOL*1.E-21_r8
    !
    !--- QAUT0 is the threshold cloud content for autoconversion to rain
    !      needed for droplets to reach a diameter of 20 microns (following
    !      Manton and Cotton, 1977; Banta and Hanson, 1987, JCAM).  It is
    !      **STRONGLY** affected by the assumed droplet number concentrations
    !     XNCW!  For example, QAUT0=1.2567, 0.8378, or 0.4189 g/m**3 for
    !     droplet number concentrations of 300, 200, and 100 cm**-3, respectively.
    !
    !--- Calculate grid-averaged XNCW based on simple ***ASSUMED*** (user-specified)
    !    values at DXmax and at DXmin.
    !
    !     XNCW_DXmax=50.E6          !--  50 /cm**3 at DXmax=100 km
    !     XNCW_DXmin=200.E6         !-- 200 /cm**3 at DXmin=5 km
    !
    !--- Note that XNCW is right now fixed throughout the domain!!
    !
    !     XNCW=XNCW_DXmax+(XNCW_DXmin-XNCW_DXmax)*FXtune
    !
    !     QAUT0=con_pi*RHOL*XNCW*(20.E-6_r8)**3/6.0_r8
    QAUTx=con_pi*RHOL*1.0E6_r8*(20.E-6_r8)**3/6.0_r8
    !
    !--- Based on rain lookup tables for mean diameters from 0.05 to 0.45 mm
    !    * Four different functional relationships of mean drop diameter as
    !      a function of rain rate (RR), derived based on simple fits to
    !      mass-weighted fall speeds of rain as functions of mean diameter
    !      from the lookup tables.
    !
    RR_DRmin=N0r0*RRATE(MDRmin)     ! RR for mean drop diameter of .05 mm
    RR_DR1=N0r0*RRATE(MDR1)         ! RR for mean drop diameter of .10 mm
    RR_DR2=N0r0*RRATE(MDR2)         ! RR for mean drop diameter of .20 mm
    RR_DR3=N0r0*RRATE(MDR3)         ! RR for mean drop diameter of .32 mm
    RR_DRmax=N0r0*RRATE(MDRmax)     ! RR for mean drop diameter of .45 mm
    !
    RQR_DRmin=N0r0*MASSR(MDRmin)    ! Rain content for mean drop diameter of .05 mm
    RQR_DR1=N0r0*MASSR(MDR1)        ! Rain content for mean drop diameter of .10 mm
    RQR_DR2=N0r0*MASSR(MDR2)        ! Rain content for mean drop diameter of .20 mm
    RQR_DR3=N0r0*MASSR(MDR3)        ! Rain content for mean drop diameter of .32 mm
    RQR_DRmax=N0r0*MASSR(MDRmax)    ! Rain content for mean drop diameter of .45 mm
    C_N0r0=con_pi*RHOL*N0r0
    CN0r0=1.E6_r8/C_N0r0**0.25_r8
    CN0r_DMRmin=1.0_r8/(con_pi*RHOL*DMRmin**4)
    CN0r_DMRmax=1.0_r8/(con_pi*RHOL*DMRmax**4)
    !
    !      endif                     !  If (first) then loop ends here
    !
    !     Find out what microphysics time step should be
    !
    !mic_step  = MAX(1, INT(dt/600.0_r8+0.5_r8))
    mic_step = max(1, int(dt/100.0_r8+0.5_r8))
    dtph     = dt / mic_step
    !      if (mype == 0) print *,' dt=',dt,' mic_step=',mic_step        &
    !     &,                ' dtph=',dtph
    !
    !--- Calculates coefficients for growth rates of ice nucleated in water
    !    saturated conditions, scaled by physics time step (lookup table)
    !
    CALL MY_GROWTH_RATES (DTPH)
    !
    !--- CIACW is used in calculating riming rates
    !      The assumed effective collection efficiency of cloud water rimed onto
    !      ice is =0.5_r8 below:
    !
    !Moor CIACW=DTPH*0.25_r8*con_pi*0.5_r8*(1.E5_r8)**C1   ! commented on 20050422
    !      ice is =0.1_r8 below:
    CIACW=DTPH*0.25_r8*con_pi*0.1_r8*(1.E5_r8)**C1! original
    !     CIACW = 0.0_r8      ! Brad's suggestion 20040614
    !
    !--- CIACR is used in calculating freezing of rain colliding with large ice
    !      The assumed collection efficiency is 1.0_r8
    !
    CIACR=con_pi*DTPH
    !
    !--- CRACW is used in calculating collection of cloud water by rain (an
    !      assumed collection efficiency of 1.0_r8)
    !
    !Moor CRACW=DTPH*0.25_r8*con_pi*1.0_r8                 ! commented on 20050422
    !
    !      assumed collection efficiency of 0.1_r8)
    CRACW=DTPH*0.25_r8*con_pi*0.1_r8 ! original

    !     CRACW = 0.0_r8      ! Brad's suggestion 20040614
    !
    ESW0=FPVSL(T0C)           ! Saturation vapor pressure at 0C
    RFmax=1.1_r8**Nrime          ! Maximum rime factor allowed
    !
    !------------------------------------------------------------------------
    !--------------- Constants passed through argument list -----------------
    !------------------------------------------------------------------------
    !
    !--- Important parameters for self collection (autoconversion) of
    !    cloud water to rain.
    !
    !--- CRAUT is proportional to the rate that cloud water is converted by
    !      self collection to rain (autoconversion rate)
    !
    CRAUT=1.0_r8-EXP(-1.E-3_r8*DTPH)
    !
    !     IF (MYPE == 0)
    !    & WRITE(6,"(A, A,F6.2,A, A,F5.4, A,F7.3,A, A,F6.2,A, A,F5.3,A)")
    !    &   'KEY MICROPHYSICAL PARAMETERS FOR '
    !    &  ,'DX=',DX,' KM:'
    !    &  ,'   FXtune=',FXtune
    !    &  ,'   RHgrd=',100.*RHgrd,' %'
    !    &  ,'   NCW=',1.E-6*XNCW,' /cm**3'
    !    &  ,'   QAUT0=',1.E3*QAUT0,' g/kg'
    !
    !--- For calculating snow optical depths by considering bulk density of
    !      snow based on emails from Q. Fu (6/27-28/01), where optical
    !      depth (T) = 1.5*SWP/(Reff*DENS), SWP is snow water path, Reff
    !      is effective radius, and DENS is the bulk density of snow.
    !
    !    SWP (kg/m**2)=(1.E-3 kg/g)*SWPrad, SWPrad in g/m**2 used in radiation
    !    T = 1.5*1.E3*SWPrad/(Reff*DENS)
    !
    !    See derivation for MASSI(INDEXS), note equal to RHO*QSNOW/NSNOW
    !
    !    SDENS=1.5e3/DENS, DENS=MASSI(INDEXS)/[con_pi*(1.E-6*INDEXS)**3]
    !
    DO I=MDImin,MDImax
       !MoorthiSDENS(I)=con_pi*1.5E-15_r8*FLOAT(I*I*I)/MASSI(I)
            SDENS(I)=con_pi*1.0E-15_r8*FLOAT(I*I*I)/MASSI(I)
       !        SDENS(I)=con_pi*1.5E-15_r8*FLOAT(I*I*I)/MASSI(I)

    ENDDO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE gsmconst

  !
  !
  !#######################################################################
  !--- Sets up lookup table for calculating initial ice crystal growth ---
  !#######################################################################
  !
  SUBROUTINE MY_GROWTH_RATES (DTPH)
    !
    IMPLICIT NONE
    !
    !--- Below are tabulated values for the predicted mass of ice crystals
    !    after 600 s of growth in water saturated conditions, based on 
    !    calculations from Miller and Young (JAS, 1979).  These values are
    !    crudely estimated from tabulated curves at 600 s from Fig. 6.9 of
    !    Young (1993).  Values at temperatures colder than -27C were 
    !    assumed to be invariant with temperature.  
    !
    !--- Used to normalize Miller & Young (1979) calculations of ice growth
    !    over large time steps using their tabulated values at 600 s.
    !    Assumes 3D growth with time**1.5 following eq. (6.3) in Young (1993).
    !
    REAL(kind=r8) dtph, dt_ice
    REAL(kind=r8) MY_600(MY_T1:MY_T2)
    !
    !-- 20090714: These values are in g and need to be converted to kg below
    DATA MY_600 /                                             &
         5.5e-8_r8 , 1.4E-7_r8 , 2.8E-7_r8, 6.E-7_r8  , 3.3E-6_r8 ,& !  -1 to  -5 deg C
         2.E-6_r8  , 9.E-7_r8  , 8.8E-7_r8, 8.2E-7_r8 , 9.4e-7_r8 ,& !  -6 to -10 deg C
         1.2E-6_r8 , 1.85E-6_r8, 5.5E-6_r8, 1.5E-5_r8 , 1.7E-5_r8 ,& ! -11 to -15 deg C
         1.5E-5_r8 , 1.E-5_r8  , 3.4E-6_r8, 1.85E-6_r8, 1.35E-6_r8,& ! -16 to -20 deg C
         1.05E-6_r8, 1.E-6_r8  , 9.5E-7_r8, 9.0E-7_r8, 9.5E-7_r8 ,& ! -21 to -25 deg C
         9.5E-7_r8 , 9.E-7_r8  , 9.E-7_r8 , 9.E-7_r8  , 9.E-7_r8  ,& ! -26 to -30 deg C
         9.E-7_r8  , 9.E-7_r8  , 9.E-7_r8 , 9.E-7_r8  , 9.E-7_r8/    ! -31 to -35 deg C
    !
    !-----------------------------------------------------------------------
    !
    DT_ICE=(DTPH/600.0_r8)**1.5_r8
    !     MY_GROWTH=DT_ICE*MY_600             ! original version
    MY_GROWTH=DT_ICE*MY_600*1.E-3_r8    !-- 20090714: Convert from g to kg
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE MY_GROWTH_RATES
  !

  !
  !#######################################################################
  !--------------- Creates lookup tables for ice processes ---------------
  !#######################################################################
  !
  SUBROUTINE ice_lookup()
    !
    IMPLICIT NONE
    !-----------------------------------------------------------------------------------
    !
    !---- Key diameter values in mm
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Key concepts:
    !       - Actual physical diameter of particles (D)
    !       - Ratio of actual particle diameters to mean diameter (x=D/MD)
    !       - Mean diameter of exponentially distributed particles, which is the
    !         same as 1./LAMDA of the distribution (MD)
    !       - All quantitative relationships relating ice particle characteristics as
    !         functions of their diameter (e.g., ventilation coefficients, normalized
    !         accretion rates, ice content, and mass-weighted fall speeds) are a result
    !         of using composite relationships for ice crystals smaller than 1.5 mm
    !         diameter merged with analogous relationships for larger sized aggregates.
    !         Relationships are derived as functions of mean ice particle sizes assuming
    !         exponential size spectra and assuming the properties of ice crystals at
    !         sizes smaller than 1.5 mm and aggregates at larger sizes.  
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Actual ice particle diameters for which integrated distributions are derived
    !       - DminI - minimum diameter for integration (.02 mm, 20 microns)
    !       - DmaxI - maximum diameter for integration (2 cm)
    !       - DdelI - interval for integration (1 micron)
    !
    REAL(kind=r8)   , PARAMETER :: DminI=.02e-3_r8
    REAL(kind=r8)   , PARAMETER :: DmaxI=20.e-3_r8
    REAL(kind=r8)   , PARAMETER :: DdelI=1.e-6_r8     
    REAL(kind=r8)   , PARAMETER :: XImin=1.e6_r8*DminI
    REAL(kind=r8)   , PARAMETER :: XImax=1.e6_r8*DmaxI
    INTEGER         , PARAMETER :: IDImin=INT(XImin)
    INTEGER         , PARAMETER :: IDImax=INT(XImax)
    !
    !---- Meaning of the following arrays:
    !        - diam - ice particle diameter (m)
    !        - mass - ice particle mass (kg)
    !        - vel  - ice particle fall speeds (m/s)
    !        - vent1 - 1st term in ice particle ventilation factor
    !        - vent2 - 2nd term in ice particle ventilation factor
    !
    REAL(kind=r8)               :: diam (IDImin:IDImax)
    REAL(kind=r8)               :: mass (IDImin:IDImax)
    REAL(kind=r8)               :: vel  (IDImin:IDImax)
    REAL(kind=r8)               :: vent1(IDImin:IDImax)
    REAL(kind=r8)               :: vent2(IDImin:IDImax)
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Found from trial & error that the m(D) & V(D) mass & velocity relationships
    !       between the ice crystals and aggregates overlapped & merged near a particle
    !       diameter sizes of 1.5 mm.  Thus, ice crystal relationships are used for
    !       sizes smaller than 1.5 mm and aggregate relationships for larger sizes.
    !
    REAL(kind=r8), PARAMETER :: d_crystal_max=1.5_r8
    !
    !---- The quantity xmax represents the maximum value of "x" in which the
    !       integrated values are calculated.  For xmax=20., this means that
    !       integrated ventilation, accretion, mass, and precipitation rates are
    !       calculated for ice particle sizes less than 20.*mdiam, the mean particle diameter.
    !
    REAL(kind=r8), PARAMETER :: xmax=20.0_r8
    !
    !-----------------------------------------------------------------------------------
    !
    !---- Meaning of the following arrays:
    !        - mdiam - mean diameter (m)
    !        - VENTI1 - integrated quantity associated w/ ventilation effects
    !                   (capacitance only) for calculating vapor deposition onto ice
    !        - VENTI2 - integrated quantity associated w/ ventilation effects
    !                   (with fall speed) for calculating vapor deposition onto ice
    !        - ACCRI  - integrated quantity associated w/ cloud water collection by ice
    !        - MASSI  - integrated quantity associated w/ ice mass 
    !        - VSNOWI - mass-weighted fall speed of snow, used to calculate precip rates
    !
    !--- Mean ice-particle diameters varying from 50 microns to 1000 microns (1 mm), 
    !      assuming an exponential size distribution.  
    !
    REAL(kind=r8) ::    mdiam
    !
    !-----------------------------------------------------------------------------------
    !------------- Constants & parameters for ventilation factors of ice ---------------
    !-----------------------------------------------------------------------------------
    !
    !---- These parameters are used for calculating the ventilation factors for ice
    !       crystals between 0.2 and 1.5 mm diameter (Hall and Pruppacher, JAS, 1976).  
    !       From trial & error calculations, it was determined that the ventilation
    !       factors of smaller ice crystals could be approximated by a simple linear
    !       increase in the ventilation coefficient from 1.0 at 50 microns (.05 mm) to 
    !       1.1 at 200 microns (0.2 mm), rather than using the more complex function of
    !       1.0 + .14*(Sc**.33*Re**.5)**2 recommended by Hall & Pruppacher.
    !
    REAL(kind=r8), PARAMETER :: cvent1i=.86_r8
    REAL(kind=r8), PARAMETER :: cvent2i=.28_r8
    !
    !---- These parameters are used for calculating the ventilation factors for larger
    !       aggregates, where D>=1.5 mm (see Rutledge and Hobbs, JAS, 1983; 
    !       Thorpe and Mason, 1966).
    !
    REAL(kind=r8), PARAMETER :: cvent1a=.65_r8
    REAL(kind=r8), PARAMETER :: cvent2a=.44_r8
    !
    REAL(kind=r8)            :: m_agg
    REAL(kind=r8)            :: m_bullet
    REAL(kind=r8)            :: m_column
    REAL(kind=r8)            :: m_ice
    REAL(kind=r8)            :: m_plate
    !
    !---- Various constants
    !
    REAL(kind=r8), PARAMETER :: c1=2.0_r8/3.0_r8
    REAL(kind=r8), PARAMETER :: cexp=1.0_r8/3.0_r8
    !
    !LOGICAL            :: iprint
    !LOGICAL, PARAMETER :: print_diag=.FALSE.
    !
    !-----------------------------------------------------------------------------------
    !- Constants & parameters for calculating the increase in fall speed of rimed ice --
    !-----------------------------------------------------------------------------------
    !
    !---- Constants & arrays for estimating increasing fall speeds of rimed ice.
    !     Based largely on theory and results from Bohm (JAS, 1989, 2419-2427).
    !
    !-------------------- Standard atmosphere conditions at 1000 mb --------------------
    !
    REAL(kind=r8), PARAMETER :: t_std=288.0_r8
    REAL(kind=r8), PARAMETER :: dens_std=1000.e2_r8/(287.04_r8*288.0_r8)
    !
    !---- These "bulk densities" are the actual ice densities in the ice portion of the 
    !     lattice.  They are based on text associated w/ (12) on p. 2425 of Bohm (JAS, 
    !     1989).  Columns, plates, & bullets are assumed to have an average bulk density 
    !     of 850 kg/m**3.  Aggregates were assumed to have a slightly larger bulk density 
    !     of 600 kg/m**3 compared with dendrites (i.e., the least dense, most "lacy" & 
    !     tenous ice crystal, which was assumed to be ~500 kg/m**3 in Bohm).  
    !
    REAL(kind=r8), PARAMETER :: dens_crystal=850.0_r8
    REAL(kind=r8), PARAMETER :: dens_agg=600.0_r8
    !
    !--- A value of Nrime=40 for a logarithmic ratio of 1.1 yields a maximum rime factor
    !      of 1.1**40 = 45.26 that is resolved in these tables.  This allows the largest
    !      ice particles with a mean diameter of MDImax=1000 microns to achieve bulk 
    !      densities of 900 kg/m**3 for rimed ice.  
    !
    !     integer, parameter :: Nrime=40
    REAL(kind=r8)  :: m_rime
    REAL(kind=r8)  :: rime_factor(0:Nrime)
    REAL(kind=r8)  :: rime_vel   (0:Nrime)
    REAL(kind=r8)  :: vel_rime   (IDImin:IDImax,Nrime)
    REAL(kind=r8)  :: ivel_rime  (MDImin:MDImax,Nrime)
    !
    INTEGER :: i, j, jj, k, icount
    REAL(kind=r8)  :: c2
    REAL(kind=r8)  :: cbulk
    REAL(kind=r8)  :: cbulk_ice
    REAL(kind=r8)  :: px
    REAL(kind=r8)  :: dynvis_std
    REAL(kind=r8)  :: crime1
    REAL(kind=r8)  :: crime2
    REAL(kind=r8)  :: crime3
    REAL(kind=r8)  :: crime4
    REAL(kind=r8)  :: crime5
    REAL(kind=r8)  :: d
    REAL(kind=r8)  :: c_avg
    REAL(kind=r8)  :: c_agg
    REAL(kind=r8)  :: c_bullet
    REAL(kind=r8)  :: c_column
    REAL(kind=r8)  :: c_plate
    REAL(kind=r8)  :: cl_agg
    REAL(kind=r8)  :: cl_bullet
    REAL(kind=r8)  :: cl_column
    REAL(kind=r8)  :: cl_plate
    REAL(kind=r8)  :: v_agg
    REAL(kind=r8)  :: v_bullet
    REAL(kind=r8)  :: v_column
    REAL(kind=r8)  :: v_plate
    REAL(kind=r8)  :: wd
    REAL(kind=r8)  :: ecc_column
    REAL(kind=r8)  :: cvent1
    REAL(kind=r8)  :: cvent2
    REAL(kind=r8)  :: crime_best
    REAL(kind=r8)  :: rime_m1
    REAL(kind=r8)  :: rime_m2
    REAL(kind=r8)  :: x_rime
    REAL(kind=r8)  :: re_rime
    REAL(kind=r8)  :: smom3
    REAL(kind=r8)  :: pratei
    REAL(kind=r8)  :: expf
    REAL(kind=r8)  :: bulk_dens
    REAL(kind=r8)  :: xmass
    !REAL(kind=r8)  :: xmdiam
    REAL(kind=r8)  :: ecc_plate
    REAL(kind=r8)  :: dx
    !
    !-----------------------------------------------------------------------------------
    !----------------------------- BEGIN EXECUTION -------------------------------------
    !-----------------------------------------------------------------------------------
    !
    !
    c2=1.0_r8/SQRT(3.0_r8)
    !     pi=acos(-1._r8)
    cbulk=6.0_r8/con_pi
    cbulk_ice=900.0_r8*con_pi/6.0_r8    ! Maximum bulk ice density allowed of 900 kg/m**3
    px=0.4_r8**cexp             ! Convert fall speeds from 400 mb (Starr & Cox) to 1000 mb
    !
    !--------------------- Dynamic viscosity (1000 mb, 288 K) --------------------------
    !
    dynvis_std=1.496e-6_r8*t_std**1.5_r8/(t_std+120.0_r8)
    crime1=con_pi/24.0_r8
    crime2=8.0_r8*9.81_r8*dens_std/(con_pi*dynvis_std**2)
    crime3=crime1*dens_crystal
    crime4=crime1*dens_agg
    crime5=dynvis_std/dens_std
    DO i=0,Nrime
       rime_factor(i)=1.1_r8**i
    ENDDO
    !
    !#######################################################################
    !      Characteristics as functions of actual ice particle diameter 
    !#######################################################################
    !
    !----   M(D) & V(D) for 3 categories of ice crystals described by Starr 
    !----   & Cox (1985). 
    !
    !----   Capacitance & characteristic lengths for Reynolds Number calculations
    !----   are based on Young (1993; p. 144 & p. 150).  c-axis & a-axis 
    !----   relationships are from Heymsfield (JAS, 1972; Table 1, p. 1351).
    !
    icount=60
    !
    !      if (print_diag)                                                   & 
    !     &  write(7,"(2a)") '---- Increase in fall speeds of rimed ice',    &
    !     &    ' particles as function of ice particle diameter ----'
    DO i=IDImin,IDImax
       !        if (icount == 60 .and. print_diag) then
       !          write(6,"(/2a/3a)") 'Particle masses (mg), fall speeds ',     &
       !     &      '(m/s), and ventilation factors',                           &
       !     &      '  D(mm)  CR_mass   Mass_bull   Mass_col  Mass_plat ',      &
       !     &      '  Mass_agg   CR_vel  V_bul CR_col CR_pla Aggreg',          &
       !     &      '    Vent1      Vent2 '                               
       !          write(7,"(3a)") '        <----------------------------------',&
       !     &      '---------------  Rime Factor  --------------------------', &
       !     &      '--------------------------->'
       !          write(7,"(a,23f5.2)") '  D(mm)',(rime_factor(k), k=1,5),      &
       !     &       (rime_factor(k), k=6,40,2)
       !          icount=0
       !        endif
       d=(float(i)+0.5_r8)*1.e-3_r8    ! in mm
       c_avg=0.0_r8
       c_agg=0.0_r8
       c_bullet=0.0_r8
       c_column=0.0_r8
       c_plate=0.0_r8
       cl_agg=0.0_r8
       cl_bullet=0.0_r8
       cl_column=0.0_r8
       cl_plate=0._r8
       m_agg=0.0_r8
       m_bullet=0.0_r8
       m_column=0.0_r8
       m_plate=0.0_r8
       v_agg=0.0_r8
       v_bullet=0.0_r8
       v_column=0.0_r8
       v_plate=0.0_r8
       IF (d < d_crystal_max) THEN
          !
          !---- This block of code calculates bulk characteristics based on average
          !     characteristics of bullets, plates, & column ice crystals <1.5 mm size
          !
          !---- Mass-diameter relationships from Heymsfield (1972) & used
          !       in Starr & Cox (1985), units in mg
          !---- "d" is maximum dimension size of crystal in mm, 
          !
          ! Mass of pure ice for spherical particles, used as an upper limit for the
          !   mass of small columns (<~ 80 microns) & plates (<~ 35 microns)
          !
          m_ice=0.48_r8*d**3   ! Mass of pure ice for spherical particle
          !
          m_bullet=MIN(0.044_r8*d**3, m_ice)
          m_column=MIN(0.017_r8*d**1.7_r8, m_ice)
          m_plate=MIN(0.026_r8*d**2.5_r8, m_ice)
          !
          mass(i)=m_bullet+m_column+m_plate
          !
          !---- These relationships are from Starr & Cox (1985), applicable at 400 mb
          !---- "d" is maximum dimension size of crystal in mm, dx in microns
          !
          dx=1000.0_r8*d            ! Convert from mm to microns
          IF (dx <= 200.0_r8) THEN
             v_column=8.114e-5_r8*dx**1.585_r8
             v_bullet=5.666e-5_r8*dx**1.663_r8
             v_plate=1.e-3_r8*dx
          ELSE IF (dx <= 400.0_r8) THEN
             v_column=4.995e-3_r8*dx**.807_r8
             v_bullet=3.197e-3_r8*dx**.902_r8
             v_plate=1.48e-3_r8*dx**.926_r8
          ELSE IF (dx <= 600.0_r8) THEN
             v_column=2.223e-2_r8*dx**.558_r8
             v_bullet=2.977e-2_r8*dx**.529_r8
             v_plate=9.5e-4_r8*dx
          ELSE IF (dx <= 800.0_r8) THEN
             v_column=4.352e-2_r8*dx**.453_r8
             v_bullet=2.144e-2_r8*dx**.581_r8
             v_plate=3.161e-3_r8*dx**.812_r8
          ELSE 
             v_column=3.833e-2_r8*dx**.472_r8
             v_bullet=3.948e-2_r8*dx**.489_r8
             v_plate=7.109e-3_r8*dx**.691_r8
          ENDIF
          !
          !---- Reduce fall speeds from 400 mb to 1000 mb
          !
          v_column=px*v_column
          v_bullet=px*v_bullet
          v_plate=px*v_plate
          !
          !---- DIFFERENT VERSION!  CALCULATES MASS-WEIGHTED CRYSTAL FALL SPEEDS
          !
          vel(i)=(m_bullet*v_bullet+m_column*v_column+m_plate*v_plate)/mass(i)
          mass(i)=mass(i)/3.0_r8
          !
          !---- Shape factor and characteristic length of various ice habits,
          !     capacitance is equal to 4*con_pi*(Shape factor)
          !       See Young (1993, pp. 143-152 for guidance)
          !
          !---- Bullets:
          !
          !---- Shape factor for bullets (Heymsfield, 1975)
          c_bullet=0.5_r8*d
          !---- Length-width functions for bullets from Heymsfield (JAS, 1972)
          IF (d > 0.3_r8) THEN
             wd=0.25_r8*d**0.7856_r8     ! Width (mm); a-axis
          ELSE
             wd=0.185_r8*d**0.552_r8
          ENDIF
          !---- Characteristic length for bullets (see first multiplicative term on right
          !       side of eq. 7 multiplied by crystal width on p. 821 of Heymsfield, 1975)
          cl_bullet=0.5_r8*con_pi*wd*(0.25_r8*wd+d)/(d+wd)
          !
          !---- Plates:
          !
          !---- Length-width function for plates from Heymsfield (JAS, 1972)
          wd=0.0449_r8*d**0.449_r8      ! Width or thickness (mm); c-axis
          !---- Eccentricity & shape factor for thick plates following Young (1993, p. 144)
          ecc_plate=SQRT(1.0_r8-wd*wd/(d*d))         ! Eccentricity
          c_plate=d*ecc_plate/ASIN(ecc_plate)    ! Shape factor
          !---- Characteristic length for plates following Young (1993, p. 150, eq. 6.6)
          cl_plate=d+2.0_r8*wd      ! Characteristic lengths for plates
          !
          !---- Columns:
          !
          !---- Length-width function for columns from Heymsfield (JAS, 1972)
          IF (d > 0.2_r8) THEN
             wd=0.1973_r8*d**0.414_r8    ! Width (mm); a-axis
          ELSE
             wd=0.5_r8*d             ! Width (mm); a-axis
          ENDIF
          !---- Eccentricity & shape factor for columns following Young (1993, p. 144)
          ecc_column=SQRT(1.0_r8-wd*wd/(d*d))                     ! Eccentricity
          c_column=ecc_column*d/LOG((1.0_r8+ecc_column)*d/wd)     ! Shape factor
          !---- Characteristic length for columns following Young (1993, p. 150, eq. 6.7)
          cl_column=(wd+2.0_r8*d)/(c1+c2*d/wd)       ! Characteristic lengths for columns
          !
          !---- Convert shape factor & characteristic lengths from mm to m for 
          !       ventilation calculations
          !
          c_bullet=0.001_r8*c_bullet
          c_plate=0.001_r8*c_plate
          c_column=0.001_r8*c_column
          cl_bullet=0.001_r8*cl_bullet
          cl_plate=0.001_r8*cl_plate
          cl_column=0.001_r8*cl_column
          !
          !---- Make a smooth transition between a ventilation coefficient of 1.0 at 50 microns
          !       to 1.1 at 200 microns
          !
          IF (d > 0.2_r8) THEN
             cvent1=cvent1i
             cvent2=cvent2i/3.0_r8
          ELSE
             cvent1=1.0_r8+0.1_r8*MAX(0.0_r8, d-0.05_r8)/0.15_r8
             cvent2=0.0_r8
          ENDIF
          !
          !---- Ventilation factors for ice crystals:
          !
          vent1(i)=cvent1*(c_bullet+c_plate+c_column)/3.0_r8
          vent2(i)=cvent2*(c_bullet*SQRT(v_bullet*cl_bullet)            &
               &                    +c_plate*SQRT(v_plate*cl_plate)               &
               &                    +c_column*SQRT(v_column*cl_column) )
          crime_best=crime3     ! For calculating Best No. of rimed ice crystals
       ELSE
          !
          !---- This block of code calculates bulk characteristics based on average
          !     characteristics of unrimed aggregates >= 1.5 mm using Locatelli & 
          !     Hobbs (JGR, 1974, 2185-2197) data.
          !
          !----- This category is a composite of aggregates of unrimed radiating 
          !-----   assemblages of dendrites or dendrites; aggregates of unrimed
          !-----   radiating assemblages of plates, side planes, bullets, & columns;
          !-----   aggregates of unrimed side planes (mass in mg, velocity in m/s)
          !
          m_agg=(0.073_r8*d**1.4_r8+0.037_r8*d**1.9_r8+0.04_r8*d**1.4_r8)/3.0_r8
          v_agg=(0.8_r8*d**0.16_r8+0.69_r8*d**0.41_r8+0.82_r8*d**0.12_r8)/3.0_r8
          mass(i)=m_agg
          vel(i)=v_agg
          !
          !---- Assume spherical aggregates
          !
          !---- Shape factor is the same as for bullets, = D/2
          c_agg=0.001_r8*0.5_r8*d         ! Units of m
          !---- Characteristic length is surface area divided by perimeter
          !       (.25*con_pi*D**2)/(con_pi*D**2) = D/4
          cl_agg=0.5_r8*c_agg         ! Units of m
          !
          !---- Ventilation factors for aggregates:
          !
          vent1(i)=cvent1a*c_agg
          vent2(i)=cvent2a*c_agg*SQRT(v_agg*cl_agg)
          crime_best=crime4     ! For calculating Best No. of rimed aggregates
       ENDIF
       !
       !---- Convert from shape factor to capacitance for ventilation factors
       !
       vent1(i)=4.0_r8*con_pi*vent1(i)
       vent2(i)=4.0_r8*con_pi*vent2(i)
       diam(i)=1.e-3_r8*d             ! Convert from mm to m
       mass(i)=1.e-6_r8*mass(i)       ! Convert from mg to kg
       !
       !---- Calculate increase in fall speeds of individual rimed ice particles
       !
       DO k=0,Nrime
          !---- Mass of rimed ice particle associated with rime_factor(k)
          rime_m1=rime_factor(k)*mass(i)
          rime_m2=cbulk_ice*diam(i)**3
          m_rime=MIN(rime_m1, rime_m2)
          !---- Best Number (X) of rimed ice particle combining eqs. (8) & (12) in Bohm
          x_rime=crime2*m_rime*(crime_best/m_rime)**0.25_r8
          !---- Reynolds Number for rimed ice particle using eq. (11) in Bohm
          re_rime=8.5_r8*(SQRT(1.0+0.1519_r8*SQRT(x_rime))-1.0_r8)**2
          rime_vel(k)=crime5*re_rime/diam(i)
       ENDDO
       DO k=1,Nrime
          vel_rime(i,k)=rime_vel(k)/rime_vel(0)
       ENDDO
       !IF (print_diag) THEN
       !   !
       !   !---- Determine if statistics should be printed out.
       !   !
       !   iprint=.FALSE.
       !   IF (d <= 1.) THEN
       !      IF (MOD(i,10) == 0) iprint=.TRUE.
       !   ELSE
       !      IF (MOD(i,100) == 0) iprint=.TRUE.
       !   ENDIF
       !   !IF (iprint) THEN
       !   !   WRITE(6,"(f7.4,5e11.4,1x,5f7.4,1x,2e11.4)")                 & 
       !   !        &        d,1.e6*mass(i),m_bullet,m_column,m_plate,m_agg,           &
       !   !        &        vel(i),v_bullet,v_column,v_plate,v_agg,                   &
       !   !        &        vent1(i),vent2(i)
       !   !   WRITE(7,"(f7.4,23f5.2)") d,(vel_rime(i,k), k=1,5),          &
       !   !        &        (vel_rime(i,k), k=6,40,2)
       !   !   icount=icount+1
       !   !ENDIF
       !ENDIF
    ENDDO
    !
    !#######################################################################
    !      Characteristics as functions of mean particle diameter
    !#######################################################################
    !
    VENTI1=0.0_r8
    VENTI2=0.0_r8
    ACCRI=0.0_r8
    MASSI=0.0_r8
    VSNOWI=0.0_r8
    VEL_RF=0.0_r8
    ivel_rime=0.0_r8
    icount=0
    !      if (print_diag) then
    !        icount=60
    !        write(6,"(/2a)") '------------- Statistics as functions of ',   &
    !     &    'mean particle diameter -------------'
    !        write(7,"(/2a)") '------ Increase in fall speeds of rimed ice', &
    !     &    ' particles as functions of mean particle diameter -----'
    !      endif
    DO j=MDImin,MDImax
       !        if (icount == 60 .AND. print_diag) then
       !          write(6,"(/2a)") 'D(mm)    Vent1      Vent2    ',             &
       !     &       'Accrete       Mass     Vel  Dens  '
       !          write(7,"(3a)") '      <----------------------------------',  &
       !     &      '---------------  Rime Factor  --------------------------', &
       !     &      '--------------------------->'
       !          write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),        &
       !     &       (rime_factor(k), k=6,40,2)
       !          icount=0
       !        endif
       mdiam=DelDMI*float(j)       ! in m
       smom3=0.0_r8
       pratei=0.0_r8
       rime_vel=0._r8                 ! Note that this array is being reused!
       DO i=IDImin,IDImax
          dx=diam(i)/mdiam
          IF (dx <= xmax) THEN      ! To prevent arithmetic underflows
             expf=EXP(-dx)*DdelI
             VENTI1(J)=VENTI1(J)+vent1(i)*expf
             VENTI2(J)=VENTI2(J)+vent2(i)*expf
             ACCRI(J)=ACCRI(J)+diam(i)*diam(i)*vel(i)*expf
             xmass=mass(i)*expf
             DO k=1,Nrime
                rime_vel(k)=rime_vel(k)+xmass*vel_rime(i,k)
             ENDDO
             MASSI(J)=MASSI(J)+xmass
             pratei=pratei+xmass*vel(i)
             smom3=smom3+diam(i)**3*expf
          ELSE
             EXIT
          ENDIF
       ENDDO
       !
       !--- Increased fall velocities functions of mean diameter (j),
       !      normalized by ice content, and rime factor (k) 
       !
       DO k=1,Nrime
          ivel_rime(j,k)=rime_vel(k)/MASSI(J)
       ENDDO
       !
       !--- Increased fall velocities functions of ice content at 0.1 mm
       !      intervals (j_100) and rime factor (k); accumulations here
       !
       jj=j/100
       IF (jj >= 2 .AND. jj <= 9) THEN
          DO k=1,Nrime
             VEL_RF(jj,k)=VEL_RF(jj,k)+ivel_rime(j,k)
          ENDDO
       ENDIF
       bulk_dens=cbulk*MASSI(J)/smom3
       VENTI1(J)=VENTI1(J)/mdiam
       VENTI2(J)=VENTI2(J)/mdiam
       ACCRI(J)=ACCRI(J)/mdiam
       VSNOWI(J)=pratei/MASSI(J)
       MASSI(J)=MASSI(J)/mdiam
       !        if (mod(j,10) == 0 .AND. print_diag) then
       !          xmdiam=1.e3*mdiam
       !          write(6,"(f5.3,4e11.4,f6.3,f8.3)") xmdiam,VENTI1(j),VENTI2(j),&
       !     &      ACCRI(j),MASSI(j),VSNOWI(j),bulk_dens
       !          write(7,"(f5.3,23f5.2)") xmdiam,(ivel_rime(j,k), k=1,5),      &
       !     &       (ivel_rime(j,k), k=6,40,2)
       !          icount=icount+1
       !        endif
    ENDDO
    !
    !--- Average increase in fall velocities rimed ice as functions of mean
    !      particle diameter (j, only need 0.1 mm intervals) and rime factor (k)
    !
    !      if (print_diag) then
    !        write(7,"(/2a)") ' ------- Increase in fall speeds of rimed ',  &
    !     &    'ice particles at reduced, 0.1-mm intervals  --------'
    !        write(7,"(3a)") '        <----------------------------------',  &
    !     &    '---------------  Rime Factor  --------------------------',   &
    !     &    '--------------------------->'
    !        write(7,"(a,23f5.2)") 'D(mm)',(rime_factor(k), k=1,5),          &
    !     &    (rime_factor(k), k=6,40,2)
    !      endif
    DO j=2,9
       VEL_RF(j,0)=1.0_r8
       DO k=1,Nrime
          VEL_RF(j,k)=0.01_r8*VEL_RF(j,k)
       ENDDO
       !        if (print_diag) write(7,"(f3.1,2x,23f5.2)") 0.1*j,              &
       !     &    (VEL_RF(j,k), k=1,5),(VEL_RF(j,k), k=6,40,2)
    ENDDO
    !
    !-----------------------------------------------------------------------------------
    !
  END SUBROUTINE ice_lookup

  !
  !#######################################################################
  !-------------- Creates lookup tables for rain processes ---------------
  !#######################################################################
  !
  SUBROUTINE rain_lookup()
    IMPLICIT NONE
    !
    !--- Parameters & arrays for fall speeds of rain as a function of rain drop
    !      diameter.  These quantities are integrated over exponential size
    !      distributions of rain drops at 1 micron intervals (DdelR) from minimum 
    !      drop sizes of .05 mm (50 microns, DminR) to maximum drop sizes of 10 mm 
    !      (DmaxR). 
    !
    REAL(kind=r8), PARAMETER :: DminR=0.05e-3_r8
    REAL(kind=r8), PARAMETER :: DmaxR=10.e-3_r8
    REAL(kind=r8), PARAMETER :: DdelR=1.e-6_r8
    REAL(kind=r8), PARAMETER :: XRmin=1.e6_r8*DminR
    REAL(kind=r8), PARAMETER :: XRmax=1.e6_r8*DmaxR
    INTEGER      , PARAMETER :: IDRmin=INT(XRmin)
    INTEGER      , PARAMETER :: IDRmax=INT(XRmax)
    REAL(kind=r8)            :: diam(IDRmin:IDRmax)
    REAL(kind=r8)            :: vel (IDRmin:IDRmax)
    !
    !--- Parameters rain lookup tables, which establish the range of mean drop
    !      diameters; from a minimum mean diameter of 0.05 mm (DMRmin) to a 
    !      maximum mean diameter of 0.45 mm (DMRmax).  The tables store solutions
    !      at 1 micron intervals (DelDMR) of mean drop diameter.  
    !
    REAL(kind=r8) :: mdiam
    !REAL(kind=r8) :: mass
    !
    !LOGICAL, PARAMETER :: print_diag=.FALSE.
    !
    REAL(kind=r8) :: d, cmass, pi2, expf
    INTEGER       :: i, j, i1, i2
    !
    !-----------------------------------------------------------------------
    !------- Fall speeds of rain as function of rain drop diameter ---------
    !-----------------------------------------------------------------------
    !
    DO i=IDRmin,IDRmax
       diam(i)=float(i)*DdelR
       d=100.0_r8*diam(i)         ! Diameter in cm
       IF (d <= 0.42_r8) THEN
          !
          !--- Rutledge & Hobbs (1983); vel (m/s), d (cm)
          !
          vel(i)=MAX(0.0_r8, -0.267_r8+51.5_r8*d-102.25_r8*d*d+75.7_r8*d**3)
       ELSE IF (d > 0.42_r8 .AND. d <= 0.58_r8) THEN
          !
          !--- Linear interpolation of Gunn & Kinzer (1949) data
          !
          vel(i)=8.92_r8+0.25_r8/(0.58_r8-0.42_r8)*(d-0.42_r8)
       ELSE
          vel(i)=9.17_r8
       ENDIF
    ENDDO
    DO i=1,100
       i1=(i-1)*100+IDRmin
       i2=i1+90
       !
       !--- Print out rain fall speeds only for D<=5.8 mm (.58 cm)
       !
       IF (diam(i1) > 0.58e-2_r8) EXIT
       !        if (print_diag) then
       !          write(6,"(1x)")
       !          write(6,"('D(mm)->  ',10f7.3)") (1000.*diam(j), j=i1,i2,10)
       !          write(6,"('V(m/s)-> ',10f7.3)") (vel(j), j=i1,i2,10)
       !        endif
    ENDDO
    !
    !-----------------------------------------------------------------------
    !------------------- Lookup tables for rain processes ------------------
    !-----------------------------------------------------------------------
    !
    !     pi=acos(-1.)
    pi2=2.0_r8*con_pi
    cmass=1000.0_r8*con_pi/6.0_r8
    !      if (print_diag) then
    !        write(6,"(/'Diam - Mean diameter (mm)'                          &
    !     &          /'VENTR1 - 1st ventilation coefficient (m**2)'          &
    !     &          /'VENTR2 - 2nd ventilation coefficient (m**3/s**.5)'    &
    !     &          /'ACCRR - accretion moment (m**4/s)'                    &
    !     &          /'RHO*QR - mass content (kg/m**3) for N0r=8e6'          &
    !     &          /'RRATE - rain rate moment (m**5/s)'                    &
    !     &          /'VR - mass-weighted rain fall speed (m/s)'             &
    !     &    /' Diam      VENTR1      VENTR2       ACCRR      ',           &
    !     &    'RHO*QR       RRATE    VR')")
    !      endif
    DO j=MDRmin,MDRmax
       mdiam=float(j)*DelDMR
       VENTR2(J)=0.0_r8
       ACCRR(J)=0.0_r8
       MASSR(J)=0.0_r8
       RRATE(J)=0.0_r8
       DO i=IDRmin,IDRmax
          expf=EXP(-diam(i)/mdiam)*DdelR
          VENTR2(J)=VENTR2(J)+diam(i)**1.5_r8*vel(i)**0.5_r8*expf
          ACCRR(J)=ACCRR(J)+diam(i)*diam(i)*vel(i)*expf
          MASSR(J)=MASSR(J)+diam(i)**3*expf
          RRATE(J)=RRATE(J)+diam(i)**3*vel(i)*expf
       ENDDO
       !
       !--- Derived based on ventilation, F(D)=0.78+.31*Schmidt**(1/3)*Reynold**.5,
       !      where Reynold=(V*D*rho/dyn_vis), V is velocity, D is particle diameter,
       !      rho is air density, & dyn_vis is dynamic viscosity.  Only terms 
       !      containing velocity & diameter are retained in these tables.  
       !
       VENTR1(J)=0.78_r8*pi2*mdiam**2
       VENTR2(J)=0.31_r8*pi2*VENTR2(J)
       !
       MASSR(J)=cmass*MASSR(J)
       RRATE(J)=cmass*RRATE(J)
       VRAIN(J)=RRATE(J)/MASSR(J)
       !       if (print_diag) write(6,"(f5.3,5g12.5,f6.3)") 1000.*mdiam,      &
       !     &    ventr1(j),ventr2(j),accrr(j),8.e6*massr(j),rrate(j),vrain(j)
    ENDDO
    !
    !-----------------------------------------------------------------------
    !
  END SUBROUTINE rain_lookup


  !-------------------------------------------------------------------------------
  SUBROUTINE gpvsl
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvsl        Compute saturation vapor pressure table over liquid
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvsl.
    !   Exact saturation vapor pressures are calculated in subprogram fpvslx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gpvsl
    !
    ! Subprograms called:
    !   (fpvslx)   inlinable function to compute saturation vapor pressure
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER  :: jx
    REAL(r8) :: xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsl-1)
    c2xpvsl=1.0_r8/xinc
    c1xpvsl=1.0_r8-xmin*c2xpvsl
    DO jx=1,nxpvsl
       x=xmin+(jx-1)*xinc
       t=x
       tbpvsl(jx)=fpvslx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvsl
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvslx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvslx       Compute saturation vapor pressure over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and liquid, and neglects the volume of the liquid.
    !   The model does account for the variation of the latent heat
    !   of condensation with temperature.  The ice option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   pvsl=fpvslx(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvslx     Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) :: fpvslx
    REAL(r8),INTENT(in):: t
    REAL(r8),PARAMETER:: dldt=con_cvap-con_cliq
    REAL(r8),PARAMETER:: heat=con_hvap
    REAL(r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(r8) tr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvslx=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvslx

  !
  !   fpvsl           Elementally compute saturation vapor pressure over liquid
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !-------------------------------------------------------------------------------
  ELEMENTAL FUNCTION fpvsl(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsl        Compute saturation vapor pressure over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvsl. See documentation for fpvslx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvsl is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   pvsl=fpvsl(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsl      Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)           :: fpvsl
    REAL(r8),INTENT(in):: t
    INTEGER  :: jx
    REAL(r8) :: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsl+c2xpvsl*t,1.0_r8),REAL(nxpvsl,r8))
    jx=INT(MIN(xj,nxpvsl-1.0_r8))
    fpvsl=tbpvsl(jx)+(xj-jx)*(tbpvsl(jx+1)-tbpvsl(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsl
!-------------------------------------------------------------------------------
  elemental function fpvsi(t)
!$$$     Subprogram Documentation Block
!
! Subprogram: fpvsi        Compute saturation vapor pressure over ice
!   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
!
! Abstract: Compute saturation vapor pressure from the temperature.
!   A linear interpolation is done between values in a lookup table
!   computed in gpvsi. See documentation for fpvsix for details.
!   Input values outside table range are reset to table extrema.
!   The interpolation accuracy is almost 6 decimal places.
!   On the Cray, fpvsi is about 4 times faster than exact calculation.
!   This function should be expanded inline in the calling routine.
!
! Program History Log:
!   91-05-07  Iredell             made into inlinable function
!   94-12-30  Iredell             expand table
! 1999-03-01  Iredell             f90 module
! 2001-02-26  Iredell             ice phase
!
! Usage:   pvsi=fpvsi(t)
!
!   Input argument list:
!     t          Real(r8) temperature in Kelvin
!
!   Output argument list:
!     fpvsi      Real(r8) saturation vapor pressure in Pascals
!
! Attributes:
!   Language: Fortran 90.
!
!$$$
    implicit none
    real(r8) ::  fpvsi
    real(r8),intent(in):: t
    integer ::jx
    real(r8):: xj
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=min(max(c1xpvsi+c2xpvsi*t,1.0_r8),real(nxpvsi,r8))
    jx=INT(min(xj,nxpvsi-1.0_r8))
    fpvsi=tbpvsi(jx)+(xj-jx)*(tbpvsi(jx+1)-tbpvsi(jx))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function fpvsi


!-------------------------------------------------------------------------------
  subroutine gpvsi()
!$$$     Subprogram Documentation Block
!
! Subprogram: gpvsi        Compute saturation vapor pressure table over ice
!   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
!
! Abstract: Computes saturation vapor pressure table as a function of
!   temperature for the table lookup function fpvsi.
!   Exact saturation vapor pressures are calculated in subprogram fpvsix.
!   The current implementation computes a table with a length
!   of 7501 for temperatures ranging from 180. to 330. Kelvin.
!
! Program History Log:
!   91-05-07  Iredell
!   94-12-30  Iredell             expand table
! 1999-03-01  Iredell             f90 module
! 2001-02-26  Iredell             ice phase
!
! Usage:  call gpvsi
!
! Subprograms called:
!   (fpvsix)   inlinable function to compute saturation vapor pressure
!
! Attributes:
!   Language: Fortran 90.
!
!$$$
    implicit none
    integer  :: jx
    real(r8) :: xmin,xmax,xinc,x,t
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsi-1)
    c2xpvsi=1.0_r8/xinc
    c1xpvsi=1.0_r8-xmin*c2xpvsi
    do jx=1,nxpvsi
      x=xmin+(jx-1)*xinc
      t=x
      tbpvsi(jx)=fpvsix(t)
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine gpvsi


!-------------------------------------------------------------------------------
  elemental function fpvsix(t)
!$$$     Subprogram Documentation Block
!
! Subprogram: fpvsix       Compute saturation vapor pressure over ice
!   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
!
! Abstract: Exactly compute saturation vapor pressure from temperature.
!   The water model assumes a perfect gas, constant specific heats
!   for gas and ice, and neglects the volume of the ice.
!   The model does account for the variation of the latent heat
!   of condensation with temperature.  The liquid option is not included.
!   The Clausius-Clapeyron equation is integrated from the triple point
!   to get the formula
!       pvsi=con_psat*(tr**xa)*exp(xb*(1.-tr))
!   where tr is ttp/t and other values are physical constants.
!   This function should be expanded inline in the calling routine.
!
! Program History Log:
!   91-05-07  Iredell             made into inlinable function
!   94-12-30  Iredell             exact computation
! 1999-03-01  Iredell             f90 module
! 2001-02-26  Iredell             ice phase
!
! Usage:   pvsi=fpvsix(t)
!
!   Input argument list:
!     t          Real(r8) temperature in Kelvin
!
!   Output argument list:
!     fpvsix     Real(r8) saturation vapor pressure in Pascals
!
! Attributes:
!   Language: Fortran 90.
!
!$$$
    implicit none
    real(r8) fpvsix
    real(r8),intent(in):: t
    real(r8),parameter:: dldt=con_cvap-con_csol
    real(r8),parameter:: heat=con_hvap+con_hfus
    real(r8),parameter:: xpona=-dldt/con_rv
    real(r8),parameter:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    real(r8) tr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvsix=con_psat*(tr**xpona)*exp(xponb*(1.0_r8-tr))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end function fpvsix
  
END  MODULE Micro_Ferrier
!PROGRAM Main
! USE Micro_Ferrier
!END PROGRAM Main
