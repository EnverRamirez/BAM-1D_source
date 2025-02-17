MODULE PhysicalFunctions
  IMPLICIT NONE
  PRIVATE


  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  ! Diagnostics: CAPE / CIN / LCL / LFC
  ! Depend on 2D/3D flag send in
  ! Code originally from RIP4

  REAL(KIND=r8) :: psadithte (1:150)
  REAL(KIND=r8) :: psadiprs  (1:150)
  REAL(KIND=r8) :: psaditmk  (150,150) 


  !$$$  Module Documentation Block
  !module module_gfs_funcphys
  !
  ! Module:    funcphys        API for basic thermodynamic physics
  !   Author: Iredell          Org: W/NX23     Date: 1999-03-01
  !
  ! Abstract: This module provides an Application Program Interface
  !   for computing basic thermodynamic physics functions, in particular
  !     (1) saturation vapor pressure as a function of temperature,
  !     (2) dewpoint temperature as a function of vapor pressure,
  !     (3) equivalent potential temperature as a function of temperature
  !         and scaled pressure to the kappa power,
  !     (4) temperature and specific humidity along a moist adiabat
  !         as functions of equivalent potential temperature and
  !         scaled pressure to the kappa power,
  !     (5) scaled pressure to the kappa power as a function of pressure, and
  !     (6) temperature at the lifting condensation level as a function
  !         of temperature and dewpoint depression.
  !   The entry points required to set up lookup tables start with a "g".
  !   All the other entry points are functions starting with an "f" or
  !   are subroutines starting with an "s".  These other functions and
  !   subroutines are elemental; that is, they return a scalar if they
  !   are passed only scalars, but they return an array if they are passed
  !   an array.  These other functions and subroutines can be inlined, too.
  !   
  ! Program History Log:
  !   1999-03-01  Mark Iredell
  !   1999-10-15  Mark Iredell  SI unit for pressure (Pascals)
  !   2001-02-26  Mark Iredell  Ice phase changes of Hong and Moorthi
  !
  ! Public Variables:
  !   r8         Integer parameter kind or length of reals (=r8)
  !
  ! Public Subprograms:
  !   gpvsl            Compute saturation vapor pressure over liquid table
  !
  !   fpvsl           Elementally compute saturation vapor pressure over liquid
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvslq          Elementally compute saturation vapor pressure over liquid
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvslx          Elementally compute saturation vapor pressure over liquid
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   gpvsi            Compute saturation vapor pressure over ice table
  !
  !   fpvsi           Elementally compute saturation vapor pressure over ice
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvsiq          Elementally compute saturation vapor pressure over ice
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvsix          Elementally compute saturation vapor pressure over ice
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   gpvs            Compute saturation vapor pressure table
  !
  !   fpvs            Elementally compute saturation vapor pressure
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvsq           Elementally compute saturation vapor pressure
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   fpvsx           Elementally compute saturation vapor pressure
  !     function result Real(r8) saturation vapor pressure in Pascals
  !     t               Real(r8) temperature in Kelvin
  !
  !   gtdpl           Compute dewpoint temperature over liquid table
  !
  !   ftdpl           Elementally compute dewpoint temperature over liquid
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdplq          Elementally compute dewpoint temperature over liquid
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdplx          Elementally compute dewpoint temperature over liquid
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdplxg         Elementally compute dewpoint temperature over liquid
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     t               Real(r8) guess dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   gtdpi           Compute dewpoint temperature table over ice
  !
  !   ftdpi           Elementally compute dewpoint temperature over ice
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpiq          Elementally compute dewpoint temperature over ice
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpix          Elementally compute dewpoint temperature over ice
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpixg         Elementally compute dewpoint temperature over ice
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     t               Real(r8) guess dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   gtdp            Compute dewpoint temperature table
  !
  !   ftdp            Elementally compute dewpoint temperature
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpq           Elementally compute dewpoint temperature
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpx           Elementally compute dewpoint temperature
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   ftdpxg          Elementally compute dewpoint temperature
  !     function result Real(r8) dewpoint temperature in Kelvin
  !     t               Real(r8) guess dewpoint temperature in Kelvin
  !     pv              Real(r8) vapor pressure in Pascals
  !
  !   gthe            Compute equivalent potential temperature table
  !
  !   fthe            Elementally compute equivalent potential temperature
  !     function result Real(r8) equivalent potential temperature in Kelvin
  !     t               Real(r8) LCL temperature in Kelvin
  !     pk              Real(r8) LCL pressure over 1e5 Pa to the kappa power
  !
  !   ftheq           Elementally compute equivalent potential temperature
  !     function result Real(r8) equivalent potential temperature in Kelvin
  !     t               Real(r8) LCL temperature in Kelvin
  !     pk              Real(r8) LCL pressure over 1e5 Pa to the kappa power
  !
  !   fthex           Elementally compute equivalent potential temperature
  !     function result Real(r8) equivalent potential temperature in Kelvin
  !     t               Real(r8) LCL temperature in Kelvin
  !     pk              Real(r8) LCL pressure over 1e5 Pa to the kappa power
  !
  !   gtma            Compute moist adiabat tables
  !
  !   stma            Elementally compute moist adiabat temperature and moisture
  !     the             Real(r8) equivalent potential temperature in Kelvin
  !     pk              Real(r8) pressure over 1e5 Pa to the kappa power
  !     tma             Real(r8) parcel temperature in Kelvin
  !     qma             Real(r8) parcel specific humidity in kg/kg
  !
  !   stmaq           Elementally compute moist adiabat temperature and moisture
  !     the             Real(r8) equivalent potential temperature in Kelvin
  !     pk              Real(r8) pressure over 1e5 Pa to the kappa power
  !     tma             Real(r8) parcel temperature in Kelvin
  !     qma             Real(r8) parcel specific humidity in kg/kg
  !
  !   stmax           Elementally compute moist adiabat temperature and moisture
  !     the             Real(r8) equivalent potential temperature in Kelvin
  !     pk              Real(r8) pressure over 1e5 Pa to the kappa power
  !     tma             Real(r8) parcel temperature in Kelvin
  !     qma             Real(r8) parcel specific humidity in kg/kg
  !
  !   stmaxg          Elementally compute moist adiabat temperature and moisture
  !     tg              Real(r8) guess parcel temperature in Kelvin
  !     the             Real(r8) equivalent potential temperature in Kelvin
  !     pk              Real(r8) pressure over 1e5 Pa to the kappa power
  !     tma             Real(r8) parcel temperature in Kelvin
  !     qma             Real(r8) parcel specific humidity in kg/kg
  !
  !   gpkap           Compute pressure to the kappa table
  !
  !   fpkap           Elementally raise pressure to the kappa power.
  !     function result Real(r8) p over 1e5 Pa to the kappa power
  !     p               Real(r8) pressure in Pascals
  !
  !   fpkapq          Elementally raise pressure to the kappa power.
  !     function result Real(r8) p over 1e5 Pa to the kappa power
  !     p               Real(r8) pressure in Pascals
  !
  !   fpkapo          Elementally raise pressure to the kappa power.
  !     function result Real(r8) p over 1e5 Pa to the kappa power
  !     p               Real(r8) surface pressure in Pascals
  !
  !   fpkapx          Elementally raise pressure to the kappa power.
  !     function result Real(r8) p over 1e5 Pa to the kappa power
  !     p               Real(r8) pressure in Pascals
  !
  !   grkap           Compute pressure to the 1/kappa table
  !
  !   frkap           Elementally raise pressure to the 1/kappa power.
  !     function result Real(r8) pressure in Pascals
  !     pkap            Real(r8) p over 1e5 Pa to the 1/kappa power
  !
  !   frkapq          Elementally raise pressure to the kappa power.
  !     function result Real(r8) pressure in Pascals
  !     pkap            Real(r8) p over 1e5 Pa to the kappa power
  !
  !   frkapx          Elementally raise pressure to the kappa power.
  !     function result Real(r8) pressure in Pascals
  !     pkap            Real(r8) p over 1e5 Pa to the kappa power
  !
  !   gtlcl           Compute LCL temperature table
  !
  !   ftlcl           Elementally compute LCL temperature.
  !     function result Real(r8) temperature at the LCL in Kelvin
  !     t               Real(r8) temperature in Kelvin
  !     tdpd            Real(r8) dewpoint depression in Kelvin
  !
  !   ftlclq          Elementally compute LCL temperature.
  !     function result Real(r8) temperature at the LCL in Kelvin
  !     t               Real(r8) temperature in Kelvin
  !     tdpd            Real(r8) dewpoint depression in Kelvin
  !
  !   ftlclo          Elementally compute LCL temperature.
  !     function result Real(r8) temperature at the LCL in Kelvin
  !     t               Real(r8) temperature in Kelvin
  !     tdpd            Real(r8) dewpoint depression in Kelvin
  !
  !   ftlclx          Elementally compute LCL temperature.
  !     function result Real(r8) temperature at the LCL in Kelvin
  !     t               Real(r8) temperature in Kelvin
  !     tdpd            Real(r8) dewpoint depression in Kelvin
  !
  !   gfuncphys       Compute all physics function tables
  !
  ! Attributes:
  !   Language: Fortran 90
  !
  !$$$

  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Private Variables
  !module module_gfs_physcons
  !  use module_gfs_machine,only:r8
  !  Physical constants as set in NMC handbook from Smithsonian tables.
  !  Physical constants are given to 5 places.
  !  1990/04/30: g and rd are made consistent with NWS usage.
  !  2001/10/22: g made consistent with SI usage.
  !  Math constants
  REAL(kind=r8),PARAMETER:: con_pi      =3.1415926535897931 ! pi
  REAL(kind=r8),PARAMETER:: con_sqrt2   =1.414214e+0 ! square root of 2
  REAL(kind=r8),PARAMETER:: con_sqrt3   =1.732051e+0 ! square root of 3
  !  Primary constants
  REAL(kind=r8),PARAMETER:: con_rerth   =6.3712e+6 ! radius of earth     (m)
  REAL(kind=r8),PARAMETER:: con_g       =9.80665e+0! gravity             (m/s2)
  REAL(kind=r8),PARAMETER:: con_omega   =7.2921e-5 ! ang vel of earth    (1/s)
  REAL(kind=r8),PARAMETER:: con_rd      =2.8705e+2 ! gas constant air    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_rv      =4.6150e+2 ! gas constant H2O    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cp      =1.0046e+3 ! spec heat air @p    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cv      =7.1760e+2 ! spec heat air @v    (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cvap    =1.8460e+3 ! spec heat H2O gas   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_cliq    =4.1855e+3 ! spec heat H2O liq   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_csol    =2.1060e+3 ! spec heat H2O ice   (J/kg/K)
  REAL(kind=r8),PARAMETER:: con_hvap    =2.5000e+6 ! lat heat H2O cond   (J/kg)
  REAL(kind=r8),PARAMETER:: con_hfus    =3.3358e+5 ! lat heat H2O fusion (J/kg)
  REAL(kind=r8),PARAMETER:: con_psat    =6.1078e+2 ! pres at H2O 3pt     (Pa)  
  REAL(kind=r8),PARAMETER:: con_sbc     =5.6730e-8 ! stefan-boltzmann    (W/m2/K4)
  REAL(kind=r8),PARAMETER:: con_solr    =1.3533e+3 ! solar constant      (W/m2)
  REAL(kind=r8),PARAMETER:: con_t0c     =2.7315e+2 ! temp at 0C          (K)
  REAL(kind=r8),PARAMETER:: con_ttp     =2.7316e+2 ! temp at H2O 3pt     (K)
  REAL(kind=r8),PARAMETER:: con_jcal    =4.1855E+0 ! JOULES PER CALORIE  ()
  !  Secondary constants
  REAL(kind=r8),PARAMETER:: con_rocp    =con_rd/con_cp
  REAL(kind=r8),PARAMETER:: con_cpor    =con_cp/con_rd
  REAL(kind=r8),PARAMETER:: con_rog     =con_rd/con_g
  REAL(kind=r8),PARAMETER:: con_fvirt   =con_rv/con_rd-1.
  REAL(kind=r8),PARAMETER:: con_eps     =con_rd/con_rv
  REAL(kind=r8),PARAMETER:: con_epsm1   =con_rd/con_rv-1.
  REAL(kind=r8),PARAMETER:: con_dldt    =con_cvap-con_cliq
  REAL(kind=r8),PARAMETER:: con_xpona   =-con_dldt/con_rv
  REAL(kind=r8),PARAMETER:: con_xponb   =-con_dldt/con_rv+con_hvap/(con_rv*con_ttp)
  !end module module_gfs_physcons

  REAL(r8),PARAMETER :: psatb=con_psat*1.e-5
  INTEGER ,PARAMETER :: nxpvsl=7501
  REAL(r8)           :: c1xpvsl,c2xpvsl,tbpvsl(nxpvsl)
  INTEGER ,PARAMETER :: nxpvsi=7501
  REAL(r8)           :: c1xpvsi,c2xpvsi,tbpvsi(nxpvsi)
  INTEGER ,PARAMETER :: nxpvs=7501
  REAL(r8)           :: c1xpvs,c2xpvs,tbpvs(nxpvs)
  INTEGER ,PARAMETER :: nxtdpl=5001
  REAL(r8)           :: c1xtdpl,c2xtdpl,tbtdpl(nxtdpl)
  INTEGER ,PARAMETER :: nxtdpi=5001
  REAL(r8)           :: c1xtdpi,c2xtdpi,tbtdpi(nxtdpi)
  INTEGER ,PARAMETER :: nxtdp=5001
  REAL(r8)           :: c1xtdp,c2xtdp,tbtdp(nxtdp)
  INTEGER ,PARAMETER :: nxthe=241,nythe=151
  REAL(r8)           :: c1xthe,c2xthe,c1ythe,c2ythe,tbthe(nxthe,nythe)
  INTEGER ,PARAMETER :: nxma=151,nyma=121
  REAL(r8)           :: c1xma,c2xma,c1yma,c2yma,tbtma(nxma,nyma),tbqma(nxma,nyma)
  ! integer,parameter  :: nxpkap=5501
  INTEGER ,PARAMETER :: nxpkap=11001
  REAL(r8)           :: c1xpkap,c2xpkap,tbpkap(nxpkap)
  INTEGER ,PARAMETER :: nxrkap=5501
  REAL(r8)           :: c1xrkap,c2xrkap,tbrkap(nxrkap)
  INTEGER ,PARAMETER :: nxtlcl=151,nytlcl=61
  REAL(r8)           :: c1xtlcl,c2xtlcl,c1ytlcl,c2ytlcl,tbtlcl(nxtlcl,nytlcl)


  INTEGER :: k850
  INTEGER :: k500
  INTEGER :: k700
  REAL(KIND=r8)            :: ae  (2)
  REAL(KIND=r8)            :: be  (2)
  REAL(KIND=r8)            :: ht  (2)

  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ! Public Subprograms
  !
  !
  PUBLIC :: gpvsl,fpvsl,fpvslq,fpvslx
  PUBLIC :: gpvsi,fpvsi,fpvsiq,fpvsix
  PUBLIC :: gpvs,fpvs,fpvsq,fpvsx,fpvs2es5
  PUBLIC :: gtdpl,ftdpl,ftdplq,ftdplx,ftdplxg
  PUBLIC :: gtdpi,ftdpi,ftdpiq,ftdpix,ftdpixg
  PUBLIC :: gtdp,ftdp,ftdpq,ftdpx,ftdpxg
  PUBLIC :: gthe,fthe,ftheq,fthex
  PUBLIC :: gtma,stma,stmaq,stmax,stmaxg
  PUBLIC :: gpkap,fpkap,fpkapq,fpkapo,fpkapx
  PUBLIC :: grkap,frkap,frkapq,frkapx
  PUBLIC :: gtlcl,ftlcl,ftlclq,ftlclo,ftlclx
  PUBLIC :: gfuncphys  
  PUBLIC :: InitPhysicalFunctions
  PUBLIC :: calc_cape
  PUBLIC :: SWEAT_index
CONTAINS

  SUBROUTINE InitPhysicalFunctions(si,sl,kMax)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: si(kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: sl(kMax)
    !
    ht(1)=con_hvap/con_cp
    
    ht(2)=2.834e6_r8/con_cp
    
    be(1)=0.622_r8*ht(1)/0.286_r8
    
    ae(1)=be(1)/273.0_r8+LOG(610.71_r8)
    
    be(2)=0.622_r8*ht(2)/0.286_r8
    
    ae(2)=be(2)/273.0_r8+LOG(610.71_r8)

    CALL gfuncphys()
    CALL RD()
    CALL INDEXLEVEL(sl,kMax)
  END SUBROUTINE InitPhysicalFunctions


  SUBROUTINE INDEXLEVEL(sl,kMax)
    IMPLICIT NONE
    INTEGER      , INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: sl(kMax)
    REAL(KIND=r8) :: Dsl(kMax)
    INTEGER       :: level(1)
    INTEGER       :: k
    Dsl=1.0e20_r8
    DO k=1,kMax
       Dsl(k)=abs(sl(k)-0.850_r8)
    END DO
    level=MINLOC(Dsl,DIM=1)
    k850=level(1)
    Dsl=1.0e20_r8
    DO k=1,kMax
       Dsl(k)=abs(sl(k)-0.500_r8)
    END DO
    level=MINLOC(Dsl,DIM=1)
    k500=level(1)
    Dsl=1.0e20_r8
    DO k=1,kMax
       Dsl(k)=abs(sl(k)-0.700_r8)
    END DO
    level=MINLOC(Dsl,DIM=1)
    k700=level(1)
  END SUBROUTINE INDEXLEVEL

  SUBROUTINE gfuncphys()
    IMPLICIT NONE
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gfuncphys    Compute all physics function tables
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute all physics function tables.  Lookup tables are
    !   set up for computing saturation vapor pressure, dewpoint temperature,
    !   equivalent potential temperature, moist adiabatic temperature and humidity,
    !   pressure to the kappa, and lifting condensation level temperature.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gfuncphys
    !
    ! Subprograms called:
    !   gpvsl       compute saturation vapor pressure over liquid table
    !   gpvsi       compute saturation vapor pressure over ice table
    !   gpvs        compute saturation vapor pressure table
    !   gtdpl       compute dewpoint temperature over liquid table
    !   gtdpi       compute dewpoint temperature over ice table
    !   gtdp        compute dewpoint temperature table
    !   gthe        compute equivalent potential temperature table
    !   gtma        compute moist adiabat tables
    !   gpkap       compute pressure to the kappa table
    !   grkap       compute pressure to the 1/kappa table
    !   gtlcl       compute LCL temperature table
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    CALL gpvsl()
    CALL gpvsi()
    CALL gpvs()
    CALL gtdpl()
    CALL gtdpi()
    CALL gtdp()
    CALL gthe()
    CALL gtma()
    CALL gpkap()
    CALL grkap()
    CALL gtlcl()

  END SUBROUTINE gfuncphys


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
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsl-1)
    !   c1xpvsl=1.-xmin/xinc
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
  ! elemental function fpvsl(t)
  FUNCTION fpvsl(t)
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
    REAL(r8)             :: fpvsl
    REAL(r8),INTENT(in)  :: t
    INTEGER              :: jx
    REAL(r8)             :: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsl+c2xpvsl*t,1.0_r8),REAL(nxpvsl,r8))
    jx=MIN(xj,nxpvsl-1.0_r8)
    fpvsl=tbpvsl(jx)+(xj-jx)*(tbpvsl(jx+1)-tbpvsl(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsl
  !-------------------------------------------------------------------------------
  ! elemental function fpvslq(t)
  FUNCTION fpvslq(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvslq       Compute saturation vapor pressure over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gpvsl. See documentation for fpvslx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 9 decimal places.
    !   On the Cray, fpvslq is about 3 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   pvsl=fpvslq(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvslq     Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)            :: fpvslq
    REAL(r8),INTENT(in) :: t
    INTEGER             :: jx
    REAL(r8)            :: xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsl+c2xpvsl*t,1.0_r8),REAL(nxpvsl,r8))
    jx=MIN(MAX(NINT(xj),2),nxpvsl-1)
    dxj=xj-jx
    fj1=tbpvsl(jx-1)
    fj2=tbpvsl(jx)
    fj3=tbpvsl(jx+1)
    fpvslq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvslq
  !-------------------------------------------------------------------------------
  ! elemental function fpvslx(t)
  FUNCTION fpvslx(t)
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
    REAL(r8) fpvslx
    REAL(r8),INTENT(in):: t
    REAL(r8),PARAMETER:: dldt=con_cvap-con_cliq
    REAL(r8),PARAMETER:: heat=con_hvap
    REAL(r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(r8) tr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvslx=con_psat*(tr**xpona)*EXP(xponb*(1.0-tr))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvslx
  !-------------------------------------------------------------------------------
  SUBROUTINE gpvsi
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
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,x,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180.0_r8
    xmax=330.0_r8
    xinc=(xmax-xmin)/(nxpvsi-1)
    !   c1xpvsi=1.-xmin/xinc
    c2xpvsi=1.0/xinc
    c1xpvsi=1.0-xmin*c2xpvsi
    DO jx=1,nxpvsi
       x=xmin+(jx-1)*xinc
       t=x
       tbpvsi(jx)=fpvsix(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvsi
  !-------------------------------------------------------------------------------
  ! elemental function fpvsi(t)
  FUNCTION fpvsi(t)
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
    IMPLICIT NONE
    REAL(r8)            :: fpvsi
    REAL(r8),INTENT(in) :: t
    INTEGER             :: jx
    REAL(r8)            :: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsi+c2xpvsi*t,1.0_r8),REAL(nxpvsi,r8))
    jx=MIN(xj,nxpvsi-1.0_r8)
    fpvsi=tbpvsi(jx)+(xj-jx)*(tbpvsi(jx+1)-tbpvsi(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsi
  !-------------------------------------------------------------------------------
  ! elemental function fpvsiq(t)
  FUNCTION fpvsiq(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsiq       Compute saturation vapor pressure over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gpvsi. See documentation for fpvsix for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 9 decimal places.
    !   On the Cray, fpvsiq is about 3 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvsi=fpvsiq(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsiq     Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fpvsiq
    REAL(r8),INTENT(in):: t
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvsi+c2xpvsi*t,1.0_r8),REAL(nxpvsi,r8))
    jx=MIN(MAX(NINT(xj),2),nxpvsi-1)
    dxj=xj-jx
    fj1=tbpvsi(jx-1)
    fj2=tbpvsi(jx)
    fj3=tbpvsi(jx+1)
    fpvsiq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsiq
  !-------------------------------------------------------------------------------
  ! elemental function fpvsix(t)
  FUNCTION fpvsix(t)
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
    IMPLICIT NONE
    REAL(r8)            :: fpvsix
    REAL(r8),INTENT(in) :: t
    REAL(r8),PARAMETER  :: dldt=con_cvap-con_csol
    REAL(r8),PARAMETER  :: heat=con_hvap+con_hfus
    REAL(r8),PARAMETER  :: xpona=-dldt/con_rv
    REAL(r8),PARAMETER  :: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(r8) tr
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    fpvsix=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsix
  !-------------------------------------------------------------------------------
  SUBROUTINE gpvs
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gpvs         Compute saturation vapor pressure table
    !   Author: N Phillips            W/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Computes saturation vapor pressure table as a function of
    !   temperature for the table lookup function fpvs.
    !   Exact saturation vapor pressures are calculated in subprogram fpvsx.
    !   The current implementation computes a table with a length
    !   of 7501 for temperatures ranging from 180. to 330. Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gpvs
    !
    ! Subprograms called:
    !   (fpvsx)    inlinable function to compute saturation vapor pressure
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
    xinc=(xmax-xmin)/(nxpvs-1)
    !   c1xpvs=1.-xmin/xinc
    c2xpvs=1.0_r8/xinc
    c1xpvs=1.0_r8-xmin*c2xpvs
    DO jx=1,nxpvs
       x=xmin+(jx-1)*xinc
       t=x
       tbpvs(jx)=fpvsx(t)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpvs

  !---------------------------------
  REAL(KIND=r8) FUNCTION es5(t)
    REAL(KIND=r8), INTENT(IN) :: t !K
    REAL(KIND=r8)   , PARAMETER :: tcrit    =   273.15_r8
    !
    IF (t <= tcrit) THEN
       es5 = EXP(ae(2)-be(2)/t)!Pa
    ELSE
       es5 = EXP(ae(1)-be(1)/t)!Pa
    END IF
  END FUNCTION es5
  !-------------------------------------------------------------------------------
  ! elemental function fpvs(t)
  FUNCTION fpvs2es5(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvs         Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvs is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvs(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvs       Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)            :: fpvs2es5
    REAL(r8),INTENT(in) :: t
    INTEGER             :: jx
    REAL(r8)            :: xj
    IF(t>270.0_r8)THEN
!erg 4/8/2015    IF(t>237.0_r8)THEN  !t>237 it's above than spontaneous freezing level
      fpvs2es5=es5(t)
    ELSE
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=MIN(xj,nxpvs-1.0_r8)
    fpvs2es5=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    END IF
  END FUNCTION fpvs2es5

  !-------------------------------------------------------------------------------
  ! elemental function fpvs(t)
  FUNCTION fpvs(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvs         Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 6 decimal places.
    !   On the Cray, fpvs is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvs(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvs       Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)            :: fpvs
    REAL(r8),INTENT(in) :: t
    INTEGER             :: jx
    REAL(r8)            :: xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=MIN(xj,nxpvs-1.0_r8)
    fpvs=tbpvs(jx)+(xj-jx)*(tbpvs(jx+1)-tbpvs(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvs
  !-------------------------------------------------------------------------------
  ! elemental function fpvsq(t)
  FUNCTION fpvsq(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsq        Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute saturation vapor pressure from the temperature.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gpvs. See documentation for fpvsx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is almost 9 decimal places.
    !   On the Cray, fpvsq is about 3 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvsq(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsq      Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)           :: fpvsq
    REAL(r8),INTENT(in):: t
    INTEGER            :: jx
    REAL(r8)           :: xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpvs+c2xpvs*t,1.0_r8),REAL(nxpvs,r8))
    jx=MIN(MAX(NINT(xj),2),nxpvs-1)
    dxj=xj-jx
    fj1=tbpvs(jx-1)
    fj2=tbpvs(jx)
    fj3=tbpvs(jx+1)
    fpvsq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsq
  !-------------------------------------------------------------------------------
  ! elemental function fpvsx(t)
  FUNCTION fpvsx(t)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpvsx        Compute saturation vapor pressure
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute saturation vapor pressure from temperature.
    !   The saturation vapor pressure over either liquid and ice is computed
    !   over liquid for temperatures above the triple point,
    !   over ice for temperatures 20 degress below the triple point,
    !   and a linear combination of the two for temperatures in between.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas, liquid and ice, and neglects the volume of the condensate.
    !   The model does account for the variation of the latent heat
    !   of condensation and sublimation with temperature.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are ph:ysical constants.
    !   The reference for this computation is Emanuel(1994), pages 116-117.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   pvs=fpvsx(t)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !
    !   Output argument list:
    !     fpvsx      Real(r8) saturation vapor pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8)            :: fpvsx
    REAL(r8),INTENT(in) :: t
    REAL(r8),PARAMETER  :: tliq=con_ttp ! 2.7316e+2  temp at H2O 3pt     (K)
    REAL(r8),PARAMETER  :: tice=con_ttp-20.0_r8
    REAL(r8),PARAMETER  :: dldtl=con_cvap-con_cliq
    REAL(r8),PARAMETER  :: heatl=con_hvap
    REAL(r8),PARAMETER  :: xponal=-dldtl/con_rv
    REAL(r8),PARAMETER  :: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
    REAL(r8),PARAMETER  :: dldti=con_cvap-con_csol
    REAL(r8),PARAMETER  :: heati=con_hvap+con_hfus
    REAL(r8),PARAMETER  :: xponai=-dldti/con_rv
    REAL(r8),PARAMETER  :: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
    REAL(r8) tr,w,pvl,pvi
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/t
    IF(t.GE.tliq) THEN
       fpvsx=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
    ELSEIF(t.LT.tice) THEN
       fpvsx=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
    ELSE
       w=(t-tice)/(tliq-tice)
       pvl=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
       pvi=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
       fpvsx=w*pvl+(1.0_r8-w)*pvi
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpvsx
  !-------------------------------------------------------------------------------
  SUBROUTINE gtdpl
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gtdpl        Compute dewpoint temperature over liquid table
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature table as a function of
    !   vapor pressure for inlinable function ftdpl.
    !   Exact dewpoint temperatures are calculated in subprogram ftdplxg.
    !   The current implementation computes a table with a length
    !   of 5001 for vapor pressures ranging from 1 to 10001 Pascals
    !   giving a dewpoint temperature range of 208 to 319 Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gtdpl
    !
    ! Subprograms called:
    !   (ftdplxg)  inlinable function to compute dewpoint temperature over liquid
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,t,x,pv
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=1
    xmax=10001
    xinc=(xmax-xmin)/(nxtdpl-1)
    c1xtdpl=1.0_r8-xmin/xinc
    c2xtdpl=1.0_r8/xinc
    t=208.00_r8
    DO jx=1,nxtdpl
       x=xmin+(jx-1)*xinc
       pv=x
       t=ftdplxg(t,pv)
       tbtdpl(jx)=t
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gtdpl
  !-------------------------------------------------------------------------------
  ! elemental function ftdpl(pv)
  FUNCTION ftdpl(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpl        Compute dewpoint temperature over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gtdpl. See documentation for ftdplxg for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.0005 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.02 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdpl is about 75 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tdpl=ftdpl(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpl      Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpl
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdpl+c2xtdpl*pv,1.0_r8),REAL(nxtdpl,r8))
    jx=MIN(xj,nxtdpl-1.0_r8)
    ftdpl=tbtdpl(jx)+(xj-jx)*(tbtdpl(jx+1)-tbtdpl(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpl
  !-------------------------------------------------------------------------------
  ! elemental function ftdplq(pv)
  FUNCTION ftdplq(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdplq       Compute dewpoint temperature over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gtdpl. see documentation for ftdplxg for details.
    !   Input values outside table range are reset to table extrema.
    !   the interpolation accuracy is better than 0.00001 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.002 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdplq is about 60 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tdpl=ftdplq(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdplq     Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdplq
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdpl+c2xtdpl*pv,1.0_r8),REAL(nxtdpl,r8))
    jx=MIN(MAX(NINT(xj),2),nxtdpl-1)
    dxj=xj-jx
    fj1=tbtdpl(jx-1)
    fj2=tbtdpl(jx)
    fj3=tbtdpl(jx+1)
    ftdplq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdplq
  !-------------------------------------------------------------------------------
  ! elemental function ftdplx(pv)
  FUNCTION ftdplx(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdplx       Compute dewpoint temperature over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: exactly compute dewpoint temperature from vapor pressure.
    !   An approximate dewpoint temperature for function ftdplxg
    !   is obtained using ftdpl so gtdpl must be already called.
    !   See documentation for ftdplxg for details.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tdpl=ftdplx(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdplx     Real(r8) dewpoint temperature in Kelvin
    !
    ! Subprograms called:
    !   (ftdpl)    inlinable function to compute dewpoint temperature over liquid
    !   (ftdplxg)  inlinable function to compute dewpoint temperature over liquid
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdplx
    REAL(r8),INTENT(in):: pv
    REAL(r8) tg
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tg=ftdpl(pv)
    ftdplx=ftdplxg(tg,pv)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdplx
  !-------------------------------------------------------------------------------
  ! elemental function ftdplxg(tg,pv)
  FUNCTION ftdplxg(tg,pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdplxg      Compute dewpoint temperature over liquid
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute dewpoint temperature from vapor pressure.
    !   A guess dewpoint temperature must be provided.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and liquid, and neglects the volume of the liquid.
    !   The model does account for the variation of the latent heat
    !   of condensation with temperature.  The ice option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvs=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The formula is inverted by iterating Newtonian approximations
    !   for each pvs until t is found to within 1.e-6 Kelvin.
    !   This function can be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tdpl=ftdplxg(tg,pv)
    !
    !   Input argument list:
    !     tg         Real(r8) guess dewpoint temperature in Kelvin
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdplxg    Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdplxg
    REAL(r8),INTENT(in):: tg,pv
    REAL(r8),PARAMETER:: terrm=1.0e-6_r8
    REAL(r8),PARAMETER:: dldt=con_cvap-con_cliq
    REAL(r8),PARAMETER:: heat=con_hvap
    REAL(r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(r8) t,tr,pvt,el,dpvt,terr
    INTEGER i
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t=tg
    DO i=1,100
       tr=con_ttp/t
       pvt=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
       el=heat+dldt*(t-con_ttp)
       dpvt=el*pvt/(con_rv*t**2)
       terr=(pvt-pv)/dpvt
       t=t-terr
       IF(ABS(terr).LE.terrm) EXIT
    ENDDO
    ftdplxg=t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdplxg
  !-------------------------------------------------------------------------------
  SUBROUTINE gtdpi
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gtdpi        Compute dewpoint temperature over ice table
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature table as a function of
    !   vapor pressure for inlinable function ftdpi.
    !   Exact dewpoint temperatures are calculated in subprogram ftdpixg.
    !   The current implementation computes a table with a length
    !   of 5001 for vapor pressures ranging from 0.1 to 1000.1 Pascals
    !   giving a dewpoint temperature range of 197 to 279 Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gtdpi
    !
    ! Subprograms called:
    !   (ftdpixg)  inlinable function to compute dewpoint temperature over ice
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,t,x,pv
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=0.10_r8
    xmax=1000.10_r8
    xinc=(xmax-xmin)/(nxtdpi-1)
    c1xtdpi=1.0_r8-xmin/xinc
    c2xtdpi=1.0_r8/xinc
    t=197.00_r8
    DO jx=1,nxtdpi
       x=xmin+(jx-1)*xinc
       pv=x
       t=ftdpixg(t,pv)
       tbtdpi(jx)=t
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gtdpi
  !-------------------------------------------------------------------------------
  ! elemental function ftdpi(pv)
  FUNCTION ftdpi(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpi        Compute dewpoint temperature over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gtdpi. See documentation for ftdpixg for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.0005 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.02 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdpi is about 75 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdpi=ftdpi(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpi      Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpi
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdpi+c2xtdpi*pv,1.0_r8),REAL(nxtdpi,r8))
    jx=MIN(xj,nxtdpi-1.0_r8)
    ftdpi=tbtdpi(jx)+(xj-jx)*(tbtdpi(jx+1)-tbtdpi(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpi
  !-------------------------------------------------------------------------------
  ! elemental function ftdpiq(pv)
  FUNCTION ftdpiq(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpiq       Compute dewpoint temperature over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gtdpi. see documentation for ftdpixg for details.
    !   Input values outside table range are reset to table extrema.
    !   the interpolation accuracy is better than 0.00001 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.002 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdpiq is about 60 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdpi=ftdpiq(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpiq     Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpiq
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdpi+c2xtdpi*pv,1.0_r8),REAL(nxtdpi,r8))
    jx=MIN(MAX(NINT(xj),2),nxtdpi-1)
    dxj=xj-jx
    fj1=tbtdpi(jx-1)
    fj2=tbtdpi(jx)
    fj3=tbtdpi(jx+1)
    ftdpiq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpiq
  !-------------------------------------------------------------------------------
  ! elemental function ftdpix(pv)
  FUNCTION ftdpix(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpix       Compute dewpoint temperature over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: exactly compute dewpoint temperature from vapor pressure.
    !   An approximate dewpoint temperature for function ftdpixg
    !   is obtained using ftdpi so gtdpi must be already called.
    !   See documentation for ftdpixg for details.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdpi=ftdpix(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpix     Real(r8) dewpoint temperature in Kelvin
    !
    ! Subprograms called:
    !   (ftdpi)    inlinable function to compute dewpoint temperature over ice
    !   (ftdpixg)  inlinable function to compute dewpoint temperature over ice
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpix
    REAL(r8),INTENT(in):: pv
    REAL(r8) tg
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tg=ftdpi(pv)
    ftdpix=ftdpixg(tg,pv)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpix
  !-------------------------------------------------------------------------------
  ! elemental function ftdpixg(tg,pv)
  FUNCTION ftdpixg(tg,pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpixg      Compute dewpoint temperature over ice
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute dewpoint temperature from vapor pressure.
    !   A guess dewpoint temperature must be provided.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and ice, and neglects the volume of the ice.
    !   The model does account for the variation of the latent heat
    !   of sublimation with temperature.  The liquid option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvs=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The formula is inverted by iterating Newtonian approximations
    !   for each pvs until t is found to within 1.e-6 Kelvin.
    !   This function can be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdpi=ftdpixg(tg,pv)
    !
    !   Input argument list:
    !     tg         Real(r8) guess dewpoint temperature in Kelvin
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpixg    Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpixg
    REAL(r8),INTENT(in):: tg,pv
    REAL(r8),PARAMETER:: terrm=1.e-6_r8
    REAL(r8),PARAMETER:: dldt=con_cvap-con_csol
    REAL(r8),PARAMETER:: heat=con_hvap+con_hfus
    REAL(r8),PARAMETER:: xpona=-dldt/con_rv
    REAL(r8),PARAMETER:: xponb=-dldt/con_rv+heat/(con_rv*con_ttp)
    REAL(r8) t,tr,pvt,el,dpvt,terr
    INTEGER i
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t=tg
    DO i=1,100
       tr=con_ttp/t
       pvt=con_psat*(tr**xpona)*EXP(xponb*(1.0_r8-tr))
       el=heat+dldt*(t-con_ttp)
       dpvt=el*pvt/(con_rv*t**2)
       terr=(pvt-pv)/dpvt
       t=t-terr
       IF(ABS(terr).LE.terrm) EXIT
    ENDDO
    ftdpixg=t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpixg
  !-------------------------------------------------------------------------------
  SUBROUTINE gtdp
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gtdp         Compute dewpoint temperature table
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature table as a function of
    !   vapor pressure for inlinable function ftdp.
    !   Exact dewpoint temperatures are calculated in subprogram ftdpxg.
    !   The current implementation computes a table with a length
    !   of 5001 for vapor pressures ranging from 0.5 to 1000.5 Pascals
    !   giving a dewpoint temperature range of 208 to 319 Kelvin.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:  call gtdp
    !
    ! Subprograms called:
    !   (ftdpxg)   inlinable function to compute dewpoint temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,t,x,pv
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=0.50_r8
    xmax=10000.50_r8
    xinc=(xmax-xmin)/(nxtdp-1)
    c1xtdp=1.0-xmin/xinc
    c2xtdp=1.0/xinc
    t=208.00_r8
    DO jx=1,nxtdp
       x=xmin+(jx-1)*xinc
       pv=x
       t=ftdpxg(t,pv)
       tbtdp(jx)=t
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gtdp
  !-------------------------------------------------------------------------------
  ! elemental function ftdp(pv)
  FUNCTION ftdp(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdp         Compute dewpoint temperature 
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gtdp. See documentation for ftdpxg for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.0005 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.02 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdp is about 75 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdp=ftdp(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdp       Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdp
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdp+c2xtdp*pv,1.0_r8),REAL(nxtdp,r8))
    jx=MIN(xj,nxtdp-1.0_r8)
    ftdp=tbtdp(jx)+(xj-jx)*(tbtdp(jx+1)-tbtdp(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdp
  !-------------------------------------------------------------------------------
  ! elemental function ftdpq(pv)
  FUNCTION ftdpq(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpq        Compute dewpoint temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute dewpoint temperature from vapor pressure.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gtdp. see documentation for ftdpxg for details.
    !   Input values outside table range are reset to table extrema.
    !   the interpolation accuracy is better than 0.00001 Kelvin
    !   for dewpoint temperatures greater than 250 Kelvin,
    !   but decreases to 0.002 Kelvin for a dewpoint around 230 Kelvin.
    !   On the Cray, ftdpq is about 60 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdp=ftdpq(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpq      Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpq
    REAL(r8),INTENT(in):: pv
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtdp+c2xtdp*pv,1.0_r8),REAL(nxtdp,r8))
    jx=MIN(MAX(NINT(xj),2),nxtdp-1)
    dxj=xj-jx
    fj1=tbtdp(jx-1)
    fj2=tbtdp(jx)
    fj3=tbtdp(jx+1)
    ftdpq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpq
  !-------------------------------------------------------------------------------
  ! elemental function ftdpx(pv)
  FUNCTION ftdpx(pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpx        Compute dewpoint temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: exactly compute dewpoint temperature from vapor pressure.
    !   An approximate dewpoint temperature for function ftdpxg
    !   is obtained using ftdp so gtdp must be already called.
    !   See documentation for ftdpxg for details.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdp=ftdpx(pv)
    !
    !   Input argument list:
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpx      Real(r8) dewpoint temperature in Kelvin
    !
    ! Subprograms called:
    !   (ftdp)     inlinable function to compute dewpoint temperature
    !   (ftdpxg)   inlinable function to compute dewpoint temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpx
    REAL(r8),INTENT(in):: pv
    REAL(r8) tg
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tg=ftdp(pv)
    ftdpx=ftdpxg(tg,pv)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpx
  !-------------------------------------------------------------------------------
  ! elemental function ftdpxg(tg,pv)
  FUNCTION ftdpxg(tg,pv)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftdpxg       Compute dewpoint temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute dewpoint temperature from vapor pressure.
    !   A guess dewpoint temperature must be provided.
    !   The saturation vapor pressure over either liquid and ice is computed
    !   over liquid for temperatures above the triple point,
    !   over ice for temperatures 20 degress below the triple point,
    !   and a linear combination of the two for temperatures in between.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas, liquid and ice, and neglects the volume of the condensate.
    !   The model does account for the variation of the latent heat
    !   of condensation and sublimation with temperature.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formula
    !       pvsl=con_psat*(tr**xa)*exp(xb*(1.-tr))
    !   where tr is ttp/t and other values are physical constants.
    !   The reference for this decision is Emanuel(1994), pages 116-117.
    !   The formula is inverted by iterating Newtonian approximations
    !   for each pvs until t is found to within 1.e-6 Kelvin.
    !   This function can be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    ! 2001-02-26  Iredell             ice phase
    !
    ! Usage:   tdp=ftdpxg(tg,pv)
    !
    !   Input argument list:
    !     tg         Real(r8) guess dewpoint temperature in Kelvin
    !     pv         Real(r8) vapor pressure in Pascals
    !
    !   Output argument list:
    !     ftdpxg     Real(r8) dewpoint temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftdpxg
    REAL(r8),INTENT(in):: tg,pv
    REAL(r8),PARAMETER:: terrm=1.e-6_r8
    REAL(r8),PARAMETER:: tliq=con_ttp
    REAL(r8),PARAMETER:: tice=con_ttp-20.00_r8
    REAL(r8),PARAMETER:: dldtl=con_cvap-con_cliq
    REAL(r8),PARAMETER:: heatl=con_hvap
    REAL(r8),PARAMETER:: xponal=-dldtl/con_rv
    REAL(r8),PARAMETER:: xponbl=-dldtl/con_rv+heatl/(con_rv*con_ttp)
    REAL(r8),PARAMETER:: dldti=con_cvap-con_csol
    REAL(r8),PARAMETER:: heati=con_hvap+con_hfus
    REAL(r8),PARAMETER:: xponai=-dldti/con_rv
    REAL(r8),PARAMETER:: xponbi=-dldti/con_rv+heati/(con_rv*con_ttp)
    REAL(r8) t,tr,w,pvtl,pvti,pvt,ell,eli,el,dpvt,terr
    INTEGER i
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t=tg
    DO i=1,100
       tr=con_ttp/t
       IF(t.GE.tliq) THEN
          pvt=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
          el=heatl+dldtl*(t-con_ttp)
          dpvt=el*pvt/(con_rv*t**2)
       ELSEIF(t.LT.tice) THEN
          pvt=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
          el=heati+dldti*(t-con_ttp)
          dpvt=el*pvt/(con_rv*t**2)
       ELSE
          w=(t-tice)/(tliq-tice)
          pvtl=con_psat*(tr**xponal)*EXP(xponbl*(1.0_r8-tr))
          pvti=con_psat*(tr**xponai)*EXP(xponbi*(1.0_r8-tr))
          pvt=w*pvtl+(1.0_r8-w)*pvti
          ell=heatl+dldtl*(t-con_ttp)
          eli=heati+dldti*(t-con_ttp)
          dpvt=(w*ell*pvtl+(1.0_r8-w)*eli*pvti)/(con_rv*t**2)
       ENDIF
       terr=(pvt-pv)/dpvt
       t=t-terr
       IF(ABS(terr).LE.terrm) EXIT
    ENDDO
    ftdpxg=t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftdpxg
  !-------------------------------------------------------------------------------
  SUBROUTINE gthe
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gthe        Compute equivalent potential temperature table
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute equivalent potential temperature table
    !   as a function of LCL temperature and pressure over 1e5 Pa
    !   to the kappa power for function fthe.
    !   Equivalent potential temperatures are calculated in subprogram fthex
    !   the current implementation computes a table with a first dimension
    !   of 241 for temperatures ranging from 183.16 to 303.16 Kelvin
    !   and a second dimension of 151 for pressure over 1e5 Pa
    !   to the kappa power ranging from 0.04**rocp to 1.10**rocp.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gthe
    !
    ! Subprograms called:
    !   (fthex)    inlinable function to compute equiv. pot. temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx,jy
    REAL(r8) xmin,xmax,ymin,ymax,xinc,yinc,x,y,pk,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=con_ttp-90.0_r8
    xmax=con_ttp+30.0_r8
    ymin=0.04_r8**con_rocp
    ymax=1.10_r8**con_rocp
    xinc=(xmax-xmin)/(nxthe-1)
    c1xthe=1.0_r8-xmin/xinc
    c2xthe=1.0_r8/xinc
    yinc=(ymax-ymin)/(nythe-1)
    c1ythe=1.0_r8-ymin/yinc
    c2ythe=1.0_r8/yinc
    DO jy=1,nythe
       y=ymin+(jy-1)*yinc
       pk=y
       DO jx=1,nxthe
          x=xmin+(jx-1)*xinc
          t=x
          tbthe(jx,jy)=fthex(t,pk)
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gthe
  !-------------------------------------------------------------------------------
  ! elemental function fthe(t,pk)
  FUNCTION fthe(t,pk)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fthe         Compute equivalent potential temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute equivalent potential temperature at the LCL
    !   from temperature and pressure over 1e5 Pa to the kappa power.
    !   A bilinear interpolation is done between values in a lookup table
    !   computed in gthe. see documentation for fthex for details.
    !   Input values outside table range are reset to table extrema,
    !   except zero is returned for too cold or high LCLs.
    !   The interpolation accuracy is better than 0.01 Kelvin.
    !   On the Cray, fthe is almost 6 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   the=fthe(t,pk)
    !
    !   Input argument list:
    !     t          Real(r8) LCL temperature in Kelvin
    !     pk         Real(r8) LCL pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     fthe       Real(r8) equivalent potential temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fthe
    REAL(r8),INTENT(in):: t,pk
    INTEGER jx,jy
    REAL(r8) xj,yj,ftx1,ftx2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(c1xthe+c2xthe*t,REAL(nxthe,r8))
    yj=MIN(c1ythe+c2ythe*pk,REAL(nythe,r8))
    IF(xj.GE.1.0_r8.AND.yj.GE.1.0_r8) THEN
       jx=MIN(xj,nxthe-1.0_r8)
       jy=MIN(yj,nythe-1.0_r8)
       ftx1=tbthe(jx,jy)+(xj-jx)*(tbthe(jx+1,jy)-tbthe(jx,jy))
       ftx2=tbthe(jx,jy+1)+(xj-jx)*(tbthe(jx+1,jy+1)-tbthe(jx,jy+1))
       fthe=ftx1+(yj-jy)*(ftx2-ftx1)
    ELSE
       fthe=0.0_r8
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fthe
  !-------------------------------------------------------------------------------
  ! elemental function ftheq(t,pk)
  FUNCTION ftheq(t,pk)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftheq        Compute equivalent potential temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute equivalent potential temperature at the LCL
    !   from temperature and pressure over 1e5 Pa to the kappa power.
    !   A biquadratic interpolation is done between values in a lookup table
    !   computed in gthe. see documentation for fthex for details.
    !   Input values outside table range are reset to table extrema,
    !   except zero is returned for too cold or high LCLs.
    !   The interpolation accuracy is better than 0.0002 Kelvin.
    !   On the Cray, ftheq is almost 3 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   the=ftheq(t,pk)
    !
    !   Input argument list:
    !     t          Real(r8) LCL temperature in Kelvin
    !     pk         Real(r8) LCL pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     ftheq      Real(r8) equivalent potential temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftheq
    REAL(r8),INTENT(in):: t,pk
    INTEGER jx,jy
    REAL(r8) xj,yj,dxj,dyj
    REAL(r8) ft11,ft12,ft13,ft21,ft22,ft23,ft31,ft32,ft33
    REAL(r8) ftx1,ftx2,ftx3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(c1xthe+c2xthe*t,REAL(nxthe,r8))
    yj=MIN(c1ythe+c2ythe*pk,REAL(nythe,r8))
    IF(xj.GE.1.0_r8.AND.yj.GE.1.0_r8) THEN
       jx=MIN(MAX(NINT(xj),2),nxthe-1)
       jy=MIN(MAX(NINT(yj),2),nythe-1)
       dxj=xj-jx
       dyj=yj-jy
       ft11=tbthe(jx-1,jy-1)
       ft12=tbthe(jx-1,jy)
       ft13=tbthe(jx-1,jy+1)
       ft21=tbthe(jx,jy-1)
       ft22=tbthe(jx,jy)
       ft23=tbthe(jx,jy+1)
       ft31=tbthe(jx+1,jy-1)
       ft32=tbthe(jx+1,jy)
       ft33=tbthe(jx+1,jy+1)
       ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21
       ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22
       ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23
       ftheq=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2
    ELSE
       ftheq=0.0_r8
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftheq
  !-------------------------------------------------------------------------------
  ! elemental function fthex(t,pk)
  FUNCTION fthex(t,pk)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fthex        Compute equivalent potential temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute equivalent potential temperature at the LCL
    !   from temperature and pressure over 1e5 Pa to the kappa power.
    !   Equivalent potential temperature is constant for a saturated parcel
    !   rising adiabatically up a moist adiabat when the heat and mass
    !   of the condensed water are neglected.  Ice is also neglected.
    !   The formula for equivalent potential temperature (Holton) is
    !       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))
    !   where t is the temperature, pv is the saturated vapor pressure,
    !   pd is the dry pressure p-pv, el is the temperature dependent
    !   latent heat of condensation hvap+dldt*(t-ttp), and other values
    !   are physical constants defined in parameter statements in the code.
    !   Zero is returned if the input values make saturation impossible.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   the=fthex(t,pk)
    !
    !   Input argument list:
    !     t          Real(r8) LCL temperature in Kelvin
    !     pk         Real(r8) LCL pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     fthex      Real(r8) equivalent potential temperature in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fthex
    REAL(r8),INTENT(in):: t,pk
    REAL(r8) p,tr,pv,pd,el,expo
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    p=pk**con_cpor
    tr=con_ttp/t
    pv=psatb*(tr**con_xpona)*EXP(con_xponb*(1.0_r8-tr))
    pd=p-pv
    IF(pd.GT.pv) THEN
       el=con_hvap+con_dldt*(t-con_ttp)
       expo=el*con_eps*pv/(con_cp*t*pd)
       fthex=t*pd**(-con_rocp)*EXP(expo)
    ELSE
       fthex=0.0_r8
    ENDIF
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fthex
  !-------------------------------------------------------------------------------
  SUBROUTINE gtma
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gtma         Compute moist adiabat tables
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute temperature and specific humidity tables
    !   as a function of equivalent potential temperature and
    !   pressure over 1e5 Pa to the kappa power for subprogram stma.
    !   Exact parcel temperatures are calculated in subprogram stmaxg.
    !   The current implementation computes a table with a first dimension
    !   of 151 for equivalent potential temperatures ranging from 200 to 500
    !   Kelvin and a second dimension of 121 for pressure over 1e5 Pa
    !   to the kappa power ranging from 0.01**rocp to 1.10**rocp.
    !
    ! Program History Log:
    !   91-05-07  Iredell
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gtma
    !
    ! Subprograms called:
    !   (stmaxg)   inlinable subprogram to compute parcel temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx,jy
    REAL(r8) xmin,xmax,ymin,ymax,xinc,yinc,x,y,pk,the,t,q,tg
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=200.0_r8
    xmax=500.0_r8
    ymin=0.01_r8**con_rocp
    ymax=1.10_r8**con_rocp
    xinc=(xmax-xmin)/(nxma-1)
    c1xma=1.0_r8-xmin/xinc
    c2xma=1.0_r8/xinc
    yinc=(ymax-ymin)/(nyma-1)
    c1yma=1.0_r8-ymin/yinc
    c2yma=1.0_r8/yinc
    DO jy=1,nyma
       y=ymin+(jy-1)*yinc
       pk=y
       tg=xmin*y
       DO jx=1,nxma
          x=xmin+(jx-1)*xinc
          the=x
          CALL stmaxg(tg,the,pk,t,q)
          tbtma(jx,jy)=t
          tbqma(jx,jy)=q
          tg=t
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gtma
  !-------------------------------------------------------------------------------
  ! elemental subroutine stma(the,pk,tma,qma)
  SUBROUTINE stma(the,pk,tma,qma)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: stma         Compute moist adiabat temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute temperature and specific humidity of a parcel
    !   lifted up a moist adiabat from equivalent potential temperature
    !   at the LCL and pressure over 1e5 Pa to the kappa power.
    !   Bilinear interpolations are done between values in a lookup table
    !   computed in gtma. See documentation for stmaxg for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.01 Kelvin
    !   and 5.e-6 kg/kg for temperature and humidity, respectively.
    !   On the Cray, stma is about 35 times faster than exact calculation.
    !   This subprogram should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             expand table
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call stma(the,pk,tma,qma)
    !
    !   Input argument list:
    !     the        Real(r8) equivalent potential temperature in Kelvin
    !     pk         Real(r8) pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     tma        Real(r8) parcel temperature in Kelvin
    !     qma        Real(r8) parcel specific humidity in kg/kg
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8),INTENT(in):: the,pk
    REAL(r8),INTENT(out):: tma,qma
    INTEGER jx,jy
    REAL(r8) xj,yj,ftx1,ftx2,qx1,qx2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xma+c2xma*the,1.0_r8),REAL(nxma,r8))
    yj=MIN(MAX(c1yma+c2yma*pk,1.0_r8),REAL(nyma,r8))
    jx=MIN(xj,nxma-1.0_r8)
    jy=MIN(yj,nyma-1.0_r8)
    ftx1=tbtma(jx,jy)+(xj-jx)*(tbtma(jx+1,jy)-tbtma(jx,jy))
    ftx2=tbtma(jx,jy+1)+(xj-jx)*(tbtma(jx+1,jy+1)-tbtma(jx,jy+1))
    tma=ftx1+(yj-jy)*(ftx2-ftx1)
    qx1=tbqma(jx,jy)+(xj-jx)*(tbqma(jx+1,jy)-tbqma(jx,jy))
    qx2=tbqma(jx,jy+1)+(xj-jx)*(tbqma(jx+1,jy+1)-tbqma(jx,jy+1))
    qma=qx1+(yj-jy)*(qx2-qx1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE stma
  !-------------------------------------------------------------------------------
  ! elemental subroutine stmaq(the,pk,tma,qma)
  SUBROUTINE stmaq(the,pk,tma,qma)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: stmaq        Compute moist adiabat temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute temperature and specific humidity of a parcel
    !   lifted up a moist adiabat from equivalent potential temperature
    !   at the LCL and pressure over 1e5 Pa to the kappa power.
    !   Biquadratic interpolations are done between values in a lookup table
    !   computed in gtma. See documentation for stmaxg for details.
    !   Input values outside table range are reset to table extrema.
    !   the interpolation accuracy is better than 0.0005 Kelvin
    !   and 1.e-7 kg/kg for temperature and humidity, respectively.
    !   On the Cray, stmaq is about 25 times faster than exact calculation.
    !   This subprogram should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             quadratic interpolation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call stmaq(the,pk,tma,qma)
    !
    !   Input argument list:
    !     the        Real(r8) equivalent potential temperature in Kelvin
    !     pk         Real(r8) pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     tmaq       Real(r8) parcel temperature in Kelvin
    !     qma        Real(r8) parcel specific humidity in kg/kg
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8),INTENT(in):: the,pk
    REAL(r8),INTENT(out):: tma,qma
    INTEGER jx,jy
    REAL(r8) xj,yj,dxj,dyj
    REAL(r8) ft11,ft12,ft13,ft21,ft22,ft23,ft31,ft32,ft33
    REAL(r8) ftx1,ftx2,ftx3
    REAL(r8) q11,q12,q13,q21,q22,q23,q31,q32,q33,qx1,qx2,qx3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xma+c2xma*the,1.0_r8),REAL(nxma,r8))
    yj=MIN(MAX(c1yma+c2yma*pk,1.0_r8),REAL(nyma,r8))
    jx=MIN(MAX(NINT(xj),2),nxma-1)
    jy=MIN(MAX(NINT(yj),2),nyma-1)
    dxj=xj-jx
    dyj=yj-jy
    ft11=tbtma(jx-1,jy-1)
    ft12=tbtma(jx-1,jy)
    ft13=tbtma(jx-1,jy+1)
    ft21=tbtma(jx,jy-1)
    ft22=tbtma(jx,jy)
    ft23=tbtma(jx,jy+1)
    ft31=tbtma(jx+1,jy-1)
    ft32=tbtma(jx+1,jy)
    ft33=tbtma(jx+1,jy+1)
    ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21
    ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22
    ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23
    tma=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2
    q11=tbqma(jx-1,jy-1)
    q12=tbqma(jx-1,jy)
    q13=tbqma(jx-1,jy+1)
    q21=tbqma(jx,jy-1)
    q22=tbqma(jx,jy)
    q23=tbqma(jx,jy+1)
    q31=tbqma(jx+1,jy-1)
    q32=tbqma(jx+1,jy)
    q33=tbqma(jx+1,jy+1)
    qx1=(((q31+q11)/2-q21)*dxj+(q31-q11)/2)*dxj+q21
    qx2=(((q32+q12)/2-q22)*dxj+(q32-q12)/2)*dxj+q22
    qx3=(((q33+q13)/2-q23)*dxj+(q33-q13)/2)*dxj+q23
    qma=(((qx3+qx1)/2-qx2)*dyj+(qx3-qx1)/2)*dyj+qx2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE stmaq
  !-------------------------------------------------------------------------------
  ! elemental subroutine stmax(the,pk,tma,qma)
  SUBROUTINE stmax(the,pk,tma,qma)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: stmax        Compute moist adiabat temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Exactly compute temperature and humidity of a parcel
    !   lifted up a moist adiabat from equivalent potential temperature
    !   at the LCL and pressure over 1e5 Pa to the kappa power.
    !   An approximate parcel temperature for subprogram stmaxg
    !   is obtained using stma so gtma must be already called.
    !   See documentation for stmaxg for details.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call stmax(the,pk,tma,qma)
    !
    !   Input argument list:
    !     the        Real(r8) equivalent potential temperature in Kelvin
    !     pk         Real(r8) pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     tma        Real(r8) parcel temperature in Kelvin
    !     qma        Real(r8) parcel specific humidity in kg/kg
    !
    ! Subprograms called:
    !   (stma)     inlinable subprogram to compute parcel temperature
    !   (stmaxg)   inlinable subprogram to compute parcel temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8),INTENT(in):: the,pk
    REAL(r8),INTENT(out):: tma,qma
    REAL(r8) tg,qg
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    CALL stma(the,pk,tg,qg)
    CALL stmaxg(tg,the,pk,tma,qma)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE stmax
  !-------------------------------------------------------------------------------
  ! elemental subroutine stmaxg(tg,the,pk,tma,qma)
  SUBROUTINE stmaxg(tg,the,pk,tma,qma)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: stmaxg       Compute moist adiabat temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: exactly compute temperature and humidity of a parcel
    !   lifted up a moist adiabat from equivalent potential temperature
    !   at the LCL and pressure over 1e5 Pa to the kappa power.
    !   A guess parcel temperature must be provided.
    !   Equivalent potential temperature is constant for a saturated parcel
    !   rising adiabatically up a moist adiabat when the heat and mass
    !   of the condensed water are neglected.  Ice is also neglected.
    !   The formula for equivalent potential temperature (Holton) is
    !       the=t*(pd**(-rocp))*exp(el*eps*pv/(cp*t*pd))
    !   where t is the temperature, pv is the saturated vapor pressure,
    !   pd is the dry pressure p-pv, el is the temperature dependent
    !   latent heat of condensation hvap+dldt*(t-ttp), and other values
    !   are physical constants defined in parameter statements in the code.
    !   The formula is inverted by iterating Newtonian approximations
    !   for each the and p until t is found to within 1.e-4 Kelvin.
    !   The specific humidity is then computed from pv and pd.
    !   This subprogram can be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             exact computation
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call stmaxg(tg,the,pk,tma,qma)
    !
    !   Input argument list:
    !     tg         Real(r8) guess parcel temperature in Kelvin
    !     the        Real(r8) equivalent potential temperature in Kelvin
    !     pk         Real(r8) pressure over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     tma        Real(r8) parcel temperature in Kelvin
    !     qma        Real(r8) parcel specific humidity in kg/kg
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8),INTENT(in):: tg,the,pk
    REAL(r8),INTENT(out):: tma,qma
    REAL(r8),PARAMETER:: terrm=1.e-4_r8
    REAL(r8) t,p,tr,pv,pd,el,expo,thet,dthet,terr
    INTEGER i
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    t=tg
    p=pk**con_cpor
    DO i=1,100
       tr=con_ttp/t
       pv=psatb*(tr**con_xpona)*EXP(con_xponb*(1.0_r8-tr))
       pd=p-pv
       el=con_hvap+con_dldt*(t-con_ttp)
       expo=el*con_eps*pv/(con_cp*t*pd)
       thet=t*pd**(-con_rocp)*EXP(expo)
       dthet=thet/t*(1.0_r8+expo*(con_dldt*t/el+el*p/(con_rv*t*pd)))
       terr=(thet-the)/dthet
       t=t-terr
       IF(ABS(terr).LE.terrm) EXIT
    ENDDO
    tma=t
    tr=con_ttp/t
    pv=psatb*(tr**con_xpona)*EXP(con_xponb*(1.0_r8-tr))
    pd=p-pv
    qma=con_eps*pv/(pd+con_eps*pv)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE stmaxg
  !-------------------------------------------------------------------------------
  SUBROUTINE gpkap
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: gpkap        Compute coefficients for p**kappa
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: Computes pressure to the kappa table as a function of pressure
    !   for the table lookup function fpkap.
    !   Exact pressure to the kappa values are calculated in subprogram fpkapx.
    !   The current implementation computes a table with a length
    !   of 5501 for pressures ranging up to 110000 Pascals.
    !
    ! Program History Log:
    !   94-12-30  Iredell
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:  call gpkap
    !
    ! Subprograms called:
    !   fpkapx     function to compute exact pressure to the kappa
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,x,p
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=0.0_r8
    xmax=110000.0_r8
    xinc=(xmax-xmin)/(nxpkap-1)
    c1xpkap=1.0_r8-xmin/xinc
    c2xpkap=1.0_r8/xinc
    DO jx=1,nxpkap
       x=xmin+(jx-1)*xinc
       p=x
       tbpkap(jx)=fpkapx(p)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gpkap
  !-------------------------------------------------------------------------------
  ! elemental function fpkap(p)
  FUNCTION fpkap(p)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpkap        raise pressure to the kappa power.
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Raise pressure over 1e5 Pa to the kappa power.
    !   A linear interpolation is done between values in a lookup table
    !   computed in gpkap. See documentation for fpkapx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy ranges from 9 decimal places
    !   at 100000 Pascals to 5 decimal places at 1000 Pascals.
    !   On the Cray, fpkap is over 5 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             standardized kappa,
    !                                 increased range and accuracy
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:   pkap=fpkap(p)
    !
    !   Input argument list:
    !     p          Real(r8) pressure in Pascals
    !
    !   Output argument list:
    !     fpkap      Real(r8) p over 1e5 Pa to the kappa power
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fpkap
    REAL(r8),INTENT(in):: p
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpkap+c2xpkap*p,1.0_r8),REAL(nxpkap,r8))
    jx=MIN(xj,nxpkap-1.0_r8)
    fpkap=tbpkap(jx)+(xj-jx)*(tbpkap(jx+1)-tbpkap(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpkap
  !-------------------------------------------------------------------------------
  ! elemental function fpkapq(p)
  FUNCTION fpkapq(p)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: fpkapq       raise pressure to the kappa power.
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Raise pressure over 1e5 Pa to the kappa power.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in gpkap. see documentation for fpkapx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy ranges from 12 decimal places
    !   at 100000 Pascals to 7 decimal places at 1000 Pascals.
    !   On the Cray, fpkap is over 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             standardized kappa,
    !                                 increased range and accuracy
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:   pkap=fpkapq(p)
    !
    !   Input argument list:
    !     p          Real(r8) pressure in Pascals
    !
    !   Output argument list:
    !     fpkapq     Real(r8) p over 1e5 Pa to the kappa power
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fpkapq
    REAL(r8),INTENT(in):: p
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xpkap+c2xpkap*p,1.0_r8),REAL(nxpkap,r8))
    jx=MIN(MAX(NINT(xj),2),nxpkap-1)
    dxj=xj-jx
    fj1=tbpkap(jx-1)
    fj2=tbpkap(jx)
    fj3=tbpkap(jx+1)
    fpkapq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpkapq
  !-------------------------------------------------------------------------------
  FUNCTION fpkapo(p)
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: fpkapo       raise surface pressure to the kappa power.
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: Raise surface pressure over 1e5 Pa to the kappa power
    !   using a rational weighted chebyshev approximation.
    !   The numerator is of order 2 and the denominator is of order 4.
    !   The pressure range is 40000-110000 Pa and kappa is defined in fpkapx.
    !   The accuracy of this approximation is almost 8 decimal places.
    !   On the Cray, fpkap is over 10 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             standardized kappa,
    !                                 increased range and accuracy
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  pkap=fpkapo(p)
    !
    !   Input argument list:
    !     p          Real(r8) surface pressure in Pascals
    !                p should be in the range 40000 to 110000
    !
    !   Output argument list:
    !     fpkapo     Real(r8) p over 1e5 Pa to the kappa power
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fpkapo
    REAL(r8),INTENT(in):: p
    INTEGER,PARAMETER:: nnpk=2,ndpk=4
    REAL(r8):: cnpk(0:nnpk)=(/3.13198449e-1_r8,5.78544829e-2_r8,&
         8.35491871e-4_r8/)
    REAL(r8):: cdpk(0:ndpk)=(/1._r8,8.15968401e-2_r8,5.72839518e-4_r8,&
         -4.86959812e-7_r8,5.24459889e-10_r8/)
    INTEGER  :: n
    REAL(r8) :: pkpa,fnpk,fdpk
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pkpa=p*1.e-3_r8
    fnpk=cnpk(nnpk)
    DO n=nnpk-1,0,-1
       fnpk=pkpa*fnpk+cnpk(n)
    ENDDO
    fdpk=cdpk(ndpk)
    DO n=ndpk-1,0,-1
       fdpk=pkpa*fdpk+cdpk(n)
    ENDDO
    fpkapo=fnpk/fdpk
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpkapo
  !-------------------------------------------------------------------------------
  ! elemental function fpkapx(p)
  FUNCTION fpkapx(p)
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: fpkapx       raise pressure to the kappa power.
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: raise pressure over 1e5 Pa to the kappa power.
    !   Kappa is equal to rd/cp where rd and cp are physical constants.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   94-12-30  Iredell             made into inlinable function
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  pkap=fpkapx(p)
    !
    !   Input argument list:
    !     p          Real(r8) pressure in Pascals
    !
    !   Output argument list:
    !     fpkapx     Real(r8) p over 1e5 Pa to the kappa power
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) fpkapx
    REAL(r8),INTENT(in):: p
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    fpkapx=(p/1.e5_r8)**con_rocp
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION fpkapx
  !-------------------------------------------------------------------------------
  SUBROUTINE grkap
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: grkap        Compute coefficients for p**(1/kappa)
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: Computes pressure to the 1/kappa table as a function of pressure
    !   for the table lookup function frkap.
    !   Exact pressure to the 1/kappa values are calculated in subprogram frkapx.
    !   The current implementation computes a table with a length
    !   of 5501 for pressures ranging up to 110000 Pascals.
    !
    ! Program History Log:
    !   94-12-30  Iredell
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:  call grkap
    !
    ! Subprograms called:
    !   frkapx     function to compute exact pressure to the 1/kappa
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx
    REAL(r8) xmin,xmax,xinc,x,p
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=0.0_r8
    xmax=fpkapx(110000.0_r8)
    xinc=(xmax-xmin)/(nxrkap-1)
    c1xrkap=1.0_r8-xmin/xinc
    c2xrkap=1.0_r8/xinc
    DO jx=1,nxrkap
       x=xmin+(jx-1)*xinc
       p=x
       tbrkap(jx)=frkapx(p)
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE grkap
  !-------------------------------------------------------------------------------
  ! elemental function frkap(pkap)
  FUNCTION frkap(pkap)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: frkap        raise pressure to the 1/kappa power.
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Raise pressure over 1e5 Pa to the 1/kappa power.
    !   A linear interpolation is done between values in a lookup table
    !   computed in grkap. See documentation for frkapx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 7 decimal places.
    !   On the IBM, fpkap is about 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             standardized kappa,
    !                                 increased range and accuracy
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:   p=frkap(pkap)
    !
    !   Input argument list:
    !     pkap       Real(r8) p over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     frkap      Real(r8) pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) frkap
    REAL(r8),INTENT(in):: pkap
    INTEGER jx
    REAL(r8) xj
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xrkap+c2xrkap*pkap,1.0_r8),REAL(nxrkap,r8))
    jx=MIN(xj,nxrkap-1.0_r8)
    frkap=tbrkap(jx)+(xj-jx)*(tbrkap(jx+1)-tbrkap(jx))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION frkap
  !-------------------------------------------------------------------------------
  ! elemental function frkapq(pkap)
  FUNCTION frkapq(pkap)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: frkapq       raise pressure to the 1/kappa power.
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Raise pressure over 1e5 Pa to the 1/kappa power.
    !   A quadratic interpolation is done between values in a lookup table
    !   computed in grkap. see documentation for frkapx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 11 decimal places.
    !   On the IBM, fpkap is almost 4 times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    !   94-12-30  Iredell             standardized kappa,
    !                                 increased range and accuracy
    ! 1999-03-01  Iredell             f90 module
    ! 1999-03-24  Iredell             table lookup
    !
    ! Usage:   p=frkapq(pkap)
    !
    !   Input argument list:
    !     pkap       Real(r8) p over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     frkapq     Real(r8) pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) frkapq
    REAL(r8),INTENT(in):: pkap
    INTEGER jx
    REAL(r8) xj,dxj,fj1,fj2,fj3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xrkap+c2xrkap*pkap,1._r8),REAL(nxrkap,r8))
    jx=MIN(MAX(NINT(xj),2),nxrkap-1)
    dxj=xj-jx
    fj1=tbrkap(jx-1)
    fj2=tbrkap(jx)
    fj3=tbrkap(jx+1)
    frkapq=(((fj3+fj1)/2-fj2)*dxj+(fj3-fj1)/2)*dxj+fj2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION frkapq
  !-------------------------------------------------------------------------------
  ! elemental function frkapx(pkap)
  FUNCTION frkapx(pkap)
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: frkapx       raise pressure to the 1/kappa power.
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: raise pressure over 1e5 Pa to the 1/kappa power.
    !   Kappa is equal to rd/cp where rd and cp are physical constants.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   94-12-30  Iredell             made into inlinable function
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  p=frkapx(pkap)
    !
    !   Input argument list:
    !     pkap       Real(r8) p over 1e5 Pa to the kappa power
    !
    !   Output argument list:
    !     frkapx     Real(r8) pressure in Pascals
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) frkapx
    REAL(r8),INTENT(in):: pkap
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    frkapx=pkap**(1/con_rocp)*1.e5_r8
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION frkapx
  !-------------------------------------------------------------------------------
  SUBROUTINE gtlcl
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: gtlcl        Compute equivalent potential temperature table
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute lifting condensation level temperature table
    !   as a function of temperature and dewpoint depression for function ftlcl.
    !   Lifting condensation level temperature is calculated in subprogram ftlclx
    !   The current implementation computes a table with a first dimension
    !   of 151 for temperatures ranging from 180.0 to 330.0 Kelvin
    !   and a second dimension of 61 for dewpoint depression ranging from
    !   0 to 60 Kelvin.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  call gtlcl
    !
    ! Subprograms called:
    !   (ftlclx)    inlinable function to compute LCL temperature
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    INTEGER jx,jy
    REAL(r8) xmin,xmax,ymin,ymax,xinc,yinc,x,y,tdpd,t
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xmin=180._r8
    xmax=330._r8
    ymin=0._r8
    ymax=60._r8
    xinc=(xmax-xmin)/(nxtlcl-1)
    c1xtlcl=1.0_r8-xmin/xinc
    c2xtlcl=1.0_r8/xinc
    yinc=(ymax-ymin)/(nytlcl-1)
    c1ytlcl=1.0_r8-ymin/yinc
    c2ytlcl=1.0_r8/yinc
    DO jy=1,nytlcl
       y=ymin+(jy-1)*yinc
       tdpd=y
       DO jx=1,nxtlcl
          x=xmin+(jx-1)*xinc
          t=x
          tbtlcl(jx,jy)=ftlclx(t,tdpd)
       ENDDO
    ENDDO
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END SUBROUTINE gtlcl
  !-------------------------------------------------------------------------------
  ! elemental function ftlcl(t,tdpd)
  FUNCTION ftlcl(t,tdpd)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftlcl        Compute LCL temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute temperature at the lifting condensation level
    !   from temperature and dewpoint depression.
    !   A bilinear interpolation is done between values in a lookup table
    !   computed in gtlcl. See documentation for ftlclx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.0005 Kelvin.
    !   On the Cray, ftlcl is ? times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tlcl=ftlcl(t,tdpd)
    !
    !   Input argument list:
    !     t          Real(r8) LCL temperature in Kelvin
    !     tdpd       Real(r8) dewpoint depression in Kelvin
    !
    !   Output argument list:
    !     ftlcl      Real(r8) temperature at the LCL in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftlcl
    REAL(r8),INTENT(in):: t,tdpd
    INTEGER jx,jy
    REAL(r8) xj,yj,ftx1,ftx2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtlcl+c2xtlcl*t,1.0_r8),REAL(nxtlcl,r8))
    yj=MIN(MAX(c1ytlcl+c2ytlcl*tdpd,1.0_r8),REAL(nytlcl,r8))
    jx=MIN(xj,nxtlcl-1.0_r8)
    jy=MIN(yj,nytlcl-1.0_r8)
    ftx1=tbtlcl(jx,jy)+(xj-jx)*(tbtlcl(jx+1,jy)-tbtlcl(jx,jy))
    ftx2=tbtlcl(jx,jy+1)+(xj-jx)*(tbtlcl(jx+1,jy+1)-tbtlcl(jx,jy+1))
    ftlcl=ftx1+(yj-jy)*(ftx2-ftx1)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftlcl
  !-------------------------------------------------------------------------------
  ! elemental function ftlclq(t,tdpd)
  FUNCTION ftlclq(t,tdpd)
    !$$$     Subprogram Documentation Block
    !
    ! Subprogram: ftlclq       Compute LCL temperature
    !   Author: N Phillips            w/NMC2X2   Date: 30 dec 82
    !
    ! Abstract: Compute temperature at the lifting condensation level
    !   from temperature and dewpoint depression.
    !   A biquadratic interpolation is done between values in a lookup table
    !   computed in gtlcl. see documentation for ftlclx for details.
    !   Input values outside table range are reset to table extrema.
    !   The interpolation accuracy is better than 0.000003 Kelvin.
    !   On the Cray, ftlclq is ? times faster than exact calculation.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:   tlcl=ftlclq(t,tdpd)
    !
    !   Input argument list:
    !     t          Real(r8) LCL temperature in Kelvin
    !     tdpd       Real(r8) dewpoint depression in Kelvin
    !
    !   Output argument list:
    !     ftlcl      Real(r8) temperature at the LCL in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftlclq
    REAL(r8),INTENT(in):: t,tdpd
    INTEGER jx,jy
    REAL(r8) xj,yj,dxj,dyj
    REAL(r8) ft11,ft12,ft13,ft21,ft22,ft23,ft31,ft32,ft33
    REAL(r8) ftx1,ftx2,ftx3
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    xj=MIN(MAX(c1xtlcl+c2xtlcl*t,1.0_r8),REAL(nxtlcl,r8))
    yj=MIN(MAX(c1ytlcl+c2ytlcl*tdpd,1.0_r8),REAL(nytlcl,r8))
    jx=MIN(MAX(NINT(xj),2),nxtlcl-1)
    jy=MIN(MAX(NINT(yj),2),nytlcl-1)
    dxj=xj-jx
    dyj=yj-jy
    ft11=tbtlcl(jx-1,jy-1)
    ft12=tbtlcl(jx-1,jy)
    ft13=tbtlcl(jx-1,jy+1)
    ft21=tbtlcl(jx,jy-1)
    ft22=tbtlcl(jx,jy)
    ft23=tbtlcl(jx,jy+1)
    ft31=tbtlcl(jx+1,jy-1)
    ft32=tbtlcl(jx+1,jy)
    ft33=tbtlcl(jx+1,jy+1)
    ftx1=(((ft31+ft11)/2-ft21)*dxj+(ft31-ft11)/2)*dxj+ft21
    ftx2=(((ft32+ft12)/2-ft22)*dxj+(ft32-ft12)/2)*dxj+ft22
    ftx3=(((ft33+ft13)/2-ft23)*dxj+(ft33-ft13)/2)*dxj+ft23
    ftlclq=(((ftx3+ftx1)/2-ftx2)*dyj+(ftx3-ftx1)/2)*dyj+ftx2
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftlclq
  !-------------------------------------------------------------------------------
  FUNCTION ftlclo(t,tdpd)
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: ftlclo       Compute LCL temperature.
    !   Author: Phillips         org: w/NMC2X2   Date: 29 dec 82
    !
    ! Abstract: Compute temperature at the lifting condensation level
    !   from temperature and dewpoint depression.  the formula used is
    !   a polynomial taken from Phillips mstadb routine which empirically
    !   approximates the original exact implicit relationship.
    !   (This kind of approximation is customary (inman, 1969), but
    !   the original source for this particular one is not yet known. -MI)
    !   Its accuracy is about 0.03 Kelvin for a dewpoint depression of 30.
    !   This function should be expanded inline in the calling routine.
    !
    ! Program History Log:
    !   91-05-07  Iredell             made into inlinable function
    ! 1999-03-01  Iredell             f90 module
    !
    ! Usage:  tlcl=ftlclo(t,tdpd)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !     tdpd       Real(r8) dewpoint depression in Kelvin
    !
    !   Output argument list:
    !     ftlclo     Real(r8) temperature at the LCL in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftlclo
    REAL(r8),INTENT(in):: t,tdpd
    REAL(r8),PARAMETER:: clcl1= 0.954442e+0_r8,clcl2= 0.967772e-3_r8,&
         clcl3=-0.710321e-3_r8,clcl4=-0.270742e-5_r8
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    ftlclo=t-tdpd*(clcl1+clcl2*t+tdpd*(clcl3+clcl4*t))
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftlclo
  !-------------------------------------------------------------------------------
  ! elemental function ftlclx(t,tdpd)
  FUNCTION ftlclx(t,tdpd)
    !$$$   Subprogram  documentation  block
    !
    ! Subprogram: ftlclx       Compute LCL temperature.
    !   Author: Iredell          org: w/NMC2X2   Date: 25 March 1999
    !
    ! Abstract: Compute temperature at the lifting condensation level
    !   from temperature and dewpoint depression.  A parcel lifted
    !   adiabatically becomes saturated at the lifting condensation level.
    !   The water model assumes a perfect gas, constant specific heats
    !   for gas and liquid, and neglects the volume of the liquid.
    !   The model does account for the variation of the latent heat
    !   of condensation with temperature.  The ice option is not included.
    !   The Clausius-Clapeyron equation is integrated from the triple point
    !   to get the formulas
    !       pvlcl=con_psat*(trlcl**xa)*exp(xb*(1.-trlcl))
    !       pvdew=con_psat*(trdew**xa)*exp(xb*(1.-trdew))
    !   where pvlcl is the saturated parcel vapor pressure at the LCL,
    !   pvdew is the unsaturated parcel vapor pressure initially,
    !   trlcl is ttp/tlcl and trdew is ttp/tdew.  The adiabatic lifting
    !   of the parcel is represented by the following formula
    !       pvdew=pvlcl*(t/tlcl)**(1/kappa)
    !   This formula is inverted by iterating Newtonian approximations
    !   until tlcl is found to within 1.e-6 Kelvin.  Note that the minimum
    !   returned temperature is 180 Kelvin.
    !
    ! Program History Log:
    ! 1999-03-25  Iredell
    !
    ! Usage:  tlcl=ftlclx(t,tdpd)
    !
    !   Input argument list:
    !     t          Real(r8) temperature in Kelvin
    !     tdpd       Real(r8) dewpoint depression in Kelvin
    !
    !   Output argument list:
    !     ftlclx     Real(r8) temperature at the LCL in Kelvin
    !
    ! Attributes:
    !   Language: Fortran 90.
    !
    !$$$
    IMPLICIT NONE
    REAL(r8) ftlclx
    REAL(r8),INTENT(in):: t,tdpd
    REAL(r8),PARAMETER:: terrm=1.e-4_r8,tlmin=180.0_r8,tlminx=tlmin-5.0_r8
    REAL(r8) tr,pvdew,tlcl,ta,pvlcl,el,dpvlcl,terr
    INTEGER i
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tr=con_ttp/(t-tdpd)
    pvdew=con_psat*(tr**con_xpona)*EXP(con_xponb*(1.0_r8-tr))
    tlcl=t-tdpd
    DO i=1,100
       tr=con_ttp/tlcl
       ta=t/tlcl
       pvlcl=con_psat*(tr**con_xpona)*EXP(con_xponb*(1.0_r8-tr))*ta**(1/con_rocp)
       el=con_hvap+con_dldt*(tlcl-con_ttp)
       dpvlcl=(el/(con_rv*t**2)+1/(con_rocp*tlcl))*pvlcl
       terr=(pvlcl-pvdew)/dpvlcl
       tlcl=tlcl-terr
       IF(ABS(terr).LE.terrm.OR.tlcl.LT.tlminx) EXIT
    ENDDO
    ftlclx=MAX(tlcl,tlmin)
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  END FUNCTION ftlclx
  !-------------------------------------------------------------------------------
  !
  !  MODULE module_calc_cape  
  !
  SUBROUTINE calc_cape( &
       nCols                , &!INTEGER      , INTENT(IN   ) :: nCols
       kMax                 , &!INTEGER      , INTENT(IN   ) :: kMax
       si                   , &!REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
       sl                   , &!REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
       PSFC                 , &!REAL(KIND=r8), INTENT(IN   ) :: PSFC (nCols)   ! psfc is pressure in hPa
       HGT                  , &!REAL(KIND=r8), INTENT(IN   ) :: HGT  (nCols)   ! topography m
       TK                   , &!REAL(KIND=r8), INTENT(IN   ) :: TK    (nCols,kMax)     ! TK is temp in K, T is theta-300
       QV                   , &!REAL(KIND=r8), INTENT(IN   ) :: QV    (nCols,kMax)
       SCRa                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRa  (nCols,kMax)
       SCRb                 , &!REAL(KIND=r8), INTENT(OUT  ) :: SCRb  (nCols,kMax)
       i3dflag                )!INTEGER      , INTENT(IN   ) :: i3dflag
    !   If i3dflag=1, this routine calculates CAPE and CIN (in m**2/s**2,
    !   or J/kg) for every grid point in the entire 3D domain (treating
    !   each grid point as a parcel).  If i3dflag=0, then it
    !   calculates CAPE and CIN only for the parcel with max theta-e in
    !   the column, (i.e. something akin to Colman's MCAPE).  By "parcel",
    !   we mean a 500-m deep parcel, with actual temperature and moisture
    !   averaged over that depth.
    !
    !   In the case of i3dflag=0,
    !   MCAPE, MCIN, LCL and LFC (2D fields are calculated)


    !  USE constants_module
    !  USE module_model_basics

    IMPLICIT NONE
    REAL(KIND=r8), PARAMETER :: EPS = 0.622_r8
    REAL(KIND=r8), PARAMETER :: G = 9.81_r8
    REAL(KIND=r8), PARAMETER :: ESLCON1 = 17.67_r8
    REAL(KIND=r8), PARAMETER :: ESLCON2 = 29.65_r8
    REAL(KIND=r8), PARAMETER :: Rm = .608_r8 
    REAL(KIND=r8), PARAMETER :: CPMD = 0.887_r8
    REAL(KIND=r8), PARAMETER :: EZERO = 6.112_r8
    REAL(KIND=r8), PARAMETER :: Rd = 287.04_r8
    REAL(KIND=r8), PARAMETER :: Cp = 7.0_r8*Rd/2.0_r8
    REAL(KIND=r8), PARAMETER :: GAMMA_RIP = Rd/Cp 
    REAL(KIND=r8), PARAMETER :: GAMMAMD = Rm-CPMD
    REAL(KIND=r8), PARAMETER :: CELKEL = 273.15_r8
    REAL(KIND=r8), PARAMETER :: THTECON1 = 3376.0_r8
    REAL(KIND=r8), PARAMETER :: THTECON2 = 2.54_r8
    REAL(KIND=r8), PARAMETER :: THTECON3 = 0.81_r8
    REAL(KIND=r8), PARAMETER :: TLCLC1 = 2840.0_r8
    REAL(KIND=r8), PARAMETER :: TLCLC2 = 3.5_r8
    REAL(KIND=r8), PARAMETER :: TLCLC3 = 4.805_r8
    REAL(KIND=r8), PARAMETER :: TLCLC4 = 55.0_r8

    !Arguments
    INTEGER, INTENT(IN   ) :: nCols
    INTEGER, INTENT(IN   ) :: kMax
    REAL(KIND=r8), INTENT(IN   ) :: si   (kMax+1)
    REAL(KIND=r8), INTENT(IN   ) :: sl   (kMax)
    REAL(KIND=r8), INTENT(IN   ) :: PSFC (nCols)   ! psfc is pressure in hPa
    REAL(KIND=r8), INTENT(IN   ) :: HGT  (nCols)   ! topography m

    REAL(KIND=r8), INTENT(IN   ) :: TK    (nCols,kMax)     ! TK is temp in K, T is theta-300
    REAL(KIND=r8), INTENT(IN   ) :: QV    (nCols,kMax)

    REAL(KIND=r8), INTENT(OUT  ) :: SCRa  (nCols,kMax)
    REAL(KIND=r8), INTENT(OUT  ) :: SCRb  (nCols,kMax)
    INTEGER            , INTENT(IN    ) :: i3dflag

    ! Local variables
    CHARACTER(len=128) :: cname
    CHARACTER(len=128) :: cdesc
    CHARACTER(len=128) :: cunits

    INTEGER :: i, k, kk
    INTEGER :: kpar, kpar1, kpar2, kLevel, klev, kel
    INTEGER :: ilcl, klcl, klfc
    REAL(KIND=r8), DIMENSION(nCols):: ter
    REAL(KIND=r8), DIMENSION(nCols,kMax)  :: prs, tmk, ght, qvp,delz
    REAL(KIND=r8), DIMENSION(nCols,kMax)  :: cape, cin
    REAL(KIND=r8), DIMENSION(nCols,kMax+1)  :: prsf
    REAL(KIND=r8):: ethpari, zlcl, tvenv
    REAL(KIND=r8):: p1, p2, pp1, pp2, pup, pdn,tv1,tv2
    REAL(KIND=r8):: totprs, totqvp, totthe
    REAL(KIND=r8):: eslift, ghtlift, qvplift, tmklift, tvlift
    REAL(KIND=r8):: ghtpari, prspari, qvppari, tmkpari
    REAL(KIND=r8):: tmkenv, qvpenv, tlcl
    REAL(KIND=r8):: fac1, fac2, facden, th, deltap
    REAL(KIND=r8):: benamin, davg, pavg, pressure, temp
    REAL(KIND=r8):: e, eth, ethmax, q, dz, cpm
    REAL(KIND=r8), DIMENSION(0:150):: buoy, zrel, benaccum
    REAL (KIND=r8), PARAMETER   :: gasr  =                  287.05_r8! gas constant of dry air        (j/kg/k)
    REAL (KIND=r8), PARAMETER   :: grav  =                   9.8e0_r8! gravity constant               (m/s**2)

  !p-> pressure in mbar
    delz=0.0_r8
    buoy(0:150)=0.0_r8
    zrel(0:150)=0.0_r8
    benaccum(0:150)=0.0_r8
    !! Get fields we want from the ones we have
    ter      = MAX(HGT,0.0_r8)
    qvp      = QV
    !prs      = PRES * 0.01_r8             ! pressure in hPa
    tmk      = TK                         ! temperature in K
    !ght      = GEOPT / G                  ! height in m
    DO k = 1, kMax
       DO i=1,nCols
            prs(i,k)      =    PSFC(i)*sl(k)          ! pressure in hPa
       END DO
    END DO
    DO i=1,nCols
        ght(i,1)=tmk(i,1)*(gasr/grav)*(PSFC(i)-prs(i,1))/PSFC(i)
        ght(i,1)=MAX(0.0_r8,ght(i,1))
        ght(i,1)=ter(i)+ght(i,1)
    END DO
    DO k=2,kMax
        DO i=1,nCols
           tv1=virtual(tmk(i,k-1),qvp(i,k-1))
           tv2=virtual(tmk(i,k),qvp(i,k))
           delz (i,k)=0.5_r8*gasr*(tv1+tv2)*LOG(prs(i,k-1)/prs(i,k))/grav
           ght  (i,k)=ght(i,k-1) + delz(i,k)
        END DO
    END DO
    !! First calculate a pressure array on FULL SIGMA levels
    !! Similar to the pfcalc.f routine from RIP4
    !! Top full sigma level is ommitted
    DO i=1,nCols
       prsf(i,1) = PSFC(i)             !! Lowest full sigma set to surface pressure
    END DO
    DO k = 2, kMax+1    
       DO i=1,nCols
          !prsf(:,k)=0.5_r8*(prs(:,k-1)+prs(:,k))
          prsf(i,k)=PSFC(i) *si(k)
       END DO
    END DO

    cape = 0.0_r8
    cin  = 0.0_r8

    !     DO j = 1,south_north_dim          !! BIG i/j loop
    DO i = 1,nCols

       IF ( i3dflag == 1 ) THEN       !! 3D case

          kpar1=kMax-1
          kpar2=1

       ELSE                           !! 2D case

          !      Find parcel with max theta-e in lowest 3 km AGL.
          ethmax = -1.0_r8
          DO k = 1, kMax
             IF ( ght(i,k)-ter(i) .LT. 3000.0_r8 ) THEN
                q        = MIN(MAX(qvp(i,k),1.e-12_r8),0.8_r8)
                temp     = tmk(i,k)
                pressure = prs(i,k)
                e        = MAX((q*pressure)/(EPS+q),1.e-9_r8)
                print*, 'erg temp', temp
                tlcl     = MAX(TLCLC1/MAX((LOG((temp**TLCLC2)/e)-TLCLC3),1.0e-12_r8)+TLCLC4,1.e-12_r8)
                eth  =     temp*(PSFC(i)/pressure)**( GAMMA_RIP*(1.0_r8+GAMMAMD*q) )*     &
                     EXP( (THTECON1/tlcl-THTECON2)*q*(1.0_r8+THTECON3*q) )
                IF ( eth .GT. ethmax ) THEN
                   klev=k
                   ethmax=eth
                END IF
             END IF
          END DO
          IF(klev <=1 ) klev=2  
          kpar1=klev
          kpar2=klev

          !      Establish average properties of that parcel
          !      (over depth of approximately davg meters)
          davg = 500.0_r8
          pavg = davg*prs(i,kpar1)*G /                        &
               ( Rd*virtual(tmk(i,kpar1),qvp(i,kpar1)) )
          p2 = MIN(prs(i,kpar1)+0.5_r8*pavg,prsf(i,1))
          p1 = p2-pavg
          totthe = 0.0_r8
          totqvp = 0.0_r8
          totprs = 0.0_r8
          DO k = 1,kMax-1
             IF ( prsf(i,k)   .LE. p1 ) GOTO 35
             IF ( prsf(i,k+1) .GE. p2 ) GOTO 34
             pressure = prs(i,k)
             pup      = prsf(i,k)
             pdn      = prsf(i,k+1)
             q        = MAX(qvp(i,k),1.e-15_r8)
             th       = tmk(i,k)*(PSFC(i)/pressure)**(GAMMA_RIP*(1.0_r8+GAMMAMD*q))
             pp1      = MAX(p1,pdn)
             pp2      = MIN(p2,pup)
             IF ( pp2 .GT. pp1 ) THEN
                deltap = pp2-pp1
                totqvp = totqvp+q*deltap
                totthe = totthe+th*deltap
                totprs = totprs+deltap
             END IF
34           CONTINUE
          END DO
35        CONTINUE
          print*, 'erg: q, deltap', q, deltap
          print*, 'erg: totqvp, totprs', totqvp, totprs
          qvppari = totqvp/totprs
          tmkpari = (totthe/totprs)*(prs(i,kpar1)/PSFC(i))**    &
               (GAMMA_RIP*(1.0_r8+GAMMAMD*qvp(i,kpar1)))
       ENDIF

!!!   END of 2D / 3D specific bits 


       DO kpar = kpar1,kpar2,-1                      !! This loop is done for both 2D / 3D

          !   Calculate temperature and moisture properties of parcel

          IF ( i3dflag == 1 ) THEN    ! (Note, qvppari and tmkpari already calculated above for 2D case.)
             !print*, 'erg inside i2dflag'
             qvppari = qvp(i,kpar)
             tmkpari = tmk(i,kpar)
          END IF
          prspari = prs(i,kpar)
          ghtpari = ght(i,kpar)
          cpm     = cp*(1.0_r8+CPMD*qvppari)
          e       = MAX(1.e-9_r8,qvppari*prspari/(EPS+qvppari))
          tlcl    = MAX(TLCLC1/(LOG(tmkpari**TLCLC2/e)-TLCLC3)+TLCLC4,1.0e-12_r8)
          !print*, 'erg', e
          !print*, 'erg', THTECON1, tlcl, THTECON2, qvppari, THTECON3, THTECON3*qvppari
          !print*, 'erg', kpar1,kpar2,kpar, (THTECON1/tlcl-THTECON2), qvppari, (1.0_r8+THTECON3*qvppari)
          !stop
          ethpari = tmkpari*(PSFC(i)/prspari)**(GAMMA_RIP*(1.0_r8+GAMMAMD*qvppari))*   &
               EXP((THTECON1/tlcl-THTECON2)*qvppari*                    &
               (1.0_r8+THTECON3*qvppari))
          zlcl    = ghtpari+(tmkpari-tlcl)/(G/cpm)

          !   Calculate buoyancy and relative height of lifted parcel at
          !   all levels, and store in bottom up arrays.  Add a level at the LCL,
          !   and at all points where buoyancy is zero.

          kk = 0                    ! for arrays that go bottom to top
          ilcl = 0
          IF ( ghtpari .GE. zlcl ) THEN
             !! initial parcel already saturated or supersaturated.
             ilcl = 2
             klcl = 1
          END IF

          DO k = kpar,kMax
33           kk=kk+1                ! for arrays that go bottom to top

             IF ( ght(i,k) .LT. zlcl ) THEN ! model level is below LCL
                qvplift = qvppari
                tmklift = tmkpari-G/cpm*(ght(i,k)-ghtpari)
                tvenv   = virtual(tmk(i,k),qvp(i,k))
                tvlift  = virtual(tmklift,qvplift)
                ghtlift = ght(i,k)
             ELSE IF ( ght(i,k) .GE. zlcl .AND. ilcl .EQ. 0 ) THEN
                !! This model level and previous model level straddle the LCL,
                !! so first create a new level in the bottom-up array, at the LCL.
                tmklift = tlcl
                qvplift = qvppari
                facden  = ght(i,k)-ght(i,k-1)
                fac1    = (zlcl-ght(i,k-1))/facden
                fac2    = (ght(i,k)-zlcl)/facden
                tmkenv  = tmk(i,k-1)*fac2+tmk(i,k)*fac1
                qvpenv  = qvp(i,k-1)*fac2+qvp(i,k)*fac1
                tvenv   = virtual(tmkenv,qvpenv)
                tvlift  = virtual(tmklift,qvplift)
                ghtlift = zlcl
                ilcl    = 1
             ELSE
                tmklift = tonpsadiabat(ethpari,prs(i,k),PSFC(i))                                
                eslift  = EZERO*EXP(eslcon1*(tmklift-CELKEL)/    &
                     MAX(tmklift-eslcon2,1.0e-12_r8))
                qvplift = EPS*eslift/(prs(i,k)-eslift)
                tvenv   = virtual(tmk(i,k),qvp(i,k))
                tvlift  = virtual(tmklift,qvplift)
                ghtlift = ght(i,k)
             END IF

             buoy(kk) = G*(tvlift-tvenv)/tvenv  ! buoyancy
             zrel(kk) = ghtlift-ghtpari

             IF ( (buoy(kk)*buoy(kk-1).LT.0.0_r8) .AND. (kk.GT.1) ) THEN
                !! Parcel ascent curve crosses sounding curve, so create a new level
                !! in the bottom-up array at the crossing.
                kk = kk+1
                buoy(kk)   = buoy(kk-1)
                zrel(kk)   = zrel(kk-1)
                buoy(kk-1) = 0.0_r8
                zrel(kk-1) = zrel(kk-2) + buoy(kk-2)/(buoy(kk-2)-buoy(kk))*  &
                     (zrel(kk)-zrel(kk-2))
             END IF

             IF (ilcl == 1) THEN
                klcl = kk
                ilcl = 2
                GOTO 33
             END IF

          END DO         !! END DO k = kpar,kLevel

          kLevel = kk
          IF (kLevel .GT. 150) THEN
             PRINT*, 'in calc_cape: kLevel got too big. kLevel=',kLevel
             STOP
          ENDIF


          !        Get the accumulated buoyant energy from the parcel's starting
          !        point, at all levels up to the top level.

          benaccum(1) = 0.0_r8
          benamin     = 9e9_r8
          DO k = 2,kLevel
             dz          = zrel(k)-zrel(k-1)
             benaccum(k) = benaccum(k-1)+0.5_r8*dz*(buoy(k-1)+buoy(k))
             IF ( benaccum(k) .LT. benamin ) THEN
                benamin = benaccum(k)
             END IF
          END DO


          !        Determine equilibrium level (EL), which we define as the highest
          !        level of non-negative buoyancy above the LCL. Note, this may be
          !        the top level if the parcel is still buoyant there.

          DO k = kLevel,klcl,-1
             IF ( buoy(k) .GE. 0.0_r8 ) THEN
                kel = k   ! k of equilibrium level
                GOTO 50
             END IF
          END DO


          !        If we got through that loop, then there is no non-negative
          !        buoyancy above the LCL in the sounding.  In these situations,
          !        both CAPE and CIN will be set to -0.1 J/kg.  Also, where CAPE is
          !        non-zero, CAPE and CIN will be set to a minimum of +0.1 J/kg, so
          !        that the zero contour in either the CIN or CAPE fields will
          !        circumscribe regions of non-zero CAPE.

          cape(i,kpar) = -0.1_r8
          cin(i,kpar)  = -0.1_r8
          klfc = kLevel

          GOTO 102

50        CONTINUE



          !        If there is an equilibrium level, then CAPE is positive.  We'll
          !        define the level of free convection (LFC) as the point below the
          !        EL, but at or above the LCL, where accumulated buoyant energy is a
          !        minimum.  The net positive area (accumulated buoyant energy) from
          !        the LFC up to the EL will be defined as the CAPE, and the net
          !        negative area (negative of accumulated buoyant energy) from the
          !        parcel starting point to the LFC will be defined as the convective
          !        inhibition (CIN).

          !        First get the LFC according to the above definition.

          benamin = 9e9_r8
          klfc = kLevel
          DO k = klcl,kel
             IF ( benaccum(k) .LT. benamin ) THEN
                benamin = benaccum(k)
                klfc = k
             END IF
          END DO

          !        Now we can assign values to cape and cin

          cape(i,kpar) = MAX(benaccum(kel)-benamin,0.1_r8)
          cin(i,kpar)  = MAX(-benamin,0.1_r8)

          !        CIN is uninteresting when CAPE is small (< 100 J/kg), so set
          !        CIN to -.1 in that case.

          IF ( cape(i,kpar) .LT. 100.0_r8 ) cin(i,kpar) = -0.1_r8

102       CONTINUE

       ENDDO          !! END of BIG 2D/3D loop


       IF ( i3dflag == 0 ) THEN
          SCRa(i,1) = cape(i,kpar1)
          SCRa(i,2) = cin(i,kpar1)
          SCRa(i,3) = zrel(klcl)+ghtpari-ter(i)   ! meters AGL (LCL)
          SCRa(i,4) = zrel(klfc)+ghtpari-ter(i)   ! meters AGL (LFC)
       ENDIF


    END DO
    !     END DO                !! END BIG i/j loop


    !! These will be set by module_diagnostics as we have more than 1 field

    IF ( i3dflag == 1 ) THEN
       SCRa = cape
       SCRb = cin
    ENDIF

    cname    = " "
    cdesc    = " "
    cunits   = " "


  END SUBROUTINE calc_cape


  !*********************************************************************c
  !*********************************************************************c
  REAL(KIND=r8) FUNCTION tonpsadiabat (thte,prs,PSFC)

    !  USE constants_module
    REAL(KIND=r8), PARAMETER :: Rd = 287.04_r8
    REAL(KIND=r8), PARAMETER :: Cp = 7.0_r8*Rd/2.0_r8
    REAL(KIND=r8), PARAMETER :: GAMMA_RIP = Rd/Cp 
    REAL(KIND=r8), INTENT(IN  ) ::  thte,prs,PSFC
    REAL(KIND=r8) :: fracip
    REAL(KIND=r8) :: fracip2
    REAL(KIND=r8) :: fracjt
    REAL(KIND=r8) :: fracjt2
    INTEGER       :: jtch
    INTEGER       :: jt
    INTEGER       :: ipch
    INTEGER       :: ip

    !   This function gives the temperature (in K) on a moist adiabat
    !   (specified by thte in K) given pressure in hPa.  It uses a
    !   lookup table, with data that was generated by the Bolton (1980)
    !   formula for theta_e.


    !     First check if pressure is less than min pressure in lookup table.
    !     If it is, assume parcel is so dry that the given theta-e value can
    !     be interpretted as theta, and get temperature from the simple dry
    !     theta formula.

    IF ( prs .LE. psadiprs(150) .or. prs .GE. psadiprs(1) ) THEN
       tonpsadiabat = thte*(prs/PSFC)**GAMMA_RIP
       RETURN
    ENDIF

    IF ( thte.GE.psadithte(150) .or. thte.LE.psadithte(1)) THEN
       tonpsadiabat = thte*(prs/PSFC)**GAMMA_RIP
       RETURN
    ENDIF

    !     Otherwise, look for the given thte/prs point in the lookup table.

    DO jtch = 1,150-1
       IF ( thte.GE.psadithte(jtch) .AND. thte.LT.psadithte(jtch+1) ) THEN
          jt = jtch
          GOTO 213
       END IF
    END DO
    jt = -1
213 CONTINUE

    DO ipch = 1,150-1
       IF ( prs.LE.psadiprs(ipch) .AND. prs.GT.psadiprs(ipch+1) ) THEN
          ip = ipch
          GOTO 215
       END IF
    ENDDO
    ip = -1
215 CONTINUE


    IF ( jt.EQ.-1 .OR. ip.EQ.-1 ) THEN
       PRINT*, 'Outside of lookup table bounds. prs,thte=',prs,thte
       STOP 
    ENDIF


    fracjt  = (thte-psadithte(jt))/(psadithte(jt+1)-psadithte(jt))
    fracjt2 = 1.0_r8-fracjt
    fracip  = (psadiprs(ip)-prs)/(psadiprs(ip)-psadiprs(ip+1))
    fracip2 = 1.0_r8-fracip

    IF ( psaditmk(ip,jt  ).GT.1e9_r8 .OR. psaditmk(ip+1,jt  ).GT.1e9_r8 .OR.   &
         psaditmk(ip,jt+1).GT.1e9_r8 .OR. psaditmk(ip+1,jt+1).GT.1e9_r8 ) THEN
       PRINT*, 'Tried to access missing tmperature in lookup table.'
       PRINT*, 'Prs and Thte probably unreasonable. prs,thte=',prs,thte
       tonpsadiabat = thte*(prs/PSFC)**GAMMA_RIP
       RETURN
    ENDIF

    tonpsadiabat = fracip2*fracjt2*psaditmk(ip  ,jt  )+   &
         fracip *fracjt2*psaditmk(ip+1,jt  )+   &
         fracip2*fracjt *psaditmk(ip  ,jt+1)+   &
         fracip *fracjt *psaditmk(ip+1,jt+1)


  END FUNCTION tonpsadiabat

  FUNCTION virtual (tmp,rmix)
    !      This function returns virtual temperature in K, given temperature
    !      in K and mixing ratio in kg/kg.

    REAL(KIND=r8)                              :: tmp, rmix, virtual

    virtual=tmp*(0.622_r8+rmix)/(0.622_r8*(1.0_r8+rmix))

  END FUNCTION virtual



  SUBROUTINE RD()
    IMPLICIT NONE


    psadithte(1:150)= (/ &
         0.2000000E+03_r8,     0.2016653E+03_r8,     0.2033444E+03_r8,     0.2050374E+03_r8,     0.2067446E+03_r8,&  
         0.2084660E+03_r8,     0.2102018E+03_r8,     0.2119519E+03_r8,     0.2137168E+03_r8,     0.2154962E+03_r8,&  
         0.2172905E+03_r8,     0.2190997E+03_r8,     0.2209240E+03_r8,     0.2227634E+03_r8,     0.2246182E+03_r8,&  
         0.2264884E+03_r8,     0.2283742E+03_r8,     0.2302757E+03_r8,     0.2321930E+03_r8,     0.2341263E+03_r8,&  
         0.2360757E+03_r8,     0.2380413E+03_r8,     0.2400232E+03_r8,     0.2420218E+03_r8,     0.2440370E+03_r8,&  
         0.2460689E+03_r8,     0.2481177E+03_r8,     0.2501835E+03_r8,     0.2522666E+03_r8,     0.2543671E+03_r8,&  
         0.2564850E+03_r8,     0.2586205E+03_r8,     0.2607738E+03_r8,     0.2629451E+03_r8,     0.2651344E+03_r8,&  
         0.2673420E+03_r8,     0.2695679E+03_r8,     0.2718124E+03_r8,     0.2740757E+03_r8,     0.2763577E+03_r8,&  
         0.2786587E+03_r8,     0.2809789E+03_r8,     0.2833184E+03_r8,     0.2856773E+03_r8,     0.2880559E+03_r8,&  
         0.2904543E+03_r8,     0.2928727E+03_r8,     0.2953112E+03_r8,     0.2977700E+03_r8,     0.3002493E+03_r8,&  
         0.3027493E+03_r8,     0.3052700E+03_r8,     0.3078117E+03_r8,     0.3103748E+03_r8,     0.3129590E+03_r8,&  
         0.3155648E+03_r8,     0.3181922E+03_r8,     0.3208416E+03_r8,     0.3235130E+03_r8,     0.3262066E+03_r8,&  
         0.3289226E+03_r8,     0.3316613E+03_r8,     0.3344228E+03_r8,     0.3372073E+03_r8,     0.3400149E+03_r8,&  
         0.3428459E+03_r8,     0.3457006E+03_r8,     0.3485789E+03_r8,     0.3514814E+03_r8,     0.3544079E+03_r8,&  
         0.3573588E+03_r8,     0.3603342E+03_r8,     0.3633344E+03_r8,     0.3663596E+03_r8,     0.3694100E+03_r8,&  
         0.3724858E+03_r8,     0.3755872E+03_r8,     0.3787144E+03_r8,     0.3818676E+03_r8,     0.3850471E+03_r8,&  
         0.3882531E+03_r8,     0.3914858E+03_r8,     0.3947454E+03_r8,     0.3980323E+03_r8,     0.4013464E+03_r8,&  
         0.4046881E+03_r8,     0.4080576E+03_r8,     0.4114552E+03_r8,     0.4148810E+03_r8,     0.4183354E+03_r8,&  
         0.4218185E+03_r8,     0.4253307E+03_r8,     0.4288721E+03_r8,     0.4324429E+03_r8,     0.4360435E+03_r8,&  
         0.4396741E+03_r8,     0.4433349E+03_r8,     0.4470262E+03_r8,     0.4507485E+03_r8,     0.4545015E+03_r8,&  
         0.4582858E+03_r8,     0.4621015E+03_r8,     0.4659491E+03_r8,     0.4698286E+03_r8,     0.4737405E+03_r8,&  
         0.4776850E+03_r8,     0.4816623E+03_r8,     0.4856727E+03_r8,     0.4897165E+03_r8,     0.4937940E+03_r8,&  
         0.4979054E+03_r8,     0.5020511E+03_r8,     0.5062312E+03_r8,     0.5104465E+03_r8,     0.5146965E+03_r8,&  
         0.5189820E+03_r8,     0.5233032E+03_r8,     0.5276603E+03_r8,     0.5320536E+03_r8,     0.5364836E+03_r8,&  
         0.5409505E+03_r8,     0.5454546E+03_r8,     0.5499962E+03_r8,     0.5545755E+03_r8,     0.5591930E+03_r8,&  
         0.5638489E+03_r8,     0.5685437E+03_r8,     0.5732775E+03_r8,     0.5780510E+03_r8,     0.5828636E+03_r8,&  
         0.5877170E+03_r8,     0.5926104E+03_r8,     0.5975446E+03_r8,     0.6025199E+03_r8,     0.6075366E+03_r8,&  
         0.6125950E+03_r8,     0.6176956E+03_r8,     0.6228387E+03_r8,     0.6280245E+03_r8,     0.6332536E+03_r8,&  
         0.6385262E+03_r8,     0.6438427E+03_r8,     0.6492034E+03_r8,     0.6546091E+03_r8,     0.6600593E+03_r8,&  
         0.6655554E+03_r8,     0.6710969E+03_r8,     0.6766846E+03_r8,     0.6823188E+03_r8,     0.6879999E+03_r8/)  

    psadiprs(1:150)= (/ &
         0.1100000E+04_r8,     0.1075803E+04_r8,     0.1052138E+04_r8,     0.1028994E+04_r8,     0.1006359E+04_r8,&  
         0.9842219E+03_r8,     0.9625715E+03_r8,     0.9413974E+03_r8,     0.9206895E+03_r8,     0.9004366E+03_r8,&  
         0.8806293E+03_r8,     0.8612581E+03_r8,     0.8423126E+03_r8,     0.8237839E+03_r8,     0.8056631E+03_r8,&  
         0.7879406E+03_r8,     0.7706078E+03_r8,     0.7536568E+03_r8,     0.7370782E+03_r8,     0.7208644E+03_r8,&  
         0.7050075E+03_r8,     0.6894991E+03_r8,     0.6743322E+03_r8,     0.6594986E+03_r8,     0.6449913E+03_r8,&  
         0.6308035E+03_r8,     0.6169274E+03_r8,     0.6033565E+03_r8,     0.5900845E+03_r8,     0.5771041E+03_r8,&  
         0.5644093E+03_r8,     0.5519940E+03_r8,     0.5398515E+03_r8,     0.5279761E+03_r8,     0.5163622E+03_r8,&  
         0.5050036E+03_r8,     0.4938948E+03_r8,     0.4830306E+03_r8,     0.4724051E+03_r8,     0.4620134E+03_r8,&  
         0.4518505E+03_r8,     0.4419109E+03_r8,     0.4321902E+03_r8,     0.4226831E+03_r8,     0.4133852E+03_r8,&  
         0.4042919E+03_r8,     0.3953985E+03_r8,     0.3867007E+03_r8,     0.3781945E+03_r8,     0.3698752E+03_r8,&  
         0.3617389E+03_r8,     0.3537817E+03_r8,     0.3459994E+03_r8,     0.3383883E+03_r8,     0.3309447E+03_r8,&  
         0.3236648E+03_r8,     0.3165450E+03_r8,     0.3095819E+03_r8,     0.3027719E+03_r8,     0.2961117E+03_r8,&  
         0.2895981E+03_r8,     0.2832277E+03_r8,     0.2769974E+03_r8,     0.2709043E+03_r8,     0.2649451E+03_r8,&  
         0.2591170E+03_r8,     0.2534172E+03_r8,     0.2478426E+03_r8,     0.2423907E+03_r8,     0.2370589E+03_r8,&  
         0.2318442E+03_r8,     0.2267442E+03_r8,     0.2217565E+03_r8,     0.2168784E+03_r8,     0.2121076E+03_r8,&  
         0.2074419E+03_r8,     0.2028787E+03_r8,     0.1984159E+03_r8,     0.1940513E+03_r8,     0.1897827E+03_r8,&  
         0.1856079E+03_r8,     0.1815251E+03_r8,     0.1775320E+03_r8,     0.1736268E+03_r8,     0.1698075E+03_r8,&  
         0.1660722E+03_r8,     0.1624190E+03_r8,     0.1588463E+03_r8,     0.1553521E+03_r8,     0.1519347E+03_r8,&  
         0.1485926E+03_r8,     0.1453240E+03_r8,     0.1421272E+03_r8,     0.1390008E+03_r8,     0.1359432E+03_r8,&  
         0.1329527E+03_r8,     0.1300282E+03_r8,     0.1271679E+03_r8,     0.1243705E+03_r8,     0.1216347E+03_r8,&  
         0.1189591E+03_r8,     0.1163423E+03_r8,     0.1137831E+03_r8,     0.1112802E+03_r8,     0.1088323E+03_r8,&  
         0.1064383E+03_r8,     0.1040970E+03_r8,     0.1018071E+03_r8,     0.9956760E+02_r8,     0.9737740E+02_r8,&  
         0.9523535E+02_r8,     0.9314041E+02_r8,     0.9109160E+02_r8,     0.8908781E+02_r8,     0.8712811E+02_r8,&  
         0.8521156E+02_r8,     0.8333711E+02_r8,     0.8150391E+02_r8,     0.7971107E+02_r8,     0.7795763E+02_r8,&  
         0.7624275E+02_r8,     0.7456564E+02_r8,     0.7292538E+02_r8,     0.7132121E+02_r8,     0.6975236E+02_r8,&  
         0.6821799E+02_r8,     0.6671736E+02_r8,     0.6524978E+02_r8,     0.6381445E+02_r8,     0.6241069E+02_r8,&  
         0.6103785E+02_r8,     0.5969517E+02_r8,     0.5838202E+02_r8,     0.5709779E+02_r8,     0.5584179E+02_r8,&  
         0.5461341E+02_r8,     0.5341208E+02_r8,     0.5223714E+02_r8,     0.5108807E+02_r8,     0.4996428E+02_r8,&  
         0.4886519E+02_r8,     0.4779029E+02_r8,     0.4673904E+02_r8,     0.4571090E+02_r8,     0.4470538E+02_r8,&  
         0.4372198E+02_r8,     0.4276021E+02_r8,     0.4181961E+02_r8,     0.4089969E+02_r8,     0.4000000E+02_r8/)  

    psaditmk(1:150,    1)= (/ &
         0.2055137E+03_r8,     0.2042125E+03_r8,     0.2029193E+03_r8,     0.2016342E+03_r8,     0.2003570E+03_r8,&  
         0.1990878E+03_r8,     0.1978265E+03_r8,     0.1965731E+03_r8,     0.1953276E+03_r8,     0.1940898E+03_r8,&  
         0.1928599E+03_r8,     0.1916377E+03_r8,     0.1904232E+03_r8,     0.1892163E+03_r8,     0.1880171E+03_r8,&  
         0.1868254E+03_r8,     0.1856413E+03_r8,     0.1844646E+03_r8,     0.1832954E+03_r8,     0.1821335E+03_r8,&  
         0.1809791E+03_r8,     0.1798319E+03_r8,     0.1786920E+03_r8,     0.1775593E+03_r8,     0.1764338E+03_r8,&  
         0.1753154E+03_r8,     0.1742041E+03_r8,     0.1730999E+03_r8,     0.1720026E+03_r8,     0.1709123E+03_r8,&  
         0.1698289E+03_r8,     0.1687524E+03_r8,     0.1676826E+03_r8,     0.1666197E+03_r8,     0.1655635E+03_r8,&  
         0.1645140E+03_r8,     0.1634712E+03_r8,     0.1624350E+03_r8,     0.1614053E+03_r8,     0.1603821E+03_r8,&  
         0.1593655E+03_r8,     0.1583553E+03_r8,     0.1573515E+03_r8,     0.1563540E+03_r8,     0.1553629E+03_r8,&  
         0.1543781E+03_r8,     0.1533994E+03_r8,     0.1524270E+03_r8,     0.1514608E+03_r8,     0.1505007E+03_r8,&  
         0.1495467E+03_r8,     0.1485987E+03_r8,     0.1476568E+03_r8,     0.1467207E+03_r8,     0.1457907E+03_r8,&  
         0.1448665E+03_r8,     0.1439482E+03_r8,     0.1430358E+03_r8,     0.1421290E+03_r8,     0.1412281E+03_r8,&  
         0.1403329E+03_r8,     0.1394433E+03_r8,     0.1385594E+03_r8,     0.1376810E+03_r8,     0.1368083E+03_r8,&  
         0.1359410E+03_r8,     0.1350793E+03_r8,     0.1342231E+03_r8,     0.1333722E+03_r8,     0.1325268E+03_r8,&  
         0.1316867E+03_r8,     0.1308519E+03_r8,     0.1300225E+03_r8,     0.1291983E+03_r8,     0.1283793E+03_r8,&  
         0.1275655E+03_r8,     0.1267568E+03_r8,     0.1259533E+03_r8,     0.1251549E+03_r8,     0.1243616E+03_r8,&  
         0.1235732E+03_r8,     0.1227899E+03_r8,     0.1220116E+03_r8,     0.1212381E+03_r8,     0.1204696E+03_r8,&  
         0.1197059E+03_r8,     0.1189471E+03_r8,     0.1181931E+03_r8,     0.1174439E+03_r8,     0.1166994E+03_r8,&  
         0.1159597E+03_r8,     0.1152246E+03_r8,     0.1144942E+03_r8,     0.1137684E+03_r8,     0.1130473E+03_r8,&  
         0.1123306E+03_r8,     0.1116186E+03_r8,     0.1109110E+03_r8,     0.1102080E+03_r8,     0.1095094E+03_r8,&  
         0.1088152E+03_r8,     0.1081254E+03_r8,     0.1074400E+03_r8,     0.1067590E+03_r8,     0.1060822E+03_r8,&  
         0.1054098E+03_r8,     0.1047416E+03_r8,     0.1040776E+03_r8,     0.1034179E+03_r8,     0.1027623E+03_r8,&  
         0.1021109E+03_r8,     0.1014636E+03_r8,     0.1008204E+03_r8,     0.1001813E+03_r8,     0.9954629E+02_r8,&  
         0.9891528E+02_r8,     0.9828825E+02_r8,     0.9766520E+02_r8,     0.9704611E+02_r8,     0.9643093E+02_r8,&  
         0.9581966E+02_r8,     0.9521227E+02_r8,     0.9460871E+02_r8,     0.9400898E+02_r8,     0.9341307E+02_r8,&  
         0.9282092E+02_r8,     0.9223253E+02_r8,     0.9164787E+02_r8,     0.9106693E+02_r8,     0.9048965E+02_r8,&  
         0.8991605E+02_r8,     0.8934606E+02_r8,     0.8877969E+02_r8,     0.8821693E+02_r8,     0.8765772E+02_r8,&  
         0.8710206E+02_r8,     0.8654993E+02_r8,     0.8600128E+02_r8,     0.8545613E+02_r8,     0.8491442E+02_r8,&  
         0.8437615E+02_r8,     0.8384129E+02_r8,     0.8330983E+02_r8,     0.8278172E+02_r8,     0.8225697E+02_r8,&  
         0.8173555E+02_r8,     0.8121742E+02_r8,     0.8070259E+02_r8,     0.8019102E+02_r8,     0.7968269E+02_r8/)  

    psaditmk(1:150,    2)= (/ &
         0.2072220E+03_r8,     0.2059104E+03_r8,     0.2046068E+03_r8,     0.2033112E+03_r8,     0.2020237E+03_r8,&  
         0.2007441E+03_r8,     0.1994725E+03_r8,     0.1982088E+03_r8,     0.1969531E+03_r8,     0.1957052E+03_r8,&  
         0.1944651E+03_r8,     0.1932328E+03_r8,     0.1920082E+03_r8,     0.1907914E+03_r8,     0.1895822E+03_r8,&  
         0.1883807E+03_r8,     0.1871867E+03_r8,     0.1860003E+03_r8,     0.1848214E+03_r8,     0.1836499E+03_r8,&  
         0.1824858E+03_r8,     0.1813291E+03_r8,     0.1801797E+03_r8,     0.1790376E+03_r8,     0.1779028E+03_r8,&  
         0.1767751E+03_r8,     0.1756546E+03_r8,     0.1745411E+03_r8,     0.1734347E+03_r8,     0.1723353E+03_r8,&  
         0.1712429E+03_r8,     0.1701574E+03_r8,     0.1690788E+03_r8,     0.1680070E+03_r8,     0.1669420E+03_r8,&  
         0.1658838E+03_r8,     0.1648323E+03_r8,     0.1637874E+03_r8,     0.1627491E+03_r8,     0.1617175E+03_r8,&  
         0.1606924E+03_r8,     0.1596737E+03_r8,     0.1586616E+03_r8,     0.1576558E+03_r8,     0.1566565E+03_r8,&  
         0.1556634E+03_r8,     0.1546767E+03_r8,     0.1536962E+03_r8,     0.1527219E+03_r8,     0.1517538E+03_r8,&  
         0.1507918E+03_r8,     0.1498360E+03_r8,     0.1488862E+03_r8,     0.1479424E+03_r8,     0.1470046E+03_r8,&  
         0.1460727E+03_r8,     0.1451468E+03_r8,     0.1442267E+03_r8,     0.1433124E+03_r8,     0.1424040E+03_r8,&  
         0.1415013E+03_r8,     0.1406043E+03_r8,     0.1397130E+03_r8,     0.1388274E+03_r8,     0.1379474E+03_r8,&  
         0.1370729E+03_r8,     0.1362040E+03_r8,     0.1353406E+03_r8,     0.1344827E+03_r8,     0.1336302E+03_r8,&  
         0.1327831E+03_r8,     0.1319414E+03_r8,     0.1311051E+03_r8,     0.1302740E+03_r8,     0.1294482E+03_r8,&  
         0.1286276E+03_r8,     0.1278122E+03_r8,     0.1270020E+03_r8,     0.1261970E+03_r8,     0.1253970E+03_r8,&  
         0.1246021E+03_r8,     0.1238123E+03_r8,     0.1230274E+03_r8,     0.1222476E+03_r8,     0.1214727E+03_r8,&  
         0.1207026E+03_r8,     0.1199375E+03_r8,     0.1191772E+03_r8,     0.1184218E+03_r8,     0.1176711E+03_r8,&  
         0.1169252E+03_r8,     0.1161840E+03_r8,     0.1154475E+03_r8,     0.1147157E+03_r8,     0.1139885E+03_r8,&  
         0.1132659E+03_r8,     0.1125479E+03_r8,     0.1118345E+03_r8,     0.1111256E+03_r8,     0.1104212E+03_r8,&  
         0.1097212E+03_r8,     0.1090257E+03_r8,     0.1083346E+03_r8,     0.1076479E+03_r8,     0.1069655E+03_r8,&  
         0.1062874E+03_r8,     0.1056137E+03_r8,     0.1049442E+03_r8,     0.1042789E+03_r8,     0.1036179E+03_r8,&  
         0.1029611E+03_r8,     0.1023084E+03_r8,     0.1016599E+03_r8,     0.1010155E+03_r8,     0.1003751E+03_r8,&  
         0.9973886E+02_r8,     0.9910662E+02_r8,     0.9847838E+02_r8,     0.9785413E+02_r8,     0.9723383E+02_r8,&  
         0.9661746E+02_r8,     0.9600502E+02_r8,     0.9539644E+02_r8,     0.9479172E+02_r8,     0.9419085E+02_r8,&  
         0.9359377E+02_r8,     0.9300048E+02_r8,     0.9241096E+02_r8,     0.9182516E+02_r8,     0.9124308E+02_r8,&  
         0.9066470E+02_r8,     0.9008997E+02_r8,     0.8951889E+02_r8,     0.8895144E+02_r8,     0.8838757E+02_r8,&  
         0.8782729E+02_r8,     0.8727055E+02_r8,     0.8671735E+02_r8,     0.8616765E+02_r8,     0.8562144E+02_r8,&  
         0.8507868E+02_r8,     0.8453938E+02_r8,     0.8400348E+02_r8,     0.8347099E+02_r8,     0.8294186E+02_r8,&  
         0.8241610E+02_r8,     0.8189366E+02_r8,     0.8137453E+02_r8,     0.8085871E+02_r8,     0.8034615E+02_r8/)  

    psaditmk(1:150,    3)= (/ &
         0.2089439E+03_r8,     0.2076218E+03_r8,     0.2063078E+03_r8,     0.2050018E+03_r8,     0.2037038E+03_r8,&  
         0.2024139E+03_r8,     0.2011320E+03_r8,     0.1998580E+03_r8,     0.1985919E+03_r8,     0.1973338E+03_r8,&  
         0.1960835E+03_r8,     0.1948410E+03_r8,     0.1936064E+03_r8,     0.1923795E+03_r8,     0.1911603E+03_r8,&  
         0.1899489E+03_r8,     0.1887450E+03_r8,     0.1875488E+03_r8,     0.1863600E+03_r8,     0.1851788E+03_r8,&  
         0.1840051E+03_r8,     0.1828388E+03_r8,     0.1816799E+03_r8,     0.1805283E+03_r8,     0.1793840E+03_r8,&  
         0.1782469E+03_r8,     0.1771171E+03_r8,     0.1759943E+03_r8,     0.1748787E+03_r8,     0.1737702E+03_r8,&  
         0.1726687E+03_r8,     0.1715742E+03_r8,     0.1704866E+03_r8,     0.1694059E+03_r8,     0.1683320E+03_r8,&  
         0.1672650E+03_r8,     0.1662047E+03_r8,     0.1651511E+03_r8,     0.1641042E+03_r8,     0.1630640E+03_r8,&  
         0.1620303E+03_r8,     0.1610032E+03_r8,     0.1599826E+03_r8,     0.1589685E+03_r8,     0.1579608E+03_r8,&  
         0.1569595E+03_r8,     0.1559646E+03_r8,     0.1549759E+03_r8,     0.1539935E+03_r8,     0.1530173E+03_r8,&  
         0.1520474E+03_r8,     0.1510835E+03_r8,     0.1501258E+03_r8,     0.1491742E+03_r8,     0.1482286E+03_r8,&  
         0.1472890E+03_r8,     0.1463553E+03_r8,     0.1454276E+03_r8,     0.1445057E+03_r8,     0.1435897E+03_r8,&  
         0.1426795E+03_r8,     0.1417750E+03_r8,     0.1408763E+03_r8,     0.1399833E+03_r8,     0.1390959E+03_r8,&  
         0.1382142E+03_r8,     0.1373381E+03_r8,     0.1364675E+03_r8,     0.1356024E+03_r8,     0.1347429E+03_r8,&  
         0.1338887E+03_r8,     0.1330400E+03_r8,     0.1321967E+03_r8,     0.1313587E+03_r8,     0.1305260E+03_r8,&  
         0.1296986E+03_r8,     0.1288764E+03_r8,     0.1280595E+03_r8,     0.1272477E+03_r8,     0.1264411E+03_r8,&  
         0.1256396E+03_r8,     0.1248432E+03_r8,     0.1240518E+03_r8,     0.1232654E+03_r8,     0.1224841E+03_r8,&  
         0.1217076E+03_r8,     0.1209361E+03_r8,     0.1201695E+03_r8,     0.1194078E+03_r8,     0.1186508E+03_r8,&  
         0.1178987E+03_r8,     0.1171514E+03_r8,     0.1164087E+03_r8,     0.1156708E+03_r8,     0.1149376E+03_r8,&  
         0.1142090E+03_r8,     0.1134850E+03_r8,     0.1127657E+03_r8,     0.1120508E+03_r8,     0.1113406E+03_r8,&  
         0.1106348E+03_r8,     0.1099334E+03_r8,     0.1092366E+03_r8,     0.1085442E+03_r8,     0.1078561E+03_r8,&  
         0.1071724E+03_r8,     0.1064930E+03_r8,     0.1058180E+03_r8,     0.1051472E+03_r8,     0.1044807E+03_r8,&  
         0.1038184E+03_r8,     0.1031603E+03_r8,     0.1025063E+03_r8,     0.1018565E+03_r8,     0.1012109E+03_r8,&  
         0.1005693E+03_r8,     0.9993180E+02_r8,     0.9929833E+02_r8,     0.9866888E+02_r8,     0.9804343E+02_r8,&  
         0.9742192E+02_r8,     0.9680437E+02_r8,     0.9619073E+02_r8,     0.9558098E+02_r8,     0.9497511E+02_r8,&  
         0.9437305E+02_r8,     0.9377482E+02_r8,     0.9318039E+02_r8,     0.9258971E+02_r8,     0.9200278E+02_r8,&  
         0.9141959E+02_r8,     0.9084009E+02_r8,     0.9026424E+02_r8,     0.8969207E+02_r8,     0.8912351E+02_r8,&  
         0.8855856E+02_r8,     0.8799719E+02_r8,     0.8743938E+02_r8,     0.8688510E+02_r8,     0.8633434E+02_r8,&  
         0.8578706E+02_r8,     0.8524326E+02_r8,     0.8470291E+02_r8,     0.8416598E+02_r8,     0.8363245E+02_r8,&  
         0.8310231E+02_r8,     0.8257552E+02_r8,     0.8205208E+02_r8,     0.8153195E+02_r8,     0.8101512E+02_r8/)  

    psaditmk(1:150,    4)= (/ &
         0.2106792E+03_r8,     0.2093467E+03_r8,     0.2080223E+03_r8,     0.2067059E+03_r8,     0.2053975E+03_r8,&  
         0.2040972E+03_r8,     0.2028048E+03_r8,     0.2015205E+03_r8,     0.2002441E+03_r8,     0.1989757E+03_r8,&  
         0.1977151E+03_r8,     0.1964625E+03_r8,     0.1952177E+03_r8,     0.1939807E+03_r8,     0.1927515E+03_r8,&  
         0.1915300E+03_r8,     0.1903162E+03_r8,     0.1891100E+03_r8,     0.1879115E+03_r8,     0.1867205E+03_r8,&  
         0.1855370E+03_r8,     0.1843610E+03_r8,     0.1831925E+03_r8,     0.1820313E+03_r8,     0.1808775E+03_r8,&  
         0.1797310E+03_r8,     0.1785917E+03_r8,     0.1774596E+03_r8,     0.1763348E+03_r8,     0.1752170E+03_r8,&  
         0.1741063E+03_r8,     0.1730027E+03_r8,     0.1719061E+03_r8,     0.1708164E+03_r8,     0.1697336E+03_r8,&  
         0.1686576E+03_r8,     0.1675885E+03_r8,     0.1665262E+03_r8,     0.1654706E+03_r8,     0.1644217E+03_r8,&  
         0.1633794E+03_r8,     0.1623438E+03_r8,     0.1613147E+03_r8,     0.1602921E+03_r8,     0.1592760E+03_r8,&  
         0.1582664E+03_r8,     0.1572631E+03_r8,     0.1562663E+03_r8,     0.1552757E+03_r8,     0.1542914E+03_r8,&  
         0.1533133E+03_r8,     0.1523415E+03_r8,     0.1513758E+03_r8,     0.1504162E+03_r8,     0.1494628E+03_r8,&  
         0.1485153E+03_r8,     0.1475739E+03_r8,     0.1466384E+03_r8,     0.1457089E+03_r8,     0.1447852E+03_r8,&  
         0.1438674E+03_r8,     0.1429555E+03_r8,     0.1420493E+03_r8,     0.1411488E+03_r8,     0.1402541E+03_r8,&  
         0.1393650E+03_r8,     0.1384816E+03_r8,     0.1376037E+03_r8,     0.1367315E+03_r8,     0.1358647E+03_r8,&  
         0.1350035E+03_r8,     0.1341477E+03_r8,     0.1332974E+03_r8,     0.1324524E+03_r8,     0.1316128E+03_r8,&  
         0.1307785E+03_r8,     0.1299495E+03_r8,     0.1291257E+03_r8,     0.1283072E+03_r8,     0.1274939E+03_r8,&  
         0.1266857E+03_r8,     0.1258826E+03_r8,     0.1250847E+03_r8,     0.1242918E+03_r8,     0.1235039E+03_r8,&  
         0.1227210E+03_r8,     0.1219431E+03_r8,     0.1211701E+03_r8,     0.1204020E+03_r8,     0.1196387E+03_r8,&  
         0.1188804E+03_r8,     0.1181268E+03_r8,     0.1173780E+03_r8,     0.1166339E+03_r8,     0.1158946E+03_r8,&  
         0.1151599E+03_r8,     0.1144299E+03_r8,     0.1137046E+03_r8,     0.1129838E+03_r8,     0.1122676E+03_r8,&  
         0.1115559E+03_r8,     0.1108488E+03_r8,     0.1101461E+03_r8,     0.1094479E+03_r8,     0.1087541E+03_r8,&  
         0.1080647E+03_r8,     0.1073797E+03_r8,     0.1066990E+03_r8,     0.1060227E+03_r8,     0.1053506E+03_r8,&  
         0.1046828E+03_r8,     0.1040192E+03_r8,     0.1033598E+03_r8,     0.1027046E+03_r8,     0.1020536E+03_r8,&  
         0.1014067E+03_r8,     0.1007639E+03_r8,     0.1001251E+03_r8,     0.9949043E+02_r8,     0.9885976E+02_r8,&  
         0.9823308E+02_r8,     0.9761039E+02_r8,     0.9699163E+02_r8,     0.9637680E+02_r8,     0.9576588E+02_r8,&  
         0.9515882E+02_r8,     0.9455560E+02_r8,     0.9395622E+02_r8,     0.9336064E+02_r8,     0.9276882E+02_r8,&  
         0.9218077E+02_r8,     0.9159644E+02_r8,     0.9101580E+02_r8,     0.9043887E+02_r8,     0.8986557E+02_r8,&  
         0.8929591E+02_r8,     0.8872987E+02_r8,     0.8816741E+02_r8,     0.8760852E+02_r8,     0.8705318E+02_r8,&  
         0.8650134E+02_r8,     0.8595302E+02_r8,     0.8540816E+02_r8,     0.8486676E+02_r8,     0.8432879E+02_r8,&  
         0.8379424E+02_r8,     0.8326306E+02_r8,     0.8273526E+02_r8,     0.8221081E+02_r8,     0.8168967E+02_r8/)  

    psaditmk(1:150,    5)= (/ &
         0.2124279E+03_r8,     0.2110851E+03_r8,     0.2097502E+03_r8,     0.2084234E+03_r8,     0.2071046E+03_r8,&  
         0.2057939E+03_r8,     0.2044911E+03_r8,     0.2031965E+03_r8,     0.2019097E+03_r8,     0.2006310E+03_r8,&  
         0.1993601E+03_r8,     0.1980973E+03_r8,     0.1968422E+03_r8,     0.1955951E+03_r8,     0.1943557E+03_r8,&  
         0.1931241E+03_r8,     0.1919003E+03_r8,     0.1906842E+03_r8,     0.1894757E+03_r8,     0.1882748E+03_r8,&  
         0.1870816E+03_r8,     0.1858958E+03_r8,     0.1847176E+03_r8,     0.1835468E+03_r8,     0.1823834E+03_r8,&  
         0.1812273E+03_r8,     0.1800786E+03_r8,     0.1789371E+03_r8,     0.1778029E+03_r8,     0.1766758E+03_r8,&  
         0.1755559E+03_r8,     0.1744431E+03_r8,     0.1733374E+03_r8,     0.1722386E+03_r8,     0.1711468E+03_r8,&  
         0.1700619E+03_r8,     0.1689839E+03_r8,     0.1679127E+03_r8,     0.1668483E+03_r8,     0.1657907E+03_r8,&  
         0.1647397E+03_r8,     0.1636955E+03_r8,     0.1626578E+03_r8,     0.1616267E+03_r8,     0.1606022E+03_r8,&  
         0.1595841E+03_r8,     0.1585725E+03_r8,     0.1575673E+03_r8,     0.1565685E+03_r8,     0.1555760E+03_r8,&  
         0.1545898E+03_r8,     0.1536099E+03_r8,     0.1526362E+03_r8,     0.1516686E+03_r8,     0.1507072E+03_r8,&  
         0.1497519E+03_r8,     0.1488026E+03_r8,     0.1478593E+03_r8,     0.1469221E+03_r8,     0.1459907E+03_r8,&  
         0.1450653E+03_r8,     0.1441457E+03_r8,     0.1432320E+03_r8,     0.1423241E+03_r8,     0.1414219E+03_r8,&  
         0.1405254E+03_r8,     0.1396346E+03_r8,     0.1387495E+03_r8,     0.1378699E+03_r8,     0.1369960E+03_r8,&  
         0.1361275E+03_r8,     0.1352646E+03_r8,     0.1344072E+03_r8,     0.1335552E+03_r8,     0.1327086E+03_r8,&  
         0.1318674E+03_r8,     0.1310315E+03_r8,     0.1302009E+03_r8,     0.1293755E+03_r8,     0.1285554E+03_r8,&  
         0.1277405E+03_r8,     0.1269308E+03_r8,     0.1261261E+03_r8,     0.1253266E+03_r8,     0.1245322E+03_r8,&  
         0.1237428E+03_r8,     0.1229584E+03_r8,     0.1221789E+03_r8,     0.1214044E+03_r8,     0.1206349E+03_r8,&  
         0.1198702E+03_r8,     0.1191103E+03_r8,     0.1183553E+03_r8,     0.1176050E+03_r8,     0.1168596E+03_r8,&  
         0.1161188E+03_r8,     0.1153827E+03_r8,     0.1146513E+03_r8,     0.1139245E+03_r8,     0.1132024E+03_r8,&  
         0.1124848E+03_r8,     0.1117717E+03_r8,     0.1110632E+03_r8,     0.1103592E+03_r8,     0.1096596E+03_r8,&  
         0.1089645E+03_r8,     0.1082738E+03_r8,     0.1075874E+03_r8,     0.1069054E+03_r8,     0.1062278E+03_r8,&  
         0.1055544E+03_r8,     0.1048853E+03_r8,     0.1042204E+03_r8,     0.1035598E+03_r8,     0.1029033E+03_r8,&  
         0.1022510E+03_r8,     0.1016028E+03_r8,     0.1009588E+03_r8,     0.1003188E+03_r8,     0.9968288E+02_r8,&  
         0.9905098E+02_r8,     0.9842310E+02_r8,     0.9779919E+02_r8,     0.9717924E+02_r8,     0.9656324E+02_r8,&  
         0.9595113E+02_r8,     0.9534289E+02_r8,     0.9473852E+02_r8,     0.9413797E+02_r8,     0.9354123E+02_r8,&  
         0.9294828E+02_r8,     0.9235908E+02_r8,     0.9177361E+02_r8,     0.9119187E+02_r8,     0.9061380E+02_r8,&  
         0.9003940E+02_r8,     0.8946865E+02_r8,     0.8890150E+02_r8,     0.8833797E+02_r8,     0.8777799E+02_r8,&  
         0.8722157E+02_r8,     0.8666867E+02_r8,     0.8611929E+02_r8,     0.8557337E+02_r8,     0.8503092E+02_r8,&  
         0.8449192E+02_r8,     0.8395632E+02_r8,     0.8342413E+02_r8,     0.8289530E+02_r8,     0.8236983E+02_r8/)  

    psaditmk(1:150,    6)= (/ &
         0.2141900E+03_r8,     0.2128368E+03_r8,     0.2114916E+03_r8,     0.2101544E+03_r8,     0.2088253E+03_r8,&  
         0.2075041E+03_r8,     0.2061910E+03_r8,     0.2048859E+03_r8,     0.2035888E+03_r8,     0.2022997E+03_r8,&  
         0.2010186E+03_r8,     0.1997454E+03_r8,     0.1984801E+03_r8,     0.1972227E+03_r8,     0.1959732E+03_r8,&  
         0.1947315E+03_r8,     0.1934975E+03_r8,     0.1922714E+03_r8,     0.1910529E+03_r8,     0.1898421E+03_r8,&  
         0.1886390E+03_r8,     0.1874434E+03_r8,     0.1862554E+03_r8,     0.1850749E+03_r8,     0.1839018E+03_r8,&  
         0.1827361E+03_r8,     0.1815779E+03_r8,     0.1804269E+03_r8,     0.1792833E+03_r8,     0.1781469E+03_r8,&  
         0.1770176E+03_r8,     0.1758956E+03_r8,     0.1747806E+03_r8,     0.1736727E+03_r8,     0.1725718E+03_r8,&  
         0.1714778E+03_r8,     0.1703909E+03_r8,     0.1693108E+03_r8,     0.1682375E+03_r8,     0.1671711E+03_r8,&  
         0.1661114E+03_r8,     0.1650584E+03_r8,     0.1640121E+03_r8,     0.1629725E+03_r8,     0.1619394E+03_r8,&  
         0.1609129E+03_r8,     0.1598928E+03_r8,     0.1588793E+03_r8,     0.1578722E+03_r8,     0.1568714E+03_r8,&  
         0.1558770E+03_r8,     0.1548889E+03_r8,     0.1539071E+03_r8,     0.1529314E+03_r8,     0.1519620E+03_r8,&  
         0.1509987E+03_r8,     0.1500415E+03_r8,     0.1490905E+03_r8,     0.1481454E+03_r8,     0.1472063E+03_r8,&  
         0.1462731E+03_r8,     0.1453459E+03_r8,     0.1444246E+03_r8,     0.1435091E+03_r8,     0.1425994E+03_r8,&  
         0.1416954E+03_r8,     0.1407972E+03_r8,     0.1399047E+03_r8,     0.1390179E+03_r8,     0.1381366E+03_r8,&  
         0.1372610E+03_r8,     0.1363909E+03_r8,     0.1355263E+03_r8,     0.1346672E+03_r8,     0.1338136E+03_r8,&  
         0.1329653E+03_r8,     0.1321225E+03_r8,     0.1312849E+03_r8,     0.1304527E+03_r8,     0.1296258E+03_r8,&  
         0.1288041E+03_r8,     0.1279876E+03_r8,     0.1271763E+03_r8,     0.1263701E+03_r8,     0.1255691E+03_r8,&  
         0.1247731E+03_r8,     0.1239822E+03_r8,     0.1231962E+03_r8,     0.1224153E+03_r8,     0.1216393E+03_r8,&  
         0.1208682E+03_r8,     0.1201021E+03_r8,     0.1193407E+03_r8,     0.1185842E+03_r8,     0.1178326E+03_r8,&  
         0.1170856E+03_r8,     0.1163434E+03_r8,     0.1156059E+03_r8,     0.1148731E+03_r8,     0.1141449E+03_r8,&  
         0.1134213E+03_r8,     0.1127024E+03_r8,     0.1119879E+03_r8,     0.1112781E+03_r8,     0.1105727E+03_r8,&  
         0.1098717E+03_r8,     0.1091753E+03_r8,     0.1084832E+03_r8,     0.1077955E+03_r8,     0.1071122E+03_r8,&  
         0.1064333E+03_r8,     0.1057586E+03_r8,     0.1050882E+03_r8,     0.1044220E+03_r8,     0.1037601E+03_r8,&  
         0.1031024E+03_r8,     0.1024488E+03_r8,     0.1017994E+03_r8,     0.1011541E+03_r8,     0.1005129E+03_r8,&  
         0.9987570E+02_r8,     0.9924260E+02_r8,     0.9861349E+02_r8,     0.9798839E+02_r8,     0.9736725E+02_r8,&  
         0.9675004E+02_r8,     0.9613673E+02_r8,     0.9552733E+02_r8,     0.9492178E+02_r8,     0.9432008E+02_r8,&  
         0.9372219E+02_r8,     0.9312808E+02_r8,     0.9253774E+02_r8,     0.9195116E+02_r8,     0.9136827E+02_r8,&  
         0.9078909E+02_r8,     0.9021358E+02_r8,     0.8964172E+02_r8,     0.8907349E+02_r8,     0.8850885E+02_r8,&  
         0.8794779E+02_r8,     0.8739030E+02_r8,     0.8683633E+02_r8,     0.8628588E+02_r8,     0.8573891E+02_r8,&  
         0.8519542E+02_r8,     0.8465536E+02_r8,     0.8411874E+02_r8,     0.8358551E+02_r8,     0.8305566E+02_r8/)  

    psaditmk(1:150,    7)= (/ &
         0.2159651E+03_r8,     0.2146017E+03_r8,     0.2132463E+03_r8,     0.2118988E+03_r8,     0.2105593E+03_r8,&  
         0.2092278E+03_r8,     0.2079043E+03_r8,     0.2065888E+03_r8,     0.2052813E+03_r8,     0.2039818E+03_r8,&  
         0.2026904E+03_r8,     0.2014068E+03_r8,     0.2001312E+03_r8,     0.1988636E+03_r8,     0.1976039E+03_r8,&  
         0.1963520E+03_r8,     0.1951079E+03_r8,     0.1938716E+03_r8,     0.1926431E+03_r8,     0.1914223E+03_r8,&  
         0.1902092E+03_r8,     0.1890038E+03_r8,     0.1878059E+03_r8,     0.1866156E+03_r8,     0.1854328E+03_r8,&  
         0.1842575E+03_r8,     0.1830896E+03_r8,     0.1819291E+03_r8,     0.1807759E+03_r8,     0.1796301E+03_r8,&  
         0.1784914E+03_r8,     0.1773600E+03_r8,     0.1762358E+03_r8,     0.1751186E+03_r8,     0.1740086E+03_r8,&  
         0.1729056E+03_r8,     0.1718096E+03_r8,     0.1707205E+03_r8,     0.1696383E+03_r8,     0.1685630E+03_r8,&  
         0.1674945E+03_r8,     0.1664327E+03_r8,     0.1653777E+03_r8,     0.1643294E+03_r8,     0.1632877E+03_r8,&  
         0.1622527E+03_r8,     0.1612241E+03_r8,     0.1602021E+03_r8,     0.1591866E+03_r8,     0.1581776E+03_r8,&  
         0.1571749E+03_r8,     0.1561785E+03_r8,     0.1551885E+03_r8,     0.1542048E+03_r8,     0.1532273E+03_r8,&  
         0.1522560E+03_r8,     0.1512908E+03_r8,     0.1503318E+03_r8,     0.1493789E+03_r8,     0.1484319E+03_r8,&  
         0.1474910E+03_r8,     0.1465561E+03_r8,     0.1456271E+03_r8,     0.1447040E+03_r8,     0.1437867E+03_r8,&  
         0.1428752E+03_r8,     0.1419695E+03_r8,     0.1410696E+03_r8,     0.1401754E+03_r8,     0.1392868E+03_r8,&  
         0.1384039E+03_r8,     0.1375265E+03_r8,     0.1366548E+03_r8,     0.1357885E+03_r8,     0.1349277E+03_r8,&  
         0.1340724E+03_r8,     0.1332225E+03_r8,     0.1323780E+03_r8,     0.1315389E+03_r8,     0.1307051E+03_r8,&  
         0.1298765E+03_r8,     0.1290533E+03_r8,     0.1282352E+03_r8,     0.1274223E+03_r8,     0.1266146E+03_r8,&  
         0.1258120E+03_r8,     0.1250145E+03_r8,     0.1242220E+03_r8,     0.1234345E+03_r8,     0.1226521E+03_r8,&  
         0.1218746E+03_r8,     0.1211021E+03_r8,     0.1203344E+03_r8,     0.1195716E+03_r8,     0.1188137E+03_r8,&  
         0.1180605E+03_r8,     0.1173121E+03_r8,     0.1165685E+03_r8,     0.1158295E+03_r8,     0.1150953E+03_r8,&  
         0.1143657E+03_r8,     0.1136407E+03_r8,     0.1129204E+03_r8,     0.1122046E+03_r8,     0.1114933E+03_r8,&  
         0.1107866E+03_r8,     0.1100843E+03_r8,     0.1093865E+03_r8,     0.1086931E+03_r8,     0.1080041E+03_r8,&  
         0.1073194E+03_r8,     0.1066391E+03_r8,     0.1059632E+03_r8,     0.1052915E+03_r8,     0.1046240E+03_r8,&  
         0.1039608E+03_r8,     0.1033018E+03_r8,     0.1026470E+03_r8,     0.1019963E+03_r8,     0.1013497E+03_r8,&  
         0.1007073E+03_r8,     0.1000689E+03_r8,     0.9943457E+02_r8,     0.9880426E+02_r8,     0.9817795E+02_r8,&  
         0.9755560E+02_r8,     0.9693719E+02_r8,     0.9632272E+02_r8,     0.9571213E+02_r8,     0.9510540E+02_r8,&  
         0.9450254E+02_r8,     0.9390349E+02_r8,     0.9330823E+02_r8,     0.9271676E+02_r8,     0.9212903E+02_r8,&  
         0.9154501E+02_r8,     0.9096472E+02_r8,     0.9038809E+02_r8,     0.8981513E+02_r8,     0.8924580E+02_r8,&  
         0.8868007E+02_r8,     0.8811793E+02_r8,     0.8755936E+02_r8,     0.8700431E+02_r8,     0.8645279E+02_r8,&  
         0.8590478E+02_r8,     0.8536022E+02_r8,     0.8481912E+02_r8,     0.8428146E+02_r8,     0.8374720E+02_r8/)  

    psaditmk(1:150,    8)= (/ &
         0.2177531E+03_r8,     0.2163797E+03_r8,     0.2150141E+03_r8,     0.2136564E+03_r8,     0.2123066E+03_r8,&  
         0.2109648E+03_r8,     0.2096309E+03_r8,     0.2083051E+03_r8,     0.2069873E+03_r8,     0.2056774E+03_r8,&  
         0.2043755E+03_r8,     0.2030817E+03_r8,     0.2017958E+03_r8,     0.2005179E+03_r8,     0.1992478E+03_r8,&  
         0.1979857E+03_r8,     0.1967314E+03_r8,     0.1954850E+03_r8,     0.1942464E+03_r8,     0.1930155E+03_r8,&  
         0.1917924E+03_r8,     0.1905770E+03_r8,     0.1893693E+03_r8,     0.1881691E+03_r8,     0.1869765E+03_r8,&  
         0.1857914E+03_r8,     0.1846139E+03_r8,     0.1834437E+03_r8,     0.1822810E+03_r8,     0.1811256E+03_r8,&  
         0.1799775E+03_r8,     0.1788367E+03_r8,     0.1777031E+03_r8,     0.1765767E+03_r8,     0.1754574E+03_r8,&  
         0.1743452E+03_r8,     0.1732400E+03_r8,     0.1721419E+03_r8,     0.1710507E+03_r8,     0.1699664E+03_r8,&  
         0.1688891E+03_r8,     0.1678185E+03_r8,     0.1667547E+03_r8,     0.1656976E+03_r8,     0.1646473E+03_r8,&  
         0.1636036E+03_r8,     0.1625665E+03_r8,     0.1615360E+03_r8,     0.1605120E+03_r8,     0.1594946E+03_r8,&  
         0.1584835E+03_r8,     0.1574789E+03_r8,     0.1564807E+03_r8,     0.1554887E+03_r8,     0.1545031E+03_r8,&  
         0.1535237E+03_r8,     0.1525505E+03_r8,     0.1515835E+03_r8,     0.1506226E+03_r8,     0.1496678E+03_r8,&  
         0.1487191E+03_r8,     0.1477764E+03_r8,     0.1468396E+03_r8,     0.1459088E+03_r8,     0.1449839E+03_r8,&  
         0.1440648E+03_r8,     0.1431516E+03_r8,     0.1422442E+03_r8,     0.1413425E+03_r8,     0.1404465E+03_r8,&  
         0.1395562E+03_r8,     0.1386716E+03_r8,     0.1377926E+03_r8,     0.1369191E+03_r8,     0.1360511E+03_r8,&  
         0.1351887E+03_r8,     0.1343318E+03_r8,     0.1334802E+03_r8,     0.1326341E+03_r8,     0.1317934E+03_r8,&  
         0.1309579E+03_r8,     0.1301278E+03_r8,     0.1293029E+03_r8,     0.1284832E+03_r8,     0.1276688E+03_r8,&  
         0.1268595E+03_r8,     0.1260553E+03_r8,     0.1252563E+03_r8,     0.1244623E+03_r8,     0.1236733E+03_r8,&  
         0.1228894E+03_r8,     0.1221104E+03_r8,     0.1213363E+03_r8,     0.1205672E+03_r8,     0.1198029E+03_r8,&  
         0.1190435E+03_r8,     0.1182889E+03_r8,     0.1175390E+03_r8,     0.1167939E+03_r8,     0.1160536E+03_r8,&  
         0.1153179E+03_r8,     0.1145869E+03_r8,     0.1138606E+03_r8,     0.1131388E+03_r8,     0.1124216E+03_r8,&  
         0.1117090E+03_r8,     0.1110009E+03_r8,     0.1102972E+03_r8,     0.1095981E+03_r8,     0.1089033E+03_r8,&  
         0.1082130E+03_r8,     0.1075270E+03_r8,     0.1068454E+03_r8,     0.1061681E+03_r8,     0.1054951E+03_r8,&  
         0.1048264E+03_r8,     0.1041619E+03_r8,     0.1035016E+03_r8,     0.1028455E+03_r8,     0.1021936E+03_r8,&  
         0.1015458E+03_r8,     0.1009021E+03_r8,     0.1002625E+03_r8,     0.9962691E+02_r8,     0.9899539E+02_r8,&  
         0.9836786E+02_r8,     0.9774430E+02_r8,     0.9712471E+02_r8,     0.9650903E+02_r8,     0.9589726E+02_r8,&  
         0.9528938E+02_r8,     0.9468534E+02_r8,     0.9408513E+02_r8,     0.9348873E+02_r8,     0.9289610E+02_r8,&  
         0.9230724E+02_r8,     0.9172211E+02_r8,     0.9114068E+02_r8,     0.9056294E+02_r8,     0.8998887E+02_r8,&  
         0.8941843E+02_r8,     0.8885161E+02_r8,     0.8828838E+02_r8,     0.8772872E+02_r8,     0.8717261E+02_r8,&  
         0.8662003E+02_r8,     0.8607094E+02_r8,     0.8552534E+02_r8,     0.8498320E+02_r8,     0.8444449E+02_r8/)  

    psaditmk(1:150,    9)= (/ &
         0.2195538E+03_r8,     0.2181705E+03_r8,     0.2167950E+03_r8,     0.2154272E+03_r8,     0.2140672E+03_r8,&  
         0.2127152E+03_r8,     0.2113710E+03_r8,     0.2100349E+03_r8,     0.2087067E+03_r8,     0.2073865E+03_r8,&  
         0.2060743E+03_r8,     0.2047701E+03_r8,     0.2034738E+03_r8,     0.2021856E+03_r8,     0.2009053E+03_r8,&  
         0.1996329E+03_r8,     0.1983683E+03_r8,     0.1971118E+03_r8,     0.1958630E+03_r8,     0.1946220E+03_r8,&  
         0.1933888E+03_r8,     0.1921634E+03_r8,     0.1909456E+03_r8,     0.1897355E+03_r8,     0.1885331E+03_r8,&  
         0.1873382E+03_r8,     0.1861509E+03_r8,     0.1849710E+03_r8,     0.1837986E+03_r8,     0.1826336E+03_r8,&  
         0.1814760E+03_r8,     0.1803257E+03_r8,     0.1791827E+03_r8,     0.1780469E+03_r8,     0.1769183E+03_r8,&  
         0.1757969E+03_r8,     0.1746825E+03_r8,     0.1735753E+03_r8,     0.1724750E+03_r8,     0.1713817E+03_r8,&  
         0.1702953E+03_r8,     0.1692158E+03_r8,     0.1681432E+03_r8,     0.1670773E+03_r8,     0.1660182E+03_r8,&  
         0.1649659E+03_r8,     0.1639202E+03_r8,     0.1628811E+03_r8,     0.1618486E+03_r8,     0.1608226E+03_r8,&  
         0.1598031E+03_r8,     0.1587902E+03_r8,     0.1577836E+03_r8,     0.1567834E+03_r8,     0.1557896E+03_r8,&  
         0.1548020E+03_r8,     0.1538207E+03_r8,     0.1528457E+03_r8,     0.1518768E+03_r8,     0.1509140E+03_r8,&  
         0.1499574E+03_r8,     0.1490068E+03_r8,     0.1480623E+03_r8,     0.1471237E+03_r8,     0.1461911E+03_r8,&  
         0.1452644E+03_r8,     0.1443436E+03_r8,     0.1434286E+03_r8,     0.1425194E+03_r8,     0.1416160E+03_r8,&  
         0.1407183E+03_r8,     0.1398262E+03_r8,     0.1389399E+03_r8,     0.1380592E+03_r8,     0.1371840E+03_r8,&  
         0.1363144E+03_r8,     0.1354503E+03_r8,     0.1345917E+03_r8,     0.1337385E+03_r8,     0.1328907E+03_r8,&  
         0.1320484E+03_r8,     0.1312113E+03_r8,     0.1303796E+03_r8,     0.1295531E+03_r8,     0.1287319E+03_r8,&  
         0.1279158E+03_r8,     0.1271050E+03_r8,     0.1262993E+03_r8,     0.1254986E+03_r8,     0.1247031E+03_r8,&  
         0.1239126E+03_r8,     0.1231272E+03_r8,     0.1223466E+03_r8,     0.1215711E+03_r8,     0.1208005E+03_r8,&  
         0.1200347E+03_r8,     0.1192738E+03_r8,     0.1185177E+03_r8,     0.1177664E+03_r8,     0.1170199E+03_r8,&  
         0.1162782E+03_r8,     0.1155411E+03_r8,     0.1148086E+03_r8,     0.1140809E+03_r8,     0.1133577E+03_r8,&  
         0.1126392E+03_r8,     0.1119251E+03_r8,     0.1112157E+03_r8,     0.1105107E+03_r8,     0.1098101E+03_r8,&  
         0.1091141E+03_r8,     0.1084224E+03_r8,     0.1077351E+03_r8,     0.1070522E+03_r8,     0.1063736E+03_r8,&  
         0.1056993E+03_r8,     0.1050292E+03_r8,     0.1043634E+03_r8,     0.1037019E+03_r8,     0.1030445E+03_r8,&  
         0.1023913E+03_r8,     0.1017423E+03_r8,     0.1010973E+03_r8,     0.1004565E+03_r8,     0.9981970E+02_r8,&  
         0.9918694E+02_r8,     0.9855819E+02_r8,     0.9793344E+02_r8,     0.9731264E+02_r8,     0.9669577E+02_r8,&  
         0.9608283E+02_r8,     0.9547376E+02_r8,     0.9486855E+02_r8,     0.9426719E+02_r8,     0.9366962E+02_r8,&  
         0.9307585E+02_r8,     0.9248586E+02_r8,     0.9189958E+02_r8,     0.9131704E+02_r8,     0.9073818E+02_r8,&  
         0.9016299E+02_r8,     0.8959145E+02_r8,     0.8902354E+02_r8,     0.8845921E+02_r8,     0.8789847E+02_r8,&  
         0.8734129E+02_r8,     0.8678763E+02_r8,     0.8623748E+02_r8,     0.8569083E+02_r8,     0.8514764E+02_r8/)  

    psaditmk(1:150,   10)= (/ &
         0.2213667E+03_r8,     0.2199737E+03_r8,     0.2185884E+03_r8,     0.2172107E+03_r8,     0.2158407E+03_r8,&  
         0.2144786E+03_r8,     0.2131243E+03_r8,     0.2117778E+03_r8,     0.2104394E+03_r8,     0.2091089E+03_r8,&  
         0.2077863E+03_r8,     0.2064718E+03_r8,     0.2051652E+03_r8,     0.2038666E+03_r8,     0.2025760E+03_r8,&  
         0.2012933E+03_r8,     0.2000185E+03_r8,     0.1987517E+03_r8,     0.1974927E+03_r8,     0.1962415E+03_r8,&  
         0.1949982E+03_r8,     0.1937627E+03_r8,     0.1925349E+03_r8,     0.1913148E+03_r8,     0.1901024E+03_r8,&  
         0.1888977E+03_r8,     0.1877005E+03_r8,     0.1865108E+03_r8,     0.1853288E+03_r8,     0.1841541E+03_r8,&  
         0.1829869E+03_r8,     0.1818270E+03_r8,     0.1806745E+03_r8,     0.1795293E+03_r8,     0.1783913E+03_r8,&  
         0.1772605E+03_r8,     0.1761369E+03_r8,     0.1750204E+03_r8,     0.1739110E+03_r8,     0.1728086E+03_r8,&  
         0.1717132E+03_r8,     0.1706247E+03_r8,     0.1695432E+03_r8,     0.1684684E+03_r8,     0.1674005E+03_r8,&  
         0.1663394E+03_r8,     0.1652850E+03_r8,     0.1642372E+03_r8,     0.1631962E+03_r8,     0.1621617E+03_r8,&  
         0.1611337E+03_r8,     0.1601123E+03_r8,     0.1590973E+03_r8,     0.1580888E+03_r8,     0.1570867E+03_r8,&  
         0.1560910E+03_r8,     0.1551015E+03_r8,     0.1541183E+03_r8,     0.1531414E+03_r8,     0.1521706E+03_r8,&  
         0.1512060E+03_r8,     0.1502475E+03_r8,     0.1492951E+03_r8,     0.1483487E+03_r8,     0.1474083E+03_r8,&  
         0.1464739E+03_r8,     0.1455454E+03_r8,     0.1446228E+03_r8,     0.1437060E+03_r8,     0.1427951E+03_r8,&  
         0.1418899E+03_r8,     0.1409905E+03_r8,     0.1400968E+03_r8,     0.1392087E+03_r8,     0.1383262E+03_r8,&  
         0.1374494E+03_r8,     0.1365781E+03_r8,     0.1357123E+03_r8,     0.1348521E+03_r8,     0.1339972E+03_r8,&  
         0.1331478E+03_r8,     0.1323038E+03_r8,     0.1314651E+03_r8,     0.1306318E+03_r8,     0.1298037E+03_r8,&  
         0.1289809E+03_r8,     0.1281633E+03_r8,     0.1273509E+03_r8,     0.1265436E+03_r8,     0.1257414E+03_r8,&  
         0.1249444E+03_r8,     0.1241524E+03_r8,     0.1233653E+03_r8,     0.1225833E+03_r8,     0.1218063E+03_r8,&  
         0.1210341E+03_r8,     0.1202669E+03_r8,     0.1195045E+03_r8,     0.1187470E+03_r8,     0.1179943E+03_r8,&  
         0.1172463E+03_r8,     0.1165031E+03_r8,     0.1157646E+03_r8,     0.1150308E+03_r8,     0.1143016E+03_r8,&  
         0.1135770E+03_r8,     0.1128571E+03_r8,     0.1121417E+03_r8,     0.1114308E+03_r8,     0.1107244E+03_r8,&  
         0.1100226E+03_r8,     0.1093251E+03_r8,     0.1086321E+03_r8,     0.1079435E+03_r8,     0.1072592E+03_r8,&  
         0.1065793E+03_r8,     0.1059037E+03_r8,     0.1052324E+03_r8,     0.1045653E+03_r8,     0.1039025E+03_r8,&  
         0.1032439E+03_r8,     0.1025894E+03_r8,     0.1019391E+03_r8,     0.1012929E+03_r8,     0.1006508E+03_r8,&  
         0.1000128E+03_r8,     0.9937881E+02_r8,     0.9874886E+02_r8,     0.9812289E+02_r8,     0.9750088E+02_r8,&  
         0.9688284E+02_r8,     0.9626869E+02_r8,     0.9565844E+02_r8,     0.9505207E+02_r8,     0.9444953E+02_r8,&  
         0.9385081E+02_r8,     0.9325591E+02_r8,     0.9266476E+02_r8,     0.9207736E+02_r8,     0.9149368E+02_r8,&  
         0.9091370E+02_r8,     0.9033741E+02_r8,     0.8976476E+02_r8,     0.8919574E+02_r8,     0.8863033E+02_r8,&  
         0.8806851E+02_r8,     0.8751024E+02_r8,     0.8695551E+02_r8,     0.8640431E+02_r8,     0.8585659E+02_r8/)  

    psaditmk(1:150,   11)= (/ &
         0.2231913E+03_r8,     0.2217890E+03_r8,     0.2203941E+03_r8,     0.2190067E+03_r8,     0.2176269E+03_r8,&  
         0.2162548E+03_r8,     0.2148905E+03_r8,     0.2135339E+03_r8,     0.2121853E+03_r8,     0.2108445E+03_r8,&  
         0.2095117E+03_r8,     0.2081869E+03_r8,     0.2068700E+03_r8,     0.2055610E+03_r8,     0.2042601E+03_r8,&  
         0.2029671E+03_r8,     0.2016820E+03_r8,     0.2004049E+03_r8,     0.1991356E+03_r8,     0.1978743E+03_r8,&  
         0.1966208E+03_r8,     0.1953751E+03_r8,     0.1941373E+03_r8,     0.1929071E+03_r8,     0.1916847E+03_r8,&  
         0.1904700E+03_r8,     0.1892629E+03_r8,     0.1880635E+03_r8,     0.1868716E+03_r8,     0.1856872E+03_r8,&  
         0.1845103E+03_r8,     0.1833408E+03_r8,     0.1821787E+03_r8,     0.1810240E+03_r8,     0.1798766E+03_r8,&  
         0.1787364E+03_r8,     0.1776034E+03_r8,     0.1764777E+03_r8,     0.1753590E+03_r8,     0.1742474E+03_r8,&  
         0.1731429E+03_r8,     0.1720454E+03_r8,     0.1709548E+03_r8,     0.1698711E+03_r8,     0.1687943E+03_r8,&  
         0.1677244E+03_r8,     0.1666612E+03_r8,     0.1656047E+03_r8,     0.1645550E+03_r8,     0.1635118E+03_r8,&  
         0.1624753E+03_r8,     0.1614454E+03_r8,     0.1604220E+03_r8,     0.1594051E+03_r8,     0.1583947E+03_r8,&  
         0.1573906E+03_r8,     0.1563929E+03_r8,     0.1554015E+03_r8,     0.1544164E+03_r8,     0.1534376E+03_r8,&  
         0.1524650E+03_r8,     0.1514985E+03_r8,     0.1505381E+03_r8,     0.1495839E+03_r8,     0.1486357E+03_r8,&  
         0.1476935E+03_r8,     0.1467572E+03_r8,     0.1458270E+03_r8,     0.1449026E+03_r8,     0.1439840E+03_r8,&  
         0.1430713E+03_r8,     0.1421644E+03_r8,     0.1412632E+03_r8,     0.1403678E+03_r8,     0.1394780E+03_r8,&  
         0.1385938E+03_r8,     0.1377153E+03_r8,     0.1368423E+03_r8,     0.1359749E+03_r8,     0.1351129E+03_r8,&  
         0.1342564E+03_r8,     0.1334054E+03_r8,     0.1325597E+03_r8,     0.1317194E+03_r8,     0.1308845E+03_r8,&  
         0.1300548E+03_r8,     0.1292304E+03_r8,     0.1284112E+03_r8,     0.1275972E+03_r8,     0.1267884E+03_r8,&  
         0.1259847E+03_r8,     0.1251861E+03_r8,     0.1243925E+03_r8,     0.1236040E+03_r8,     0.1228205E+03_r8,&  
         0.1220419E+03_r8,     0.1212683E+03_r8,     0.1204996E+03_r8,     0.1197357E+03_r8,     0.1189767E+03_r8,&  
         0.1182225E+03_r8,     0.1174731E+03_r8,     0.1167284E+03_r8,     0.1159885E+03_r8,     0.1152533E+03_r8,&  
         0.1145227E+03_r8,     0.1137967E+03_r8,     0.1130754E+03_r8,     0.1123586E+03_r8,     0.1116463E+03_r8,&  
         0.1109386E+03_r8,     0.1102354E+03_r8,     0.1095366E+03_r8,     0.1088423E+03_r8,     0.1081523E+03_r8,&  
         0.1074667E+03_r8,     0.1067855E+03_r8,     0.1061086E+03_r8,     0.1054360E+03_r8,     0.1047676E+03_r8,&  
         0.1041035E+03_r8,     0.1034436E+03_r8,     0.1027879E+03_r8,     0.1021363E+03_r8,     0.1014889E+03_r8,&  
         0.1008455E+03_r8,     0.1002063E+03_r8,     0.9957105E+02_r8,     0.9893987E+02_r8,     0.9831269E+02_r8,&  
         0.9768951E+02_r8,     0.9707024E+02_r8,     0.9645491E+02_r8,     0.9584350E+02_r8,     0.9523594E+02_r8,&  
         0.9463224E+02_r8,     0.9403238E+02_r8,     0.9343630E+02_r8,     0.9284402E+02_r8,     0.9225548E+02_r8,&  
         0.9167067E+02_r8,     0.9108958E+02_r8,     0.9051216E+02_r8,     0.8993841E+02_r8,     0.8936829E+02_r8,&  
         0.8880179E+02_r8,     0.8823887E+02_r8,     0.8767953E+02_r8,     0.8712373E+02_r8,     0.8657145E+02_r8/)  

    psaditmk(1:150,   12)= (/ &
         0.2250271E+03_r8,     0.2236158E+03_r8,     0.2222117E+03_r8,     0.2208149E+03_r8,     0.2194255E+03_r8,&  
         0.2180437E+03_r8,     0.2166694E+03_r8,     0.2153029E+03_r8,     0.2139442E+03_r8,     0.2125933E+03_r8,&  
         0.2112503E+03_r8,     0.2099152E+03_r8,     0.2085880E+03_r8,     0.2072688E+03_r8,     0.2059575E+03_r8,&  
         0.2046542E+03_r8,     0.2033588E+03_r8,     0.2020714E+03_r8,     0.2007919E+03_r8,     0.1995203E+03_r8,&  
         0.1982566E+03_r8,     0.1970007E+03_r8,     0.1957528E+03_r8,     0.1945125E+03_r8,     0.1932801E+03_r8,&  
         0.1920553E+03_r8,     0.1908383E+03_r8,     0.1896289E+03_r8,     0.1884272E+03_r8,     0.1872330E+03_r8,&  
         0.1860463E+03_r8,     0.1848671E+03_r8,     0.1836954E+03_r8,     0.1825311E+03_r8,     0.1813741E+03_r8,&  
         0.1802244E+03_r8,     0.1790821E+03_r8,     0.1779470E+03_r8,     0.1768190E+03_r8,     0.1756982E+03_r8,&  
         0.1745845E+03_r8,     0.1734778E+03_r8,     0.1723782E+03_r8,     0.1712855E+03_r8,     0.1701997E+03_r8,&  
         0.1691209E+03_r8,     0.1680488E+03_r8,     0.1669835E+03_r8,     0.1659250E+03_r8,     0.1648733E+03_r8,&  
         0.1638281E+03_r8,     0.1627896E+03_r8,     0.1617577E+03_r8,     0.1607323E+03_r8,     0.1597135E+03_r8,&  
         0.1587010E+03_r8,     0.1576950E+03_r8,     0.1566954E+03_r8,     0.1557021E+03_r8,     0.1547151E+03_r8,&  
         0.1537344E+03_r8,     0.1527599E+03_r8,     0.1517915E+03_r8,     0.1508293E+03_r8,     0.1498732E+03_r8,&  
         0.1489232E+03_r8,     0.1479792E+03_r8,     0.1470411E+03_r8,     0.1461090E+03_r8,     0.1451829E+03_r8,&  
         0.1442625E+03_r8,     0.1433481E+03_r8,     0.1424394E+03_r8,     0.1415365E+03_r8,     0.1406393E+03_r8,&  
         0.1397478E+03_r8,     0.1388619E+03_r8,     0.1379817E+03_r8,     0.1371070E+03_r8,     0.1362379E+03_r8,&  
         0.1353743E+03_r8,     0.1345161E+03_r8,     0.1336635E+03_r8,     0.1328161E+03_r8,     0.1319743E+03_r8,&  
         0.1311376E+03_r8,     0.1303064E+03_r8,     0.1294804E+03_r8,     0.1286596E+03_r8,     0.1278440E+03_r8,&  
         0.1270336E+03_r8,     0.1262284E+03_r8,     0.1254282E+03_r8,     0.1246331E+03_r8,     0.1238431E+03_r8,&  
         0.1230580E+03_r8,     0.1222780E+03_r8,     0.1215029E+03_r8,     0.1207326E+03_r8,     0.1199673E+03_r8,&  
         0.1192069E+03_r8,     0.1184512E+03_r8,     0.1177003E+03_r8,     0.1169543E+03_r8,     0.1162129E+03_r8,&  
         0.1154762E+03_r8,     0.1147442E+03_r8,     0.1140168E+03_r8,     0.1132941E+03_r8,     0.1125759E+03_r8,&  
         0.1118623E+03_r8,     0.1111532E+03_r8,     0.1104486E+03_r8,     0.1097485E+03_r8,     0.1090528E+03_r8,&  
         0.1083615E+03_r8,     0.1076746E+03_r8,     0.1069921E+03_r8,     0.1063139E+03_r8,     0.1056399E+03_r8,&  
         0.1049703E+03_r8,     0.1043049E+03_r8,     0.1036437E+03_r8,     0.1029867E+03_r8,     0.1023339E+03_r8,&  
         0.1016852E+03_r8,     0.1010406E+03_r8,     0.1004001E+03_r8,     0.9976366E+02_r8,     0.9913126E+02_r8,&  
         0.9850288E+02_r8,     0.9787846E+02_r8,     0.9725801E+02_r8,     0.9664150E+02_r8,     0.9602888E+02_r8,&  
         0.9542016E+02_r8,     0.9481530E+02_r8,     0.9421426E+02_r8,     0.9361704E+02_r8,     0.9302361E+02_r8,&  
         0.9243394E+02_r8,     0.9184800E+02_r8,     0.9126578E+02_r8,     0.9068724E+02_r8,     0.9011238E+02_r8,&  
         0.8954116E+02_r8,     0.8897356E+02_r8,     0.8840955E+02_r8,     0.8784913E+02_r8,     0.8729225E+02_r8/)  

    psaditmk(1:150,   13)= (/ &
         0.2268734E+03_r8,     0.2254536E+03_r8,     0.2240407E+03_r8,     0.2226348E+03_r8,     0.2212361E+03_r8,&  
         0.2198448E+03_r8,     0.2184609E+03_r8,     0.2170846E+03_r8,     0.2157160E+03_r8,     0.2143551E+03_r8,&  
         0.2130020E+03_r8,     0.2116567E+03_r8,     0.2103193E+03_r8,     0.2089899E+03_r8,     0.2076683E+03_r8,&  
         0.2063547E+03_r8,     0.2050490E+03_r8,     0.2037513E+03_r8,     0.2024615E+03_r8,     0.2011796E+03_r8,&  
         0.1999057E+03_r8,     0.1986396E+03_r8,     0.1973814E+03_r8,     0.1961311E+03_r8,     0.1948885E+03_r8,&  
         0.1936537E+03_r8,     0.1924266E+03_r8,     0.1912073E+03_r8,     0.1899956E+03_r8,     0.1887915E+03_r8,&  
         0.1875950E+03_r8,     0.1864061E+03_r8,     0.1852247E+03_r8,     0.1840507E+03_r8,     0.1828841E+03_r8,&  
         0.1817249E+03_r8,     0.1805731E+03_r8,     0.1794285E+03_r8,     0.1782912E+03_r8,     0.1771610E+03_r8,&  
         0.1760381E+03_r8,     0.1749222E+03_r8,     0.1738134E+03_r8,     0.1727116E+03_r8,     0.1716168E+03_r8,&  
         0.1705290E+03_r8,     0.1694480E+03_r8,     0.1683739E+03_r8,     0.1673066E+03_r8,     0.1662460E+03_r8,&  
         0.1651922E+03_r8,     0.1641451E+03_r8,     0.1631046E+03_r8,     0.1620706E+03_r8,     0.1610433E+03_r8,&  
         0.1600224E+03_r8,     0.1590080E+03_r8,     0.1580001E+03_r8,     0.1569986E+03_r8,     0.1560033E+03_r8,&  
         0.1550144E+03_r8,     0.1540318E+03_r8,     0.1530554E+03_r8,     0.1520852E+03_r8,     0.1511211E+03_r8,&  
         0.1501632E+03_r8,     0.1492113E+03_r8,     0.1482654E+03_r8,     0.1473256E+03_r8,     0.1463917E+03_r8,&  
         0.1454637E+03_r8,     0.1445416E+03_r8,     0.1436254E+03_r8,     0.1427149E+03_r8,     0.1418103E+03_r8,&  
         0.1409113E+03_r8,     0.1400181E+03_r8,     0.1391305E+03_r8,     0.1382486E+03_r8,     0.1373722E+03_r8,&  
         0.1365014E+03_r8,     0.1356362E+03_r8,     0.1347764E+03_r8,     0.1339220E+03_r8,     0.1330731E+03_r8,&  
         0.1322295E+03_r8,     0.1313913E+03_r8,     0.1305585E+03_r8,     0.1297308E+03_r8,     0.1289085E+03_r8,&  
         0.1280914E+03_r8,     0.1272794E+03_r8,     0.1264725E+03_r8,     0.1256708E+03_r8,     0.1248742E+03_r8,&  
         0.1240826E+03_r8,     0.1232961E+03_r8,     0.1225145E+03_r8,     0.1217379E+03_r8,     0.1209662E+03_r8,&  
         0.1201994E+03_r8,     0.1194375E+03_r8,     0.1186804E+03_r8,     0.1179281E+03_r8,     0.1171805E+03_r8,&  
         0.1164377E+03_r8,     0.1156996E+03_r8,     0.1149662E+03_r8,     0.1142374E+03_r8,     0.1135133E+03_r8,&  
         0.1127937E+03_r8,     0.1120787E+03_r8,     0.1113682E+03_r8,     0.1106623E+03_r8,     0.1099608E+03_r8,&  
         0.1092638E+03_r8,     0.1085711E+03_r8,     0.1078829E+03_r8,     0.1071991E+03_r8,     0.1065195E+03_r8,&  
         0.1058443E+03_r8,     0.1051733E+03_r8,     0.1045066E+03_r8,     0.1038442E+03_r8,     0.1031859E+03_r8,&  
         0.1025318E+03_r8,     0.1018819E+03_r8,     0.1012361E+03_r8,     0.1005943E+03_r8,     0.9995665E+02_r8,&  
         0.9932304E+02_r8,     0.9869343E+02_r8,     0.9806780E+02_r8,     0.9744617E+02_r8,     0.9682845E+02_r8,&  
         0.9621465E+02_r8,     0.9560476E+02_r8,     0.9499872E+02_r8,     0.9439652E+02_r8,     0.9379815E+02_r8,&  
         0.9320356E+02_r8,     0.9261275E+02_r8,     0.9202568E+02_r8,     0.9144233E+02_r8,     0.9086268E+02_r8,&  
         0.9028671E+02_r8,     0.8971437E+02_r8,     0.8914568E+02_r8,     0.8858059E+02_r8,     0.8801907E+02_r8/)  

    psaditmk(1:150,   14)= (/ &
         0.2287293E+03_r8,     0.2273016E+03_r8,     0.2258803E+03_r8,     0.2244658E+03_r8,     0.2230582E+03_r8,&  
         0.2216577E+03_r8,     0.2202645E+03_r8,     0.2188786E+03_r8,     0.2175004E+03_r8,     0.2161296E+03_r8,&  
         0.2147665E+03_r8,     0.2134112E+03_r8,     0.2120637E+03_r8,     0.2107241E+03_r8,     0.2093923E+03_r8,&  
         0.2080685E+03_r8,     0.2067525E+03_r8,     0.2054445E+03_r8,     0.2041444E+03_r8,     0.2028523E+03_r8,&  
         0.2015681E+03_r8,     0.2002918E+03_r8,     0.1990234E+03_r8,     0.1977628E+03_r8,     0.1965101E+03_r8,&  
         0.1952652E+03_r8,     0.1940280E+03_r8,     0.1927986E+03_r8,     0.1915770E+03_r8,     0.1903630E+03_r8,&  
         0.1891566E+03_r8,     0.1879578E+03_r8,     0.1867666E+03_r8,     0.1855829E+03_r8,     0.1844067E+03_r8,&  
         0.1832378E+03_r8,     0.1820764E+03_r8,     0.1809223E+03_r8,     0.1797756E+03_r8,     0.1786360E+03_r8,&  
         0.1775037E+03_r8,     0.1763786E+03_r8,     0.1752606E+03_r8,     0.1741496E+03_r8,     0.1730457E+03_r8,&  
         0.1719488E+03_r8,     0.1708588E+03_r8,     0.1697758E+03_r8,     0.1686996E+03_r8,     0.1676302E+03_r8,&  
         0.1665676E+03_r8,     0.1655117E+03_r8,     0.1644626E+03_r8,     0.1634201E+03_r8,     0.1623842E+03_r8,&  
         0.1613548E+03_r8,     0.1603320E+03_r8,     0.1593156E+03_r8,     0.1583057E+03_r8,     0.1573022E+03_r8,&  
         0.1563051E+03_r8,     0.1553143E+03_r8,     0.1543298E+03_r8,     0.1533515E+03_r8,     0.1523794E+03_r8,&  
         0.1514134E+03_r8,     0.1504536E+03_r8,     0.1494999E+03_r8,     0.1485522E+03_r8,     0.1476106E+03_r8,&  
         0.1466749E+03_r8,     0.1457451E+03_r8,     0.1448212E+03_r8,     0.1439032E+03_r8,     0.1429910E+03_r8,&  
         0.1420846E+03_r8,     0.1411839E+03_r8,     0.1402890E+03_r8,     0.1393997E+03_r8,     0.1385160E+03_r8,&  
         0.1376380E+03_r8,     0.1367655E+03_r8,     0.1358985E+03_r8,     0.1350371E+03_r8,     0.1341811E+03_r8,&  
         0.1333305E+03_r8,     0.1324853E+03_r8,     0.1316455E+03_r8,     0.1308110E+03_r8,     0.1299818E+03_r8,&  
         0.1291579E+03_r8,     0.1283391E+03_r8,     0.1275256E+03_r8,     0.1267172E+03_r8,     0.1259140E+03_r8,&  
         0.1251158E+03_r8,     0.1243227E+03_r8,     0.1235346E+03_r8,     0.1227515E+03_r8,     0.1219734E+03_r8,&  
         0.1212002E+03_r8,     0.1204319E+03_r8,     0.1196685E+03_r8,     0.1189099E+03_r8,     0.1181562E+03_r8,&  
         0.1174072E+03_r8,     0.1166629E+03_r8,     0.1159234E+03_r8,     0.1151886E+03_r8,     0.1144584E+03_r8,&  
         0.1137328E+03_r8,     0.1130119E+03_r8,     0.1122955E+03_r8,     0.1115837E+03_r8,     0.1108764E+03_r8,&  
         0.1101735E+03_r8,     0.1094751E+03_r8,     0.1087812E+03_r8,     0.1080916E+03_r8,     0.1074064E+03_r8,&  
         0.1067256E+03_r8,     0.1060490E+03_r8,     0.1053768E+03_r8,     0.1047088E+03_r8,     0.1040451E+03_r8,&  
         0.1033855E+03_r8,     0.1027302E+03_r8,     0.1020790E+03_r8,     0.1014319E+03_r8,     0.1007889E+03_r8,&  
         0.1001500E+03_r8,     0.9951516E+02_r8,     0.9888433E+02_r8,     0.9825751E+02_r8,     0.9763466E+02_r8,&  
         0.9701575E+02_r8,     0.9640078E+02_r8,     0.9578969E+02_r8,     0.9518249E+02_r8,     0.9457912E+02_r8,&  
         0.9397959E+02_r8,     0.9338386E+02_r8,     0.9279190E+02_r8,     0.9220369E+02_r8,     0.9161922E+02_r8,&  
         0.9103844E+02_r8,     0.9046135E+02_r8,     0.8988792E+02_r8,     0.8931812E+02_r8,     0.8875193E+02_r8/)  

    psaditmk(1:150,   15)= (/ &
         0.2305941E+03_r8,     0.2291589E+03_r8,     0.2277299E+03_r8,     0.2263072E+03_r8,     0.2248912E+03_r8,&  
         0.2234819E+03_r8,     0.2220797E+03_r8,     0.2206846E+03_r8,     0.2192969E+03_r8,     0.2179165E+03_r8,&  
         0.2165437E+03_r8,     0.2151785E+03_r8,     0.2138211E+03_r8,     0.2124713E+03_r8,     0.2111295E+03_r8,&  
         0.2097954E+03_r8,     0.2084693E+03_r8,     0.2071510E+03_r8,     0.2058407E+03_r8,     0.2045383E+03_r8,&  
         0.2032438E+03_r8,     0.2019572E+03_r8,     0.2006786E+03_r8,     0.1994078E+03_r8,     0.1981449E+03_r8,&  
         0.1968898E+03_r8,     0.1956426E+03_r8,     0.1944031E+03_r8,     0.1931714E+03_r8,     0.1919473E+03_r8,&  
         0.1907310E+03_r8,     0.1895223E+03_r8,     0.1883213E+03_r8,     0.1871278E+03_r8,     0.1859418E+03_r8,&  
         0.1847633E+03_r8,     0.1835922E+03_r8,     0.1824286E+03_r8,     0.1812723E+03_r8,     0.1801233E+03_r8,&  
         0.1789816E+03_r8,     0.1778470E+03_r8,     0.1767197E+03_r8,     0.1755995E+03_r8,     0.1744865E+03_r8,&  
         0.1733804E+03_r8,     0.1722814E+03_r8,     0.1711893E+03_r8,     0.1701042E+03_r8,     0.1690259E+03_r8,&  
         0.1679545E+03_r8,     0.1668898E+03_r8,     0.1658319E+03_r8,     0.1647807E+03_r8,     0.1637362E+03_r8,&  
         0.1626983E+03_r8,     0.1616669E+03_r8,     0.1606421E+03_r8,     0.1596238E+03_r8,     0.1586119E+03_r8,&  
         0.1576065E+03_r8,     0.1566075E+03_r8,     0.1556147E+03_r8,     0.1546283E+03_r8,     0.1536481E+03_r8,&  
         0.1526741E+03_r8,     0.1517063E+03_r8,     0.1507447E+03_r8,     0.1497891E+03_r8,     0.1488396E+03_r8,&  
         0.1478961E+03_r8,     0.1469586E+03_r8,     0.1460271E+03_r8,     0.1451014E+03_r8,     0.1441816E+03_r8,&  
         0.1432676E+03_r8,     0.1423595E+03_r8,     0.1414570E+03_r8,     0.1405603E+03_r8,     0.1396693E+03_r8,&  
         0.1387840E+03_r8,     0.1379042E+03_r8,     0.1370301E+03_r8,     0.1361614E+03_r8,     0.1352983E+03_r8,&  
         0.1344406E+03_r8,     0.1335884E+03_r8,     0.1327416E+03_r8,     0.1319002E+03_r8,     0.1310640E+03_r8,&  
         0.1302332E+03_r8,     0.1294077E+03_r8,     0.1285874E+03_r8,     0.1277723E+03_r8,     0.1269623E+03_r8,&  
         0.1261575E+03_r8,     0.1253578E+03_r8,     0.1245632E+03_r8,     0.1237735E+03_r8,     0.1229890E+03_r8,&  
         0.1222094E+03_r8,     0.1214347E+03_r8,     0.1206649E+03_r8,     0.1199000E+03_r8,     0.1191400E+03_r8,&  
         0.1183847E+03_r8,     0.1176343E+03_r8,     0.1168886E+03_r8,     0.1161477E+03_r8,     0.1154114E+03_r8,&  
         0.1146798E+03_r8,     0.1139529E+03_r8,     0.1132305E+03_r8,     0.1125127E+03_r8,     0.1117995E+03_r8,&  
         0.1110908E+03_r8,     0.1103866E+03_r8,     0.1096869E+03_r8,     0.1089916E+03_r8,     0.1083007E+03_r8,&  
         0.1076142E+03_r8,     0.1069320E+03_r8,     0.1062542E+03_r8,     0.1055806E+03_r8,     0.1049114E+03_r8,&  
         0.1042463E+03_r8,     0.1035855E+03_r8,     0.1029289E+03_r8,     0.1022764E+03_r8,     0.1016281E+03_r8,&  
         0.1009839E+03_r8,     0.1003437E+03_r8,     0.9970766E+02_r8,     0.9907562E+02_r8,     0.9844758E+02_r8,&  
         0.9782352E+02_r8,     0.9720343E+02_r8,     0.9658725E+02_r8,     0.9597499E+02_r8,     0.9536661E+02_r8,&  
         0.9476208E+02_r8,     0.9416138E+02_r8,     0.9356451E+02_r8,     0.9297140E+02_r8,     0.9238205E+02_r8,&  
         0.9179645E+02_r8,     0.9121455E+02_r8,     0.9063634E+02_r8,     0.9006180E+02_r8,     0.8949090E+02_r8/)  

    psaditmk(1:150,   16)= (/ &
         0.2324663E+03_r8,     0.2310246E+03_r8,     0.2295885E+03_r8,     0.2281583E+03_r8,     0.2267343E+03_r8,&  
         0.2253168E+03_r8,     0.2239060E+03_r8,     0.2225020E+03_r8,     0.2211052E+03_r8,     0.2197155E+03_r8,&  
         0.2183332E+03_r8,     0.2169584E+03_r8,     0.2155911E+03_r8,     0.2142315E+03_r8,     0.2128796E+03_r8,&  
         0.2115355E+03_r8,     0.2101992E+03_r8,     0.2088708E+03_r8,     0.2075502E+03_r8,     0.2062375E+03_r8,&  
         0.2049328E+03_r8,     0.2036360E+03_r8,     0.2023471E+03_r8,     0.2010661E+03_r8,     0.1997929E+03_r8,&  
         0.1985277E+03_r8,     0.1972702E+03_r8,     0.1960206E+03_r8,     0.1947788E+03_r8,     0.1935448E+03_r8,&  
         0.1923184E+03_r8,     0.1910998E+03_r8,     0.1898888E+03_r8,     0.1886854E+03_r8,     0.1874897E+03_r8,&  
         0.1863014E+03_r8,     0.1851206E+03_r8,     0.1839473E+03_r8,     0.1827814E+03_r8,     0.1816229E+03_r8,&  
         0.1804717E+03_r8,     0.1793278E+03_r8,     0.1781911E+03_r8,     0.1770616E+03_r8,     0.1759392E+03_r8,&  
         0.1748240E+03_r8,     0.1737158E+03_r8,     0.1726147E+03_r8,     0.1715205E+03_r8,     0.1704332E+03_r8,&  
         0.1693529E+03_r8,     0.1682794E+03_r8,     0.1672126E+03_r8,     0.1661527E+03_r8,     0.1650995E+03_r8,&  
         0.1640529E+03_r8,     0.1630130E+03_r8,     0.1619797E+03_r8,     0.1609529E+03_r8,     0.1599326E+03_r8,&  
         0.1589188E+03_r8,     0.1579114E+03_r8,     0.1569104E+03_r8,     0.1559158E+03_r8,     0.1549274E+03_r8,&  
         0.1539453E+03_r8,     0.1529695E+03_r8,     0.1519998E+03_r8,     0.1510363E+03_r8,     0.1500789E+03_r8,&  
         0.1491275E+03_r8,     0.1481822E+03_r8,     0.1472429E+03_r8,     0.1463095E+03_r8,     0.1453820E+03_r8,&  
         0.1444605E+03_r8,     0.1435448E+03_r8,     0.1426348E+03_r8,     0.1417307E+03_r8,     0.1408322E+03_r8,&  
         0.1399395E+03_r8,     0.1390524E+03_r8,     0.1381710E+03_r8,     0.1372951E+03_r8,     0.1364248E+03_r8,&  
         0.1355600E+03_r8,     0.1347007E+03_r8,     0.1338468E+03_r8,     0.1329984E+03_r8,     0.1321553E+03_r8,&  
         0.1313176E+03_r8,     0.1304852E+03_r8,     0.1296580E+03_r8,     0.1288361E+03_r8,     0.1280195E+03_r8,&  
         0.1272079E+03_r8,     0.1264016E+03_r8,     0.1256003E+03_r8,     0.1248041E+03_r8,     0.1240130E+03_r8,&  
         0.1232269E+03_r8,     0.1224457E+03_r8,     0.1216696E+03_r8,     0.1208983E+03_r8,     0.1201319E+03_r8,&  
         0.1193704E+03_r8,     0.1186137E+03_r8,     0.1178618E+03_r8,     0.1171147E+03_r8,     0.1163723E+03_r8,&  
         0.1156347E+03_r8,     0.1149016E+03_r8,     0.1141733E+03_r8,     0.1134495E+03_r8,     0.1127304E+03_r8,&  
         0.1120158E+03_r8,     0.1113057E+03_r8,     0.1106002E+03_r8,     0.1098991E+03_r8,     0.1092024E+03_r8,&  
         0.1085102E+03_r8,     0.1078224E+03_r8,     0.1071389E+03_r8,     0.1064597E+03_r8,     0.1057849E+03_r8,&  
         0.1051143E+03_r8,     0.1044480E+03_r8,     0.1037859E+03_r8,     0.1031280E+03_r8,     0.1024743E+03_r8,&  
         0.1018247E+03_r8,     0.1011792E+03_r8,     0.1005378E+03_r8,     0.9990055E+02_r8,     0.9926728E+02_r8,&  
         0.9863802E+02_r8,     0.9801276E+02_r8,     0.9739146E+02_r8,     0.9677410E+02_r8,     0.9616065E+02_r8,&  
         0.9555109E+02_r8,     0.9494540E+02_r8,     0.9434354E+02_r8,     0.9374550E+02_r8,     0.9315125E+02_r8,&  
         0.9256077E+02_r8,     0.9197402E+02_r8,     0.9139100E+02_r8,     0.9081168E+02_r8,     0.9023602E+02_r8/)  

    psaditmk(1:150,   17)= (/ &
         0.2343449E+03_r8,     0.2328973E+03_r8,     0.2314549E+03_r8,     0.2300179E+03_r8,     0.2285866E+03_r8,&  
         0.2271615E+03_r8,     0.2257426E+03_r8,     0.2243302E+03_r8,     0.2229247E+03_r8,     0.2215261E+03_r8,&  
         0.2201346E+03_r8,     0.2187504E+03_r8,     0.2173735E+03_r8,     0.2160042E+03_r8,     0.2146425E+03_r8,&  
         0.2132884E+03_r8,     0.2119421E+03_r8,     0.2106036E+03_r8,     0.2092729E+03_r8,     0.2079500E+03_r8,&  
         0.2066351E+03_r8,     0.2053280E+03_r8,     0.2040289E+03_r8,     0.2027376E+03_r8,     0.2014542E+03_r8,&  
         0.2001788E+03_r8,     0.1989111E+03_r8,     0.1976514E+03_r8,     0.1963994E+03_r8,     0.1951553E+03_r8,&  
         0.1939189E+03_r8,     0.1926902E+03_r8,     0.1914692E+03_r8,     0.1902560E+03_r8,     0.1890503E+03_r8,&  
         0.1878522E+03_r8,     0.1866617E+03_r8,     0.1854787E+03_r8,     0.1843031E+03_r8,     0.1831349E+03_r8,&  
         0.1819742E+03_r8,     0.1808207E+03_r8,     0.1796746E+03_r8,     0.1785357E+03_r8,     0.1774041E+03_r8,&  
         0.1762796E+03_r8,     0.1751622E+03_r8,     0.1740518E+03_r8,     0.1729486E+03_r8,     0.1718523E+03_r8,&  
         0.1707629E+03_r8,     0.1696805E+03_r8,     0.1686049E+03_r8,     0.1675361E+03_r8,     0.1664741E+03_r8,&  
         0.1654189E+03_r8,     0.1643703E+03_r8,     0.1633283E+03_r8,     0.1622930E+03_r8,     0.1612642E+03_r8,&  
         0.1602420E+03_r8,     0.1592262E+03_r8,     0.1582169E+03_r8,     0.1572140E+03_r8,     0.1562174E+03_r8,&  
         0.1552271E+03_r8,     0.1542431E+03_r8,     0.1532654E+03_r8,     0.1522938E+03_r8,     0.1513285E+03_r8,&  
         0.1503692E+03_r8,     0.1494160E+03_r8,     0.1484689E+03_r8,     0.1475277E+03_r8,     0.1465925E+03_r8,&  
         0.1456633E+03_r8,     0.1447399E+03_r8,     0.1438224E+03_r8,     0.1429108E+03_r8,     0.1420048E+03_r8,&  
         0.1411047E+03_r8,     0.1402102E+03_r8,     0.1393214E+03_r8,     0.1384383E+03_r8,     0.1375607E+03_r8,&  
         0.1366887E+03_r8,     0.1358223E+03_r8,     0.1349613E+03_r8,     0.1341058E+03_r8,     0.1332557E+03_r8,&  
         0.1324110E+03_r8,     0.1315716E+03_r8,     0.1307376E+03_r8,     0.1299088E+03_r8,     0.1290854E+03_r8,&  
         0.1282671E+03_r8,     0.1274540E+03_r8,     0.1266461E+03_r8,     0.1258433E+03_r8,     0.1250456E+03_r8,&  
         0.1242529E+03_r8,     0.1234653E+03_r8,     0.1226826E+03_r8,     0.1219049E+03_r8,     0.1211322E+03_r8,&  
         0.1203643E+03_r8,     0.1196013E+03_r8,     0.1188432E+03_r8,     0.1180899E+03_r8,     0.1173413E+03_r8,&  
         0.1165975E+03_r8,     0.1158583E+03_r8,     0.1151239E+03_r8,     0.1143941E+03_r8,     0.1136690E+03_r8,&  
         0.1129485E+03_r8,     0.1122325E+03_r8,     0.1115210E+03_r8,     0.1108141E+03_r8,     0.1101117E+03_r8,&  
         0.1094137E+03_r8,     0.1087201E+03_r8,     0.1080309E+03_r8,     0.1073461E+03_r8,     0.1066657E+03_r8,&  
         0.1059895E+03_r8,     0.1053176E+03_r8,     0.1046500E+03_r8,     0.1039867E+03_r8,     0.1033275E+03_r8,&  
         0.1026725E+03_r8,     0.1020217E+03_r8,     0.1013749E+03_r8,     0.1007323E+03_r8,     0.1000938E+03_r8,&  
         0.9945930E+02_r8,     0.9882884E+02_r8,     0.9820236E+02_r8,     0.9757986E+02_r8,     0.9696130E+02_r8,&  
         0.9634667E+02_r8,     0.9573593E+02_r8,     0.9512907E+02_r8,     0.9452604E+02_r8,     0.9392684E+02_r8,&  
         0.9333144E+02_r8,     0.9273981E+02_r8,     0.9215194E+02_r8,     0.9156779E+02_r8,     0.9098734E+02_r8/)  

    psaditmk(1:150,   18)= (/ &
         0.2362281E+03_r8,     0.2347758E+03_r8,     0.2333280E+03_r8,     0.2318850E+03_r8,     0.2304472E+03_r8,&  
         0.2290150E+03_r8,     0.2275886E+03_r8,     0.2261685E+03_r8,     0.2247547E+03_r8,     0.2233476E+03_r8,&  
         0.2219473E+03_r8,     0.2205540E+03_r8,     0.2191679E+03_r8,     0.2177891E+03_r8,     0.2164178E+03_r8,&  
         0.2150540E+03_r8,     0.2136978E+03_r8,     0.2123493E+03_r8,     0.2110085E+03_r8,     0.2096756E+03_r8,&  
         0.2083505E+03_r8,     0.2070332E+03_r8,     0.2057239E+03_r8,     0.2044224E+03_r8,     0.2031288E+03_r8,&  
         0.2018431E+03_r8,     0.2005653E+03_r8,     0.1992953E+03_r8,     0.1980332E+03_r8,     0.1967789E+03_r8,&  
         0.1955324E+03_r8,     0.1942937E+03_r8,     0.1930627E+03_r8,     0.1918394E+03_r8,     0.1906238E+03_r8,&  
         0.1894158E+03_r8,     0.1882155E+03_r8,     0.1870227E+03_r8,     0.1858374E+03_r8,     0.1846595E+03_r8,&  
         0.1834891E+03_r8,     0.1823261E+03_r8,     0.1811705E+03_r8,     0.1800222E+03_r8,     0.1788811E+03_r8,&  
         0.1777472E+03_r8,     0.1766205E+03_r8,     0.1755010E+03_r8,     0.1743885E+03_r8,     0.1732831E+03_r8,&  
         0.1721847E+03_r8,     0.1710932E+03_r8,     0.1700087E+03_r8,     0.1689310E+03_r8,     0.1678602E+03_r8,&  
         0.1667961E+03_r8,     0.1657388E+03_r8,     0.1646882E+03_r8,     0.1636443E+03_r8,     0.1626069E+03_r8,&  
         0.1615762E+03_r8,     0.1605520E+03_r8,     0.1595342E+03_r8,     0.1585229E+03_r8,     0.1575181E+03_r8,&  
         0.1565196E+03_r8,     0.1555274E+03_r8,     0.1545415E+03_r8,     0.1535619E+03_r8,     0.1525885E+03_r8,&  
         0.1516212E+03_r8,     0.1506601E+03_r8,     0.1497051E+03_r8,     0.1487561E+03_r8,     0.1478131E+03_r8,&  
         0.1468761E+03_r8,     0.1459451E+03_r8,     0.1450199E+03_r8,     0.1441007E+03_r8,     0.1431872E+03_r8,&  
         0.1422795E+03_r8,     0.1413777E+03_r8,     0.1404815E+03_r8,     0.1395909E+03_r8,     0.1387061E+03_r8,&  
         0.1378268E+03_r8,     0.1369531E+03_r8,     0.1360850E+03_r8,     0.1352224E+03_r8,     0.1343652E+03_r8,&  
         0.1335135E+03_r8,     0.1326671E+03_r8,     0.1318261E+03_r8,     0.1309905E+03_r8,     0.1301602E+03_r8,&  
         0.1293351E+03_r8,     0.1285152E+03_r8,     0.1277006E+03_r8,     0.1268911E+03_r8,     0.1260867E+03_r8,&  
         0.1252875E+03_r8,     0.1244932E+03_r8,     0.1237041E+03_r8,     0.1229200E+03_r8,     0.1221408E+03_r8,&  
         0.1213665E+03_r8,     0.1205972E+03_r8,     0.1198327E+03_r8,     0.1190731E+03_r8,     0.1183183E+03_r8,&  
         0.1175683E+03_r8,     0.1168230E+03_r8,     0.1160825E+03_r8,     0.1153466E+03_r8,     0.1146154E+03_r8,&  
         0.1138889E+03_r8,     0.1131670E+03_r8,     0.1124496E+03_r8,     0.1117368E+03_r8,     0.1110285E+03_r8,&  
         0.1103247E+03_r8,     0.1096253E+03_r8,     0.1089304E+03_r8,     0.1082399E+03_r8,     0.1075538E+03_r8,&  
         0.1068720E+03_r8,     0.1061945E+03_r8,     0.1055214E+03_r8,     0.1048525E+03_r8,     0.1041878E+03_r8,&  
         0.1035274E+03_r8,     0.1028711E+03_r8,     0.1022190E+03_r8,     0.1015711E+03_r8,     0.1009272E+03_r8,&  
         0.1002874E+03_r8,     0.9965171E+02_r8,     0.9902001E+02_r8,     0.9839233E+02_r8,     0.9776862E+02_r8,&  
         0.9714887E+02_r8,     0.9653305E+02_r8,     0.9592113E+02_r8,     0.9531309E+02_r8,     0.9470889E+02_r8,&  
         0.9410854E+02_r8,     0.9351199E+02_r8,     0.9291921E+02_r8,     0.9233021E+02_r8,     0.9174492E+02_r8/)  

    psaditmk(1:150,   19)= (/ &
         0.2381144E+03_r8,     0.2366584E+03_r8,     0.2352062E+03_r8,     0.2337581E+03_r8,     0.2323147E+03_r8,&  
         0.2308763E+03_r8,     0.2294432E+03_r8,     0.2280158E+03_r8,     0.2265945E+03_r8,     0.2251793E+03_r8,&  
         0.2237707E+03_r8,     0.2223688E+03_r8,     0.2209738E+03_r8,     0.2195859E+03_r8,     0.2182052E+03_r8,&  
         0.2168319E+03_r8,     0.2154660E+03_r8,     0.2141077E+03_r8,     0.2127571E+03_r8,     0.2114141E+03_r8,&  
         0.2100790E+03_r8,     0.2087516E+03_r8,     0.2074321E+03_r8,     0.2061204E+03_r8,     0.2048166E+03_r8,&  
         0.2035207E+03_r8,     0.2022326E+03_r8,     0.2009524E+03_r8,     0.1996802E+03_r8,     0.1984157E+03_r8,&  
         0.1971591E+03_r8,     0.1959102E+03_r8,     0.1946692E+03_r8,     0.1934359E+03_r8,     0.1922103E+03_r8,&  
         0.1909924E+03_r8,     0.1897820E+03_r8,     0.1885794E+03_r8,     0.1873843E+03_r8,     0.1861967E+03_r8,&  
         0.1850166E+03_r8,     0.1838440E+03_r8,     0.1826788E+03_r8,     0.1815209E+03_r8,     0.1803703E+03_r8,&  
         0.1792271E+03_r8,     0.1780910E+03_r8,     0.1769622E+03_r8,     0.1758405E+03_r8,     0.1747259E+03_r8,&  
         0.1736183E+03_r8,     0.1725178E+03_r8,     0.1714242E+03_r8,     0.1703376E+03_r8,     0.1692578E+03_r8,&  
         0.1681849E+03_r8,     0.1671188E+03_r8,     0.1660594E+03_r8,     0.1650068E+03_r8,     0.1639608E+03_r8,&  
         0.1629215E+03_r8,     0.1618887E+03_r8,     0.1608625E+03_r8,     0.1598428E+03_r8,     0.1588296E+03_r8,&  
         0.1578228E+03_r8,     0.1568223E+03_r8,     0.1558282E+03_r8,     0.1548405E+03_r8,     0.1538589E+03_r8,&  
         0.1528836E+03_r8,     0.1519145E+03_r8,     0.1509515E+03_r8,     0.1499946E+03_r8,     0.1490438E+03_r8,&  
         0.1480990E+03_r8,     0.1471602E+03_r8,     0.1462274E+03_r8,     0.1453005E+03_r8,     0.1443794E+03_r8,&  
         0.1434642E+03_r8,     0.1425548E+03_r8,     0.1416511E+03_r8,     0.1407532E+03_r8,     0.1398610E+03_r8,&  
         0.1389744E+03_r8,     0.1380934E+03_r8,     0.1372181E+03_r8,     0.1363482E+03_r8,     0.1354839E+03_r8,&  
         0.1346251E+03_r8,     0.1337717E+03_r8,     0.1329237E+03_r8,     0.1320811E+03_r8,     0.1312439E+03_r8,&  
         0.1304119E+03_r8,     0.1295853E+03_r8,     0.1287638E+03_r8,     0.1279476E+03_r8,     0.1271365E+03_r8,&  
         0.1263306E+03_r8,     0.1255298E+03_r8,     0.1247341E+03_r8,     0.1239434E+03_r8,     0.1231577E+03_r8,&  
         0.1223770E+03_r8,     0.1216013E+03_r8,     0.1208305E+03_r8,     0.1200645E+03_r8,     0.1193034E+03_r8,&  
         0.1185472E+03_r8,     0.1177957E+03_r8,     0.1170490E+03_r8,     0.1163070E+03_r8,     0.1155697E+03_r8,&  
         0.1148372E+03_r8,     0.1141092E+03_r8,     0.1133859E+03_r8,     0.1126671E+03_r8,     0.1119529E+03_r8,&  
         0.1112433E+03_r8,     0.1105381E+03_r8,     0.1098374E+03_r8,     0.1091411E+03_r8,     0.1084493E+03_r8,&  
         0.1077618E+03_r8,     0.1070787E+03_r8,     0.1064000E+03_r8,     0.1057255E+03_r8,     0.1050553E+03_r8,&  
         0.1043894E+03_r8,     0.1037276E+03_r8,     0.1030701E+03_r8,     0.1024168E+03_r8,     0.1017675E+03_r8,&  
         0.1011224E+03_r8,     0.1004814E+03_r8,     0.9984447E+02_r8,     0.9921156E+02_r8,     0.9858266E+02_r8,&  
         0.9795775E+02_r8,     0.9733680E+02_r8,     0.9671979E+02_r8,     0.9610667E+02_r8,     0.9549746E+02_r8,&  
         0.9489211E+02_r8,     0.9429058E+02_r8,     0.9369287E+02_r8,     0.9309896E+02_r8,     0.9250880E+02_r8/)  

    psaditmk(1:150,   20)= (/ &
         0.2400018E+03_r8,     0.2385434E+03_r8,     0.2370879E+03_r8,     0.2356359E+03_r8,     0.2341878E+03_r8,&  
         0.2327440E+03_r8,     0.2313051E+03_r8,     0.2298713E+03_r8,     0.2284430E+03_r8,     0.2270205E+03_r8,&  
         0.2256040E+03_r8,     0.2241940E+03_r8,     0.2227906E+03_r8,     0.2213939E+03_r8,     0.2200042E+03_r8,&  
         0.2186217E+03_r8,     0.2172464E+03_r8,     0.2158786E+03_r8,     0.2145182E+03_r8,     0.2131654E+03_r8,&  
         0.2118203E+03_r8,     0.2104829E+03_r8,     0.2091533E+03_r8,     0.2078315E+03_r8,     0.2065176E+03_r8,&  
         0.2052115E+03_r8,     0.2039132E+03_r8,     0.2026228E+03_r8,     0.2013404E+03_r8,     0.2000657E+03_r8,&  
         0.1987989E+03_r8,     0.1975399E+03_r8,     0.1962888E+03_r8,     0.1950454E+03_r8,     0.1938098E+03_r8,&  
         0.1925818E+03_r8,     0.1913616E+03_r8,     0.1901490E+03_r8,     0.1889440E+03_r8,     0.1877466E+03_r8,&  
         0.1865568E+03_r8,     0.1853745E+03_r8,     0.1841996E+03_r8,     0.1830321E+03_r8,     0.1818720E+03_r8,&  
         0.1807192E+03_r8,     0.1795737E+03_r8,     0.1784355E+03_r8,     0.1773045E+03_r8,     0.1761806E+03_r8,&  
         0.1750638E+03_r8,     0.1739541E+03_r8,     0.1728515E+03_r8,     0.1717558E+03_r8,     0.1706671E+03_r8,&  
         0.1695852E+03_r8,     0.1685102E+03_r8,     0.1674421E+03_r8,     0.1663807E+03_r8,     0.1653260E+03_r8,&  
         0.1642780E+03_r8,     0.1632367E+03_r8,     0.1622019E+03_r8,     0.1611737E+03_r8,     0.1601520E+03_r8,&  
         0.1591368E+03_r8,     0.1581281E+03_r8,     0.1571257E+03_r8,     0.1561297E+03_r8,     0.1551400E+03_r8,&  
         0.1541566E+03_r8,     0.1531794E+03_r8,     0.1522084E+03_r8,     0.1512435E+03_r8,     0.1502848E+03_r8,&  
         0.1493321E+03_r8,     0.1483855E+03_r8,     0.1474449E+03_r8,     0.1465103E+03_r8,     0.1455815E+03_r8,&  
         0.1446587E+03_r8,     0.1437417E+03_r8,     0.1428306E+03_r8,     0.1419251E+03_r8,     0.1410255E+03_r8,&  
         0.1401315E+03_r8,     0.1392432E+03_r8,     0.1383606E+03_r8,     0.1374835E+03_r8,     0.1366120E+03_r8,&  
         0.1357460E+03_r8,     0.1348855E+03_r8,     0.1340305E+03_r8,     0.1331809E+03_r8,     0.1323367E+03_r8,&  
         0.1314978E+03_r8,     0.1306642E+03_r8,     0.1298359E+03_r8,     0.1290129E+03_r8,     0.1281951E+03_r8,&  
         0.1273825E+03_r8,     0.1265750E+03_r8,     0.1257726E+03_r8,     0.1249754E+03_r8,     0.1241832E+03_r8,&  
         0.1233960E+03_r8,     0.1226138E+03_r8,     0.1218365E+03_r8,     0.1210642E+03_r8,     0.1202968E+03_r8,&  
         0.1195342E+03_r8,     0.1187765E+03_r8,     0.1180236E+03_r8,     0.1172754E+03_r8,     0.1165320E+03_r8,&  
         0.1157933E+03_r8,     0.1150593E+03_r8,     0.1143299E+03_r8,     0.1136052E+03_r8,     0.1128851E+03_r8,&  
         0.1121695E+03_r8,     0.1114585E+03_r8,     0.1107519E+03_r8,     0.1100499E+03_r8,     0.1093523E+03_r8,&  
         0.1086591E+03_r8,     0.1079703E+03_r8,     0.1072859E+03_r8,     0.1066058E+03_r8,     0.1059300E+03_r8,&  
         0.1052585E+03_r8,     0.1045913E+03_r8,     0.1039283E+03_r8,     0.1032695E+03_r8,     0.1026149E+03_r8,&  
         0.1019644E+03_r8,     0.1013181E+03_r8,     0.1006758E+03_r8,     0.1000376E+03_r8,     0.9940348E+02_r8,&  
         0.9877336E+02_r8,     0.9814725E+02_r8,     0.9752509E+02_r8,     0.9690688E+02_r8,     0.9629259E+02_r8,&  
         0.9568220E+02_r8,     0.9507566E+02_r8,     0.9447298E+02_r8,     0.9387412E+02_r8,     0.9327905E+02_r8/)  

    psaditmk(1:150,   21)= (/ &
         0.2418880E+03_r8,     0.2404285E+03_r8,     0.2389712E+03_r8,     0.2375164E+03_r8,     0.2360648E+03_r8,&  
         0.2346168E+03_r8,     0.2331729E+03_r8,     0.2317335E+03_r8,     0.2302991E+03_r8,     0.2288700E+03_r8,&  
         0.2274464E+03_r8,     0.2260289E+03_r8,     0.2246175E+03_r8,     0.2232126E+03_r8,     0.2218143E+03_r8,&  
         0.2204229E+03_r8,     0.2190385E+03_r8,     0.2176614E+03_r8,     0.2162915E+03_r8,     0.2149292E+03_r8,&  
         0.2135743E+03_r8,     0.2122271E+03_r8,     0.2108875E+03_r8,     0.2095557E+03_r8,     0.2082316E+03_r8,&  
         0.2069154E+03_r8,     0.2056070E+03_r8,     0.2043064E+03_r8,     0.2030138E+03_r8,     0.2017289E+03_r8,&  
         0.2004520E+03_r8,     0.1991828E+03_r8,     0.1979215E+03_r8,     0.1966680E+03_r8,     0.1954223E+03_r8,&  
         0.1941843E+03_r8,     0.1929541E+03_r8,     0.1917315E+03_r8,     0.1905166E+03_r8,     0.1893094E+03_r8,&  
         0.1881097E+03_r8,     0.1869176E+03_r8,     0.1857330E+03_r8,     0.1845558E+03_r8,     0.1833861E+03_r8,&  
         0.1822238E+03_r8,     0.1810688E+03_r8,     0.1799211E+03_r8,     0.1787806E+03_r8,     0.1776474E+03_r8,&  
         0.1765214E+03_r8,     0.1754025E+03_r8,     0.1742906E+03_r8,     0.1731858E+03_r8,     0.1720880E+03_r8,&  
         0.1709972E+03_r8,     0.1699133E+03_r8,     0.1688362E+03_r8,     0.1677660E+03_r8,     0.1667025E+03_r8,&  
         0.1656458E+03_r8,     0.1645958E+03_r8,     0.1635524E+03_r8,     0.1625157E+03_r8,     0.1614855E+03_r8,&  
         0.1604618E+03_r8,     0.1594447E+03_r8,     0.1584339E+03_r8,     0.1574296E+03_r8,     0.1564317E+03_r8,&  
         0.1554401E+03_r8,     0.1544548E+03_r8,     0.1534757E+03_r8,     0.1525028E+03_r8,     0.1515361E+03_r8,&  
         0.1505755E+03_r8,     0.1496210E+03_r8,     0.1486726E+03_r8,     0.1477301E+03_r8,     0.1467937E+03_r8,&  
         0.1458631E+03_r8,     0.1449385E+03_r8,     0.1440198E+03_r8,     0.1431068E+03_r8,     0.1421997E+03_r8,&  
         0.1412983E+03_r8,     0.1404026E+03_r8,     0.1395126E+03_r8,     0.1386282E+03_r8,     0.1377495E+03_r8,&  
         0.1368763E+03_r8,     0.1360086E+03_r8,     0.1351464E+03_r8,     0.1342898E+03_r8,     0.1334385E+03_r8,&  
         0.1325926E+03_r8,     0.1317522E+03_r8,     0.1309170E+03_r8,     0.1300871E+03_r8,     0.1292625E+03_r8,&  
         0.1284431E+03_r8,     0.1276289E+03_r8,     0.1268198E+03_r8,     0.1260159E+03_r8,     0.1252171E+03_r8,&  
         0.1244234E+03_r8,     0.1236347E+03_r8,     0.1228509E+03_r8,     0.1220722E+03_r8,     0.1212984E+03_r8,&  
         0.1205295E+03_r8,     0.1197654E+03_r8,     0.1190062E+03_r8,     0.1182519E+03_r8,     0.1175023E+03_r8,&  
         0.1167574E+03_r8,     0.1160173E+03_r8,     0.1152819E+03_r8,     0.1145511E+03_r8,     0.1138250E+03_r8,&  
         0.1131034E+03_r8,     0.1123865E+03_r8,     0.1116741E+03_r8,     0.1109662E+03_r8,     0.1102628E+03_r8,&  
         0.1095638E+03_r8,     0.1088693E+03_r8,     0.1081792E+03_r8,     0.1074934E+03_r8,     0.1068120E+03_r8,&  
         0.1061349E+03_r8,     0.1054621E+03_r8,     0.1047936E+03_r8,     0.1041293E+03_r8,     0.1034693E+03_r8,&  
         0.1028134E+03_r8,     0.1021616E+03_r8,     0.1015140E+03_r8,     0.1008705E+03_r8,     0.1002311E+03_r8,&  
         0.9959576E+02_r8,     0.9896443E+02_r8,     0.9833710E+02_r8,     0.9771375E+02_r8,     0.9709434E+02_r8,&  
         0.9647886E+02_r8,     0.9586728E+02_r8,     0.9525957E+02_r8,     0.9465573E+02_r8,     0.9405571E+02_r8/)  

    psaditmk(1:150,   22)= (/ &
         0.2437707E+03_r8,     0.2423118E+03_r8,     0.2408539E+03_r8,     0.2393978E+03_r8,     0.2379439E+03_r8,&  
         0.2364928E+03_r8,     0.2350451E+03_r8,     0.2336012E+03_r8,     0.2321616E+03_r8,     0.2307267E+03_r8,&  
         0.2292968E+03_r8,     0.2278724E+03_r8,     0.2264537E+03_r8,     0.2250411E+03_r8,     0.2236348E+03_r8,&  
         0.2222350E+03_r8,     0.2208419E+03_r8,     0.2194558E+03_r8,     0.2180768E+03_r8,     0.2167050E+03_r8,&  
         0.2153407E+03_r8,     0.2139837E+03_r8,     0.2126344E+03_r8,     0.2112926E+03_r8,     0.2099586E+03_r8,&  
         0.2086324E+03_r8,     0.2073139E+03_r8,     0.2060032E+03_r8,     0.2047004E+03_r8,     0.2034053E+03_r8,&  
         0.2021182E+03_r8,     0.2008389E+03_r8,     0.1995674E+03_r8,     0.1983038E+03_r8,     0.1970479E+03_r8,&  
         0.1957999E+03_r8,     0.1945596E+03_r8,     0.1933270E+03_r8,     0.1921021E+03_r8,     0.1908849E+03_r8,&  
         0.1896754E+03_r8,     0.1884734E+03_r8,     0.1872790E+03_r8,     0.1860921E+03_r8,     0.1849127E+03_r8,&  
         0.1837408E+03_r8,     0.1825762E+03_r8,     0.1814190E+03_r8,     0.1802691E+03_r8,     0.1791264E+03_r8,&  
         0.1779910E+03_r8,     0.1768628E+03_r8,     0.1757417E+03_r8,     0.1746278E+03_r8,     0.1735208E+03_r8,&  
         0.1724209E+03_r8,     0.1713280E+03_r8,     0.1702419E+03_r8,     0.1691628E+03_r8,     0.1680905E+03_r8,&  
         0.1670250E+03_r8,     0.1659662E+03_r8,     0.1649142E+03_r8,     0.1638688E+03_r8,     0.1628300E+03_r8,&  
         0.1617979E+03_r8,     0.1607722E+03_r8,     0.1597531E+03_r8,     0.1587404E+03_r8,     0.1577342E+03_r8,&  
         0.1567343E+03_r8,     0.1557408E+03_r8,     0.1547536E+03_r8,     0.1537726E+03_r8,     0.1527978E+03_r8,&  
         0.1518292E+03_r8,     0.1508668E+03_r8,     0.1499104E+03_r8,     0.1489602E+03_r8,     0.1480159E+03_r8,&  
         0.1470776E+03_r8,     0.1461453E+03_r8,     0.1452189E+03_r8,     0.1442984E+03_r8,     0.1433837E+03_r8,&  
         0.1424748E+03_r8,     0.1415716E+03_r8,     0.1406742E+03_r8,     0.1397825E+03_r8,     0.1388964E+03_r8,&  
         0.1380159E+03_r8,     0.1371411E+03_r8,     0.1362717E+03_r8,     0.1354079E+03_r8,     0.1345496E+03_r8,&  
         0.1336966E+03_r8,     0.1328491E+03_r8,     0.1320070E+03_r8,     0.1311702E+03_r8,     0.1303387E+03_r8,&  
         0.1295125E+03_r8,     0.1286915E+03_r8,     0.1278758E+03_r8,     0.1270652E+03_r8,     0.1262597E+03_r8,&  
         0.1254593E+03_r8,     0.1246641E+03_r8,     0.1238738E+03_r8,     0.1230886E+03_r8,     0.1223083E+03_r8,&  
         0.1215330E+03_r8,     0.1207626E+03_r8,     0.1199971E+03_r8,     0.1192365E+03_r8,     0.1184806E+03_r8,&  
         0.1177296E+03_r8,     0.1169833E+03_r8,     0.1162417E+03_r8,     0.1155049E+03_r8,     0.1147727E+03_r8,&  
         0.1140452E+03_r8,     0.1133222E+03_r8,     0.1126039E+03_r8,     0.1118901E+03_r8,     0.1111808E+03_r8,&  
         0.1104761E+03_r8,     0.1097757E+03_r8,     0.1090799E+03_r8,     0.1083884E+03_r8,     0.1077013E+03_r8,&  
         0.1070186E+03_r8,     0.1063402E+03_r8,     0.1056662E+03_r8,     0.1049963E+03_r8,     0.1043308E+03_r8,&  
         0.1036694E+03_r8,     0.1030123E+03_r8,     0.1023593E+03_r8,     0.1017104E+03_r8,     0.1010657E+03_r8,&  
         0.1004250E+03_r8,     0.9978844E+02_r8,     0.9915588E+02_r8,     0.9852733E+02_r8,     0.9790276E+02_r8,&  
         0.9728217E+02_r8,     0.9666549E+02_r8,     0.9605273E+02_r8,     0.9544386E+02_r8,     0.9483884E+02_r8/)  

    psaditmk(1:150,   23)= (/ &
         0.2456473E+03_r8,     0.2441905E+03_r8,     0.2427338E+03_r8,     0.2412778E+03_r8,     0.2398230E+03_r8,&  
         0.2383702E+03_r8,     0.2369199E+03_r8,     0.2354726E+03_r8,     0.2340289E+03_r8,     0.2325892E+03_r8,&  
         0.2311539E+03_r8,     0.2297234E+03_r8,     0.2282982E+03_r8,     0.2268785E+03_r8,     0.2254648E+03_r8,&  
         0.2240571E+03_r8,     0.2226558E+03_r8,     0.2212612E+03_r8,     0.2198734E+03_r8,     0.2184926E+03_r8,&  
         0.2171190E+03_r8,     0.2157526E+03_r8,     0.2143937E+03_r8,     0.2130422E+03_r8,     0.2116984E+03_r8,&  
         0.2103622E+03_r8,     0.2090337E+03_r8,     0.2077130E+03_r8,     0.2064000E+03_r8,     0.2050949E+03_r8,&  
         0.2037976E+03_r8,     0.2025081E+03_r8,     0.2012265E+03_r8,     0.1999527E+03_r8,     0.1986867E+03_r8,&  
         0.1974285E+03_r8,     0.1961781E+03_r8,     0.1949355E+03_r8,     0.1937007E+03_r8,     0.1924734E+03_r8,&  
         0.1912540E+03_r8,     0.1900421E+03_r8,     0.1888378E+03_r8,     0.1876411E+03_r8,     0.1864520E+03_r8,&  
         0.1852703E+03_r8,     0.1840961E+03_r8,     0.1829293E+03_r8,     0.1817699E+03_r8,     0.1806177E+03_r8,&  
         0.1794729E+03_r8,     0.1783353E+03_r8,     0.1772049E+03_r8,     0.1760817E+03_r8,     0.1749656E+03_r8,&  
         0.1738565E+03_r8,     0.1727544E+03_r8,     0.1716594E+03_r8,     0.1705713E+03_r8,     0.1694900E+03_r8,&  
         0.1684156E+03_r8,     0.1673481E+03_r8,     0.1662873E+03_r8,     0.1652332E+03_r8,     0.1641858E+03_r8,&  
         0.1631450E+03_r8,     0.1621109E+03_r8,     0.1610832E+03_r8,     0.1600621E+03_r8,     0.1590475E+03_r8,&  
         0.1580393E+03_r8,     0.1570375E+03_r8,     0.1560421E+03_r8,     0.1550529E+03_r8,     0.1540700E+03_r8,&  
         0.1530934E+03_r8,     0.1521229E+03_r8,     0.1511586E+03_r8,     0.1502004E+03_r8,     0.1492483E+03_r8,&  
         0.1483022E+03_r8,     0.1473622E+03_r8,     0.1464280E+03_r8,     0.1454998E+03_r8,     0.1445775E+03_r8,&  
         0.1436610E+03_r8,     0.1427504E+03_r8,     0.1418455E+03_r8,     0.1409463E+03_r8,     0.1400529E+03_r8,&  
         0.1391651E+03_r8,     0.1382829E+03_r8,     0.1374063E+03_r8,     0.1365353E+03_r8,     0.1356698E+03_r8,&  
         0.1348098E+03_r8,     0.1339553E+03_r8,     0.1331061E+03_r8,     0.1322623E+03_r8,     0.1314240E+03_r8,&  
         0.1305909E+03_r8,     0.1297630E+03_r8,     0.1289405E+03_r8,     0.1281232E+03_r8,     0.1273110E+03_r8,&  
         0.1265039E+03_r8,     0.1257020E+03_r8,     0.1249052E+03_r8,     0.1241134E+03_r8,     0.1233267E+03_r8,&  
         0.1225449E+03_r8,     0.1217681E+03_r8,     0.1209962E+03_r8,     0.1202292E+03_r8,     0.1194671E+03_r8,&  
         0.1187098E+03_r8,     0.1179573E+03_r8,     0.1172096E+03_r8,     0.1164666E+03_r8,     0.1157283E+03_r8,&  
         0.1149947E+03_r8,     0.1142658E+03_r8,     0.1135414E+03_r8,     0.1128217E+03_r8,     0.1121065E+03_r8,&  
         0.1113959E+03_r8,     0.1106898E+03_r8,     0.1099881E+03_r8,     0.1092909E+03_r8,     0.1085981E+03_r8,&  
         0.1079097E+03_r8,     0.1072256E+03_r8,     0.1065459E+03_r8,     0.1058706E+03_r8,     0.1051994E+03_r8,&  
         0.1045326E+03_r8,     0.1038700E+03_r8,     0.1032115E+03_r8,     0.1025573E+03_r8,     0.1019072E+03_r8,&  
         0.1012612E+03_r8,     0.1006193E+03_r8,     0.9998147E+02_r8,     0.9934769E+02_r8,     0.9871792E+02_r8,&  
         0.9809216E+02_r8,     0.9747034E+02_r8,     0.9685248E+02_r8,     0.9623854E+02_r8,     0.9562848E+02_r8/)  

    psaditmk(1:150,   24)= (/ &
         0.2475152E+03_r8,     0.2460623E+03_r8,     0.2446083E+03_r8,     0.2431540E+03_r8,     0.2417000E+03_r8,&  
         0.2402469E+03_r8,     0.2387954E+03_r8,     0.2373461E+03_r8,     0.2358994E+03_r8,     0.2344560E+03_r8,&  
         0.2330163E+03_r8,     0.2315808E+03_r8,     0.2301499E+03_r8,     0.2287240E+03_r8,     0.2273034E+03_r8,&  
         0.2258885E+03_r8,     0.2244796E+03_r8,     0.2230770E+03_r8,     0.2216809E+03_r8,     0.2202914E+03_r8,&  
         0.2189089E+03_r8,     0.2175334E+03_r8,     0.2161652E+03_r8,     0.2148042E+03_r8,     0.2134507E+03_r8,&  
         0.2121048E+03_r8,     0.2107664E+03_r8,     0.2094357E+03_r8,     0.2081128E+03_r8,     0.2067976E+03_r8,&  
         0.2054902E+03_r8,     0.2041906E+03_r8,     0.2028988E+03_r8,     0.2016149E+03_r8,     0.2003387E+03_r8,&  
         0.1990705E+03_r8,     0.1978099E+03_r8,     0.1965572E+03_r8,     0.1953123E+03_r8,     0.1940750E+03_r8,&  
         0.1928456E+03_r8,     0.1916237E+03_r8,     0.1904096E+03_r8,     0.1892030E+03_r8,     0.1880040E+03_r8,&  
         0.1868126E+03_r8,     0.1856287E+03_r8,     0.1844522E+03_r8,     0.1832832E+03_r8,     0.1821215E+03_r8,&  
         0.1809672E+03_r8,     0.1798201E+03_r8,     0.1786803E+03_r8,     0.1775477E+03_r8,     0.1764224E+03_r8,&  
         0.1753041E+03_r8,     0.1741929E+03_r8,     0.1730887E+03_r8,     0.1719915E+03_r8,     0.1709013E+03_r8,&  
         0.1698180E+03_r8,     0.1687415E+03_r8,     0.1676719E+03_r8,     0.1666090E+03_r8,     0.1655529E+03_r8,&  
         0.1645034E+03_r8,     0.1634607E+03_r8,     0.1624245E+03_r8,     0.1613949E+03_r8,     0.1603718E+03_r8,&  
         0.1593552E+03_r8,     0.1583451E+03_r8,     0.1573414E+03_r8,     0.1563440E+03_r8,     0.1553529E+03_r8,&  
         0.1543681E+03_r8,     0.1533896E+03_r8,     0.1524173E+03_r8,     0.1514511E+03_r8,     0.1504911E+03_r8,&  
         0.1495371E+03_r8,     0.1485892E+03_r8,     0.1476473E+03_r8,     0.1467113E+03_r8,     0.1457814E+03_r8,&  
         0.1448572E+03_r8,     0.1439390E+03_r8,     0.1430266E+03_r8,     0.1421199E+03_r8,     0.1412190E+03_r8,&  
         0.1403239E+03_r8,     0.1394344E+03_r8,     0.1385505E+03_r8,     0.1376722E+03_r8,     0.1367995E+03_r8,&  
         0.1359323E+03_r8,     0.1350707E+03_r8,     0.1342145E+03_r8,     0.1333636E+03_r8,     0.1325183E+03_r8,&  
         0.1316783E+03_r8,     0.1308435E+03_r8,     0.1300141E+03_r8,     0.1291900E+03_r8,     0.1283710E+03_r8,&  
         0.1275573E+03_r8,     0.1267487E+03_r8,     0.1259453E+03_r8,     0.1251469E+03_r8,     0.1243536E+03_r8,&  
         0.1235653E+03_r8,     0.1227820E+03_r8,     0.1220037E+03_r8,     0.1212303E+03_r8,     0.1204619E+03_r8,&  
         0.1196983E+03_r8,     0.1189395E+03_r8,     0.1181855E+03_r8,     0.1174364E+03_r8,     0.1166919E+03_r8,&  
         0.1159522E+03_r8,     0.1152172E+03_r8,     0.1144869E+03_r8,     0.1137611E+03_r8,     0.1130400E+03_r8,&  
         0.1123234E+03_r8,     0.1116114E+03_r8,     0.1109039E+03_r8,     0.1102009E+03_r8,     0.1095023E+03_r8,&  
         0.1088082E+03_r8,     0.1081185E+03_r8,     0.1074331E+03_r8,     0.1067521E+03_r8,     0.1060754E+03_r8,&  
         0.1054030E+03_r8,     0.1047349E+03_r8,     0.1040709E+03_r8,     0.1034112E+03_r8,     0.1027557E+03_r8,&  
         0.1021043E+03_r8,     0.1014571E+03_r8,     0.1008140E+03_r8,     0.1001749E+03_r8,     0.9953991E+02_r8,&  
         0.9890894E+02_r8,     0.9828194E+02_r8,     0.9765894E+02_r8,     0.9703989E+02_r8,     0.9642474E+02_r8/)  

    psaditmk(1:150,   25)= (/ &
         0.2493711E+03_r8,     0.2479239E+03_r8,     0.2464746E+03_r8,     0.2450238E+03_r8,     0.2435721E+03_r8,&  
         0.2421204E+03_r8,     0.2406692E+03_r8,     0.2392191E+03_r8,     0.2377710E+03_r8,     0.2363251E+03_r8,&  
         0.2348821E+03_r8,     0.2334427E+03_r8,     0.2320070E+03_r8,     0.2305758E+03_r8,     0.2291494E+03_r8,&  
         0.2277280E+03_r8,     0.2263122E+03_r8,     0.2249022E+03_r8,     0.2234982E+03_r8,     0.2221006E+03_r8,&  
         0.2207097E+03_r8,     0.2193254E+03_r8,     0.2179482E+03_r8,     0.2165780E+03_r8,     0.2152152E+03_r8,&  
         0.2138597E+03_r8,     0.2125116E+03_r8,     0.2111711E+03_r8,     0.2098383E+03_r8,     0.2085131E+03_r8,&  
         0.2071957E+03_r8,     0.2058860E+03_r8,     0.2045841E+03_r8,     0.2032901E+03_r8,     0.2020038E+03_r8,&  
         0.2007254E+03_r8,     0.1994548E+03_r8,     0.1981919E+03_r8,     0.1969369E+03_r8,     0.1956896E+03_r8,&  
         0.1944501E+03_r8,     0.1932183E+03_r8,     0.1919942E+03_r8,     0.1907777E+03_r8,     0.1895688E+03_r8,&  
         0.1883676E+03_r8,     0.1871739E+03_r8,     0.1859876E+03_r8,     0.1848089E+03_r8,     0.1836376E+03_r8,&  
         0.1824737E+03_r8,     0.1813172E+03_r8,     0.1801679E+03_r8,     0.1790259E+03_r8,     0.1778912E+03_r8,&  
         0.1767636E+03_r8,     0.1756432E+03_r8,     0.1745298E+03_r8,     0.1734235E+03_r8,     0.1723242E+03_r8,&  
         0.1712319E+03_r8,     0.1701465E+03_r8,     0.1690679E+03_r8,     0.1679962E+03_r8,     0.1669313E+03_r8,&  
         0.1658731E+03_r8,     0.1648217E+03_r8,     0.1637769E+03_r8,     0.1627387E+03_r8,     0.1617071E+03_r8,&  
         0.1606821E+03_r8,     0.1596635E+03_r8,     0.1586514E+03_r8,     0.1576457E+03_r8,     0.1566464E+03_r8,&  
         0.1556534E+03_r8,     0.1546668E+03_r8,     0.1536863E+03_r8,     0.1527121E+03_r8,     0.1517441E+03_r8,&  
         0.1507822E+03_r8,     0.1498264E+03_r8,     0.1488766E+03_r8,     0.1479329E+03_r8,     0.1469952E+03_r8,&  
         0.1460634E+03_r8,     0.1451375E+03_r8,     0.1442174E+03_r8,     0.1433032E+03_r8,     0.1423948E+03_r8,&  
         0.1414922E+03_r8,     0.1405953E+03_r8,     0.1397041E+03_r8,     0.1388185E+03_r8,     0.1379385E+03_r8,&  
         0.1370641E+03_r8,     0.1361953E+03_r8,     0.1353320E+03_r8,     0.1344741E+03_r8,     0.1336217E+03_r8,&  
         0.1327746E+03_r8,     0.1319330E+03_r8,     0.1310966E+03_r8,     0.1302657E+03_r8,     0.1294399E+03_r8,&  
         0.1286194E+03_r8,     0.1278041E+03_r8,     0.1269939E+03_r8,     0.1261889E+03_r8,     0.1253890E+03_r8,&  
         0.1245941E+03_r8,     0.1238043E+03_r8,     0.1230196E+03_r8,     0.1222397E+03_r8,     0.1214649E+03_r8,&  
         0.1206949E+03_r8,     0.1199298E+03_r8,     0.1191696E+03_r8,     0.1184142E+03_r8,     0.1176636E+03_r8,&  
         0.1169177E+03_r8,     0.1161765E+03_r8,     0.1154401E+03_r8,     0.1147083E+03_r8,     0.1139812E+03_r8,&  
         0.1132587E+03_r8,     0.1125407E+03_r8,     0.1118273E+03_r8,     0.1111185E+03_r8,     0.1104141E+03_r8,&  
         0.1097142E+03_r8,     0.1090187E+03_r8,     0.1083276E+03_r8,     0.1076409E+03_r8,     0.1069586E+03_r8,&  
         0.1062806E+03_r8,     0.1056069E+03_r8,     0.1049374E+03_r8,     0.1042723E+03_r8,     0.1036113E+03_r8,&  
         0.1029545E+03_r8,     0.1023019E+03_r8,     0.1016534E+03_r8,     0.1010090E+03_r8,     0.1003687E+03_r8,&  
         0.9973247E+02_r8,     0.9910026E+02_r8,     0.9847207E+02_r8,     0.9784785E+02_r8,     0.9722760E+02_r8/)  

    psaditmk(1:150,   26)= (/ &
         0.2512121E+03_r8,     0.2497727E+03_r8,     0.2483298E+03_r8,     0.2468843E+03_r8,     0.2454368E+03_r8,&  
         0.2439881E+03_r8,     0.2425388E+03_r8,     0.2410897E+03_r8,     0.2396415E+03_r8,     0.2381945E+03_r8,&  
         0.2367497E+03_r8,     0.2353074E+03_r8,     0.2338683E+03_r8,     0.2324328E+03_r8,     0.2310013E+03_r8,&  
         0.2295744E+03_r8,     0.2281524E+03_r8,     0.2267358E+03_r8,     0.2253247E+03_r8,     0.2239195E+03_r8,&  
         0.2225206E+03_r8,     0.2211281E+03_r8,     0.2197423E+03_r8,     0.2183633E+03_r8,     0.2169913E+03_r8,&  
         0.2156266E+03_r8,     0.2142691E+03_r8,     0.2129190E+03_r8,     0.2115764E+03_r8,     0.2102414E+03_r8,&  
         0.2089141E+03_r8,     0.2075944E+03_r8,     0.2062825E+03_r8,     0.2049784E+03_r8,     0.2036820E+03_r8,&  
         0.2023935E+03_r8,     0.2011127E+03_r8,     0.1998398E+03_r8,     0.1985746E+03_r8,     0.1973173E+03_r8,&  
         0.1960677E+03_r8,     0.1948258E+03_r8,     0.1935917E+03_r8,     0.1923652E+03_r8,     0.1911465E+03_r8,&  
         0.1899353E+03_r8,     0.1887318E+03_r8,     0.1875358E+03_r8,     0.1863473E+03_r8,     0.1851663E+03_r8,&  
         0.1839928E+03_r8,     0.1828267E+03_r8,     0.1816679E+03_r8,     0.1805164E+03_r8,     0.1793722E+03_r8,&  
         0.1782353E+03_r8,     0.1771055E+03_r8,     0.1759829E+03_r8,     0.1748674E+03_r8,     0.1737590E+03_r8,&  
         0.1726576E+03_r8,     0.1715631E+03_r8,     0.1704756E+03_r8,     0.1693950E+03_r8,     0.1683212E+03_r8,&  
         0.1672542E+03_r8,     0.1661940E+03_r8,     0.1651405E+03_r8,     0.1640937E+03_r8,     0.1630535E+03_r8,&  
         0.1620199E+03_r8,     0.1609929E+03_r8,     0.1599724E+03_r8,     0.1589583E+03_r8,     0.1579507E+03_r8,&  
         0.1569494E+03_r8,     0.1559545E+03_r8,     0.1549659E+03_r8,     0.1539836E+03_r8,     0.1530075E+03_r8,&  
         0.1520376E+03_r8,     0.1510739E+03_r8,     0.1501162E+03_r8,     0.1491646E+03_r8,     0.1482191E+03_r8,&  
         0.1472795E+03_r8,     0.1463459E+03_r8,     0.1454182E+03_r8,     0.1444964E+03_r8,     0.1435805E+03_r8,&  
         0.1426703E+03_r8,     0.1417659E+03_r8,     0.1408673E+03_r8,     0.1399743E+03_r8,     0.1390871E+03_r8,&  
         0.1382054E+03_r8,     0.1373293E+03_r8,     0.1364588E+03_r8,     0.1355937E+03_r8,     0.1347342E+03_r8,&  
         0.1338802E+03_r8,     0.1330315E+03_r8,     0.1321882E+03_r8,     0.1313503E+03_r8,     0.1305176E+03_r8,&  
         0.1296903E+03_r8,     0.1288682E+03_r8,     0.1280513E+03_r8,     0.1272396E+03_r8,     0.1264330E+03_r8,&  
         0.1256315E+03_r8,     0.1248352E+03_r8,     0.1240438E+03_r8,     0.1232575E+03_r8,     0.1224762E+03_r8,&  
         0.1216998E+03_r8,     0.1209284E+03_r8,     0.1201618E+03_r8,     0.1194001E+03_r8,     0.1186432E+03_r8,&  
         0.1178912E+03_r8,     0.1171439E+03_r8,     0.1164013E+03_r8,     0.1156634E+03_r8,     0.1149302E+03_r8,&  
         0.1142017E+03_r8,     0.1134778E+03_r8,     0.1127584E+03_r8,     0.1120437E+03_r8,     0.1113334E+03_r8,&  
         0.1106277E+03_r8,     0.1099264E+03_r8,     0.1092296E+03_r8,     0.1085372E+03_r8,     0.1078492E+03_r8,&  
         0.1071655E+03_r8,     0.1064862E+03_r8,     0.1058112E+03_r8,     0.1051404E+03_r8,     0.1044740E+03_r8,&  
         0.1038117E+03_r8,     0.1031536E+03_r8,     0.1024998E+03_r8,     0.1018500E+03_r8,     0.1012044E+03_r8,&  
         0.1005629E+03_r8,     0.9992539E+02_r8,     0.9929197E+02_r8,     0.9866256E+02_r8,     0.9803714E+02_r8/)  

    psaditmk(1:150,   27)= (/ &
         0.2530353E+03_r8,     0.2516054E+03_r8,     0.2501709E+03_r8,     0.2487326E+03_r8,     0.2472911E+03_r8,&  
         0.2458473E+03_r8,     0.2444017E+03_r8,     0.2429552E+03_r8,     0.2415085E+03_r8,     0.2400621E+03_r8,&  
         0.2386168E+03_r8,     0.2371731E+03_r8,     0.2357317E+03_r8,     0.2342931E+03_r8,     0.2328578E+03_r8,&  
         0.2314263E+03_r8,     0.2299991E+03_r8,     0.2285766E+03_r8,     0.2271592E+03_r8,     0.2257472E+03_r8,&  
         0.2243409E+03_r8,     0.2229407E+03_r8,     0.2215468E+03_r8,     0.2201594E+03_r8,     0.2187787E+03_r8,&  
         0.2174050E+03_r8,     0.2160384E+03_r8,     0.2146789E+03_r8,     0.2133269E+03_r8,     0.2119822E+03_r8,&  
         0.2106451E+03_r8,     0.2093156E+03_r8,     0.2079938E+03_r8,     0.2066796E+03_r8,     0.2053732E+03_r8,&  
         0.2040746E+03_r8,     0.2027838E+03_r8,     0.2015007E+03_r8,     0.2002255E+03_r8,     0.1989580E+03_r8,&  
         0.1976983E+03_r8,     0.1964464E+03_r8,     0.1952022E+03_r8,     0.1939658E+03_r8,     0.1927370E+03_r8,&  
         0.1915160E+03_r8,     0.1903025E+03_r8,     0.1890966E+03_r8,     0.1878984E+03_r8,     0.1867077E+03_r8,&  
         0.1855244E+03_r8,     0.1843486E+03_r8,     0.1831802E+03_r8,     0.1820192E+03_r8,     0.1808656E+03_r8,&  
         0.1797192E+03_r8,     0.1785800E+03_r8,     0.1774481E+03_r8,     0.1763233E+03_r8,     0.1752056E+03_r8,&  
         0.1740951E+03_r8,     0.1729915E+03_r8,     0.1718950E+03_r8,     0.1708054E+03_r8,     0.1697227E+03_r8,&  
         0.1686468E+03_r8,     0.1675778E+03_r8,     0.1665155E+03_r8,     0.1654600E+03_r8,     0.1644111E+03_r8,&  
         0.1633689E+03_r8,     0.1623333E+03_r8,     0.1613044E+03_r8,     0.1602818E+03_r8,     0.1592658E+03_r8,&  
         0.1582562E+03_r8,     0.1572530E+03_r8,     0.1562562E+03_r8,     0.1552657E+03_r8,     0.1542815E+03_r8,&  
         0.1533035E+03_r8,     0.1523317E+03_r8,     0.1513661E+03_r8,     0.1504066E+03_r8,     0.1494532E+03_r8,&  
         0.1485058E+03_r8,     0.1475644E+03_r8,     0.1466290E+03_r8,     0.1456995E+03_r8,     0.1447759E+03_r8,&  
         0.1438582E+03_r8,     0.1429463E+03_r8,     0.1420401E+03_r8,     0.1411398E+03_r8,     0.1402451E+03_r8,&  
         0.1393561E+03_r8,     0.1384727E+03_r8,     0.1375949E+03_r8,     0.1367227E+03_r8,     0.1358560E+03_r8,&  
         0.1349949E+03_r8,     0.1341391E+03_r8,     0.1332888E+03_r8,     0.1324439E+03_r8,     0.1316043E+03_r8,&  
         0.1307701E+03_r8,     0.1299412E+03_r8,     0.1291175E+03_r8,     0.1282990E+03_r8,     0.1274857E+03_r8,&  
         0.1266776E+03_r8,     0.1258746E+03_r8,     0.1250767E+03_r8,     0.1242838E+03_r8,     0.1234960E+03_r8,&  
         0.1227131E+03_r8,     0.1219353E+03_r8,     0.1211623E+03_r8,     0.1203943E+03_r8,     0.1196311E+03_r8,&  
         0.1188727E+03_r8,     0.1181192E+03_r8,     0.1173705E+03_r8,     0.1166264E+03_r8,     0.1158872E+03_r8,&  
         0.1151525E+03_r8,     0.1144226E+03_r8,     0.1136973E+03_r8,     0.1129766E+03_r8,     0.1122604E+03_r8,&  
         0.1115488E+03_r8,     0.1108417E+03_r8,     0.1101390E+03_r8,     0.1094409E+03_r8,     0.1087471E+03_r8,&  
         0.1080578E+03_r8,     0.1073728E+03_r8,     0.1066922E+03_r8,     0.1060159E+03_r8,     0.1053438E+03_r8,&  
         0.1046761E+03_r8,     0.1040125E+03_r8,     0.1033532E+03_r8,     0.1026980E+03_r8,     0.1020470E+03_r8,&  
         0.1014002E+03_r8,     0.1007574E+03_r8,     0.1001187E+03_r8,     0.9948404E+02_r8,     0.9885341E+02_r8/)  

    psaditmk(1:150,   28)= (/ &
         0.2548374E+03_r8,     0.2534189E+03_r8,     0.2519948E+03_r8,     0.2505656E+03_r8,     0.2491321E+03_r8,&  
         0.2476950E+03_r8,     0.2462550E+03_r8,     0.2448129E+03_r8,     0.2433695E+03_r8,     0.2419252E+03_r8,&  
         0.2404810E+03_r8,     0.2390375E+03_r8,     0.2375952E+03_r8,     0.2361549E+03_r8,     0.2347170E+03_r8,&  
         0.2332821E+03_r8,     0.2318507E+03_r8,     0.2304234E+03_r8,     0.2290005E+03_r8,     0.2275824E+03_r8,&  
         0.2261696E+03_r8,     0.2247623E+03_r8,     0.2233609E+03_r8,     0.2219656E+03_r8,     0.2205767E+03_r8,&  
         0.2191945E+03_r8,     0.2178190E+03_r8,     0.2164506E+03_r8,     0.2150893E+03_r8,     0.2137352E+03_r8,&  
         0.2123885E+03_r8,     0.2110494E+03_r8,     0.2097177E+03_r8,     0.2083937E+03_r8,     0.2070774E+03_r8,&  
         0.2057687E+03_r8,     0.2044678E+03_r8,     0.2031747E+03_r8,     0.2018894E+03_r8,     0.2006118E+03_r8,&  
         0.1993421E+03_r8,     0.1980800E+03_r8,     0.1968258E+03_r8,     0.1955793E+03_r8,     0.1943405E+03_r8,&  
         0.1931095E+03_r8,     0.1918861E+03_r8,     0.1906704E+03_r8,     0.1894622E+03_r8,     0.1882617E+03_r8,&  
         0.1870687E+03_r8,     0.1858832E+03_r8,     0.1847051E+03_r8,     0.1835345E+03_r8,     0.1823712E+03_r8,&  
         0.1812153E+03_r8,     0.1800668E+03_r8,     0.1789254E+03_r8,     0.1777913E+03_r8,     0.1766644E+03_r8,&  
         0.1755446E+03_r8,     0.1744318E+03_r8,     0.1733262E+03_r8,     0.1722275E+03_r8,     0.1711358E+03_r8,&  
         0.1700509E+03_r8,     0.1689730E+03_r8,     0.1679019E+03_r8,     0.1668376E+03_r8,     0.1657800E+03_r8,&  
         0.1647292E+03_r8,     0.1636850E+03_r8,     0.1626474E+03_r8,     0.1616164E+03_r8,     0.1605919E+03_r8,&  
         0.1595739E+03_r8,     0.1585624E+03_r8,     0.1575572E+03_r8,     0.1565585E+03_r8,     0.1555661E+03_r8,&  
         0.1545799E+03_r8,     0.1536001E+03_r8,     0.1526264E+03_r8,     0.1516589E+03_r8,     0.1506976E+03_r8,&  
         0.1497423E+03_r8,     0.1487931E+03_r8,     0.1478499E+03_r8,     0.1469126E+03_r8,     0.1459814E+03_r8,&  
         0.1450560E+03_r8,     0.1441365E+03_r8,     0.1432228E+03_r8,     0.1423149E+03_r8,     0.1414128E+03_r8,&  
         0.1405164E+03_r8,     0.1396257E+03_r8,     0.1387406E+03_r8,     0.1378611E+03_r8,     0.1369872E+03_r8,&  
         0.1361189E+03_r8,     0.1352560E+03_r8,     0.1343986E+03_r8,     0.1335467E+03_r8,     0.1327001E+03_r8,&  
         0.1318589E+03_r8,     0.1310231E+03_r8,     0.1301925E+03_r8,     0.1293672E+03_r8,     0.1285472E+03_r8,&  
         0.1277323E+03_r8,     0.1269226E+03_r8,     0.1261181E+03_r8,     0.1253186E+03_r8,     0.1245242E+03_r8,&  
         0.1237349E+03_r8,     0.1229505E+03_r8,     0.1221711E+03_r8,     0.1213967E+03_r8,     0.1206272E+03_r8,&  
         0.1198625E+03_r8,     0.1191027E+03_r8,     0.1183477E+03_r8,     0.1175975E+03_r8,     0.1168521E+03_r8,&  
         0.1161113E+03_r8,     0.1153753E+03_r8,     0.1146439E+03_r8,     0.1139172E+03_r8,     0.1131951E+03_r8,&  
         0.1124776E+03_r8,     0.1117646E+03_r8,     0.1110561E+03_r8,     0.1103521E+03_r8,     0.1096526E+03_r8,&  
         0.1089575E+03_r8,     0.1082668E+03_r8,     0.1075805E+03_r8,     0.1068986E+03_r8,     0.1062209E+03_r8,&  
         0.1055476E+03_r8,     0.1048786E+03_r8,     0.1042137E+03_r8,     0.1035531E+03_r8,     0.1028967E+03_r8,&  
         0.1022444E+03_r8,     0.1015963E+03_r8,     0.1009523E+03_r8,     0.1003124E+03_r8,     0.9967648E+02_r8/)  

    psaditmk(1:150,   29)= (/ &
         0.2566152E+03_r8,     0.2552103E+03_r8,     0.2537984E+03_r8,     0.2523802E+03_r8,     0.2509566E+03_r8,&  
         0.2495282E+03_r8,     0.2480957E+03_r8,     0.2466599E+03_r8,     0.2452216E+03_r8,     0.2437813E+03_r8,&  
         0.2423400E+03_r8,     0.2408982E+03_r8,     0.2394566E+03_r8,     0.2380160E+03_r8,     0.2365768E+03_r8,&  
         0.2351398E+03_r8,     0.2337055E+03_r8,     0.2322745E+03_r8,     0.2308471E+03_r8,     0.2294239E+03_r8,&  
         0.2280054E+03_r8,     0.2265918E+03_r8,     0.2251836E+03_r8,     0.2237811E+03_r8,     0.2223846E+03_r8,&  
         0.2209943E+03_r8,     0.2196105E+03_r8,     0.2182334E+03_r8,     0.2168632E+03_r8,     0.2155000E+03_r8,&  
         0.2141441E+03_r8,     0.2127954E+03_r8,     0.2114541E+03_r8,     0.2101204E+03_r8,     0.2087943E+03_r8,&  
         0.2074757E+03_r8,     0.2061648E+03_r8,     0.2048617E+03_r8,     0.2035663E+03_r8,     0.2022787E+03_r8,&  
         0.2009989E+03_r8,     0.1997268E+03_r8,     0.1984625E+03_r8,     0.1972059E+03_r8,     0.1959571E+03_r8,&  
         0.1947161E+03_r8,     0.1934827E+03_r8,     0.1922570E+03_r8,     0.1910389E+03_r8,     0.1898285E+03_r8,&  
         0.1886257E+03_r8,     0.1874304E+03_r8,     0.1862426E+03_r8,     0.1850623E+03_r8,     0.1838894E+03_r8,&  
         0.1827240E+03_r8,     0.1815658E+03_r8,     0.1804151E+03_r8,     0.1792715E+03_r8,     0.1781352E+03_r8,&  
         0.1770061E+03_r8,     0.1758841E+03_r8,     0.1747692E+03_r8,     0.1736614E+03_r8,     0.1725606E+03_r8,&  
         0.1714668E+03_r8,     0.1703799E+03_r8,     0.1692999E+03_r8,     0.1682267E+03_r8,     0.1671603E+03_r8,&  
         0.1661007E+03_r8,     0.1650478E+03_r8,     0.1640016E+03_r8,     0.1629620E+03_r8,     0.1619290E+03_r8,&  
         0.1609025E+03_r8,     0.1598826E+03_r8,     0.1588691E+03_r8,     0.1578620E+03_r8,     0.1568613E+03_r8,&  
         0.1558670E+03_r8,     0.1548790E+03_r8,     0.1538972E+03_r8,     0.1529216E+03_r8,     0.1519523E+03_r8,&  
         0.1509891E+03_r8,     0.1500319E+03_r8,     0.1490809E+03_r8,     0.1481359E+03_r8,     0.1471968E+03_r8,&  
         0.1462638E+03_r8,     0.1453366E+03_r8,     0.1444153E+03_r8,     0.1434999E+03_r8,     0.1425903E+03_r8,&  
         0.1416863E+03_r8,     0.1407882E+03_r8,     0.1398958E+03_r8,     0.1390089E+03_r8,     0.1381278E+03_r8,&  
         0.1372522E+03_r8,     0.1363822E+03_r8,     0.1355176E+03_r8,     0.1346586E+03_r8,     0.1338050E+03_r8,&  
         0.1329568E+03_r8,     0.1321140E+03_r8,     0.1312765E+03_r8,     0.1304444E+03_r8,     0.1296175E+03_r8,&  
         0.1287959E+03_r8,     0.1279794E+03_r8,     0.1271682E+03_r8,     0.1263620E+03_r8,     0.1255610E+03_r8,&  
         0.1247651E+03_r8,     0.1239742E+03_r8,     0.1231883E+03_r8,     0.1224075E+03_r8,     0.1216315E+03_r8,&  
         0.1208605E+03_r8,     0.1200944E+03_r8,     0.1193331E+03_r8,     0.1185766E+03_r8,     0.1178250E+03_r8,&  
         0.1170781E+03_r8,     0.1163359E+03_r8,     0.1155985E+03_r8,     0.1148657E+03_r8,     0.1141376E+03_r8,&  
         0.1134141E+03_r8,     0.1126951E+03_r8,     0.1119808E+03_r8,     0.1112709E+03_r8,     0.1105656E+03_r8,&  
         0.1098647E+03_r8,     0.1091683E+03_r8,     0.1084763E+03_r8,     0.1077886E+03_r8,     0.1071054E+03_r8,&  
         0.1064264E+03_r8,     0.1057518E+03_r8,     0.1050814E+03_r8,     0.1044153E+03_r8,     0.1037534E+03_r8,&  
         0.1030958E+03_r8,     0.1024422E+03_r8,     0.1017928E+03_r8,     0.1011476E+03_r8,     0.1005064E+03_r8/)  

    psaditmk(1:150,   30)= (/ &
         0.2583661E+03_r8,     0.2569763E+03_r8,     0.2555784E+03_r8,     0.2541734E+03_r8,     0.2527615E+03_r8,&  
         0.2513438E+03_r8,     0.2499207E+03_r8,     0.2484932E+03_r8,     0.2470618E+03_r8,     0.2456274E+03_r8,&  
         0.2441907E+03_r8,     0.2427525E+03_r8,     0.2413133E+03_r8,     0.2398739E+03_r8,     0.2384352E+03_r8,&  
         0.2369975E+03_r8,     0.2355616E+03_r8,     0.2341281E+03_r8,     0.2326974E+03_r8,     0.2312702E+03_r8,&  
         0.2298469E+03_r8,     0.2284280E+03_r8,     0.2270138E+03_r8,     0.2256048E+03_r8,     0.2242013E+03_r8,&  
         0.2228036E+03_r8,     0.2214120E+03_r8,     0.2200268E+03_r8,     0.2186481E+03_r8,     0.2172762E+03_r8,&  
         0.2159112E+03_r8,     0.2145534E+03_r8,     0.2132027E+03_r8,     0.2118595E+03_r8,     0.2105236E+03_r8,&  
         0.2091953E+03_r8,     0.2078746E+03_r8,     0.2065616E+03_r8,     0.2052562E+03_r8,     0.2039586E+03_r8,&  
         0.2026687E+03_r8,     0.2013866E+03_r8,     0.2001122E+03_r8,     0.1988456E+03_r8,     0.1975867E+03_r8,&  
         0.1963356E+03_r8,     0.1950922E+03_r8,     0.1938565E+03_r8,     0.1926285E+03_r8,     0.1914082E+03_r8,&  
         0.1901954E+03_r8,     0.1889903E+03_r8,     0.1877928E+03_r8,     0.1866027E+03_r8,     0.1854202E+03_r8,&  
         0.1842450E+03_r8,     0.1830773E+03_r8,     0.1819170E+03_r8,     0.1807640E+03_r8,     0.1796182E+03_r8,&  
         0.1784798E+03_r8,     0.1773485E+03_r8,     0.1762243E+03_r8,     0.1751073E+03_r8,     0.1739973E+03_r8,&  
         0.1728944E+03_r8,     0.1717985E+03_r8,     0.1707095E+03_r8,     0.1696274E+03_r8,     0.1685521E+03_r8,&  
         0.1674837E+03_r8,     0.1664220E+03_r8,     0.1653671E+03_r8,     0.1643188E+03_r8,     0.1632772E+03_r8,&  
         0.1622422E+03_r8,     0.1612138E+03_r8,     0.1601918E+03_r8,     0.1591764E+03_r8,     0.1581674E+03_r8,&  
         0.1571648E+03_r8,     0.1561685E+03_r8,     0.1551786E+03_r8,     0.1541949E+03_r8,     0.1532175E+03_r8,&  
         0.1522462E+03_r8,     0.1512811E+03_r8,     0.1503222E+03_r8,     0.1493693E+03_r8,     0.1484224E+03_r8,&  
         0.1474816E+03_r8,     0.1465467E+03_r8,     0.1456177E+03_r8,     0.1446947E+03_r8,     0.1437775E+03_r8,&  
         0.1428661E+03_r8,     0.1419604E+03_r8,     0.1410606E+03_r8,     0.1401664E+03_r8,     0.1392779E+03_r8,&  
         0.1383950E+03_r8,     0.1375177E+03_r8,     0.1366460E+03_r8,     0.1357798E+03_r8,     0.1349191E+03_r8,&  
         0.1340638E+03_r8,     0.1332140E+03_r8,     0.1323696E+03_r8,     0.1315305E+03_r8,     0.1306967E+03_r8,&  
         0.1298682E+03_r8,     0.1290450E+03_r8,     0.1282270E+03_r8,     0.1274141E+03_r8,     0.1266065E+03_r8,&  
         0.1258039E+03_r8,     0.1250065E+03_r8,     0.1242140E+03_r8,     0.1234267E+03_r8,     0.1226443E+03_r8,&  
         0.1218668E+03_r8,     0.1210943E+03_r8,     0.1203267E+03_r8,     0.1195639E+03_r8,     0.1188060E+03_r8,&  
         0.1180529E+03_r8,     0.1173046E+03_r8,     0.1165610E+03_r8,     0.1158221E+03_r8,     0.1150879E+03_r8,&  
         0.1143584E+03_r8,     0.1136335E+03_r8,     0.1129131E+03_r8,     0.1121974E+03_r8,     0.1114862E+03_r8,&  
         0.1107795E+03_r8,     0.1100772E+03_r8,     0.1093795E+03_r8,     0.1086861E+03_r8,     0.1079971E+03_r8,&  
         0.1073126E+03_r8,     0.1066323E+03_r8,     0.1059564E+03_r8,     0.1052847E+03_r8,     0.1046173E+03_r8,&  
         0.1039541E+03_r8,     0.1032952E+03_r8,     0.1026404E+03_r8,     0.1019898E+03_r8,     0.1013432E+03_r8/)  

    psaditmk(1:150,   31)= (/ &
         0.2600871E+03_r8,     0.2587141E+03_r8,     0.2573321E+03_r8,     0.2559418E+03_r8,     0.2545437E+03_r8,&  
         0.2531384E+03_r8,     0.2517269E+03_r8,     0.2503095E+03_r8,     0.2488872E+03_r8,     0.2474605E+03_r8,&  
         0.2460304E+03_r8,     0.2445975E+03_r8,     0.2431626E+03_r8,     0.2417263E+03_r8,     0.2402894E+03_r8,&  
         0.2388527E+03_r8,     0.2374167E+03_r8,     0.2359821E+03_r8,     0.2345495E+03_r8,     0.2331195E+03_r8,&  
         0.2316926E+03_r8,     0.2302694E+03_r8,     0.2288502E+03_r8,     0.2274356E+03_r8,     0.2260259E+03_r8,&  
         0.2246215E+03_r8,     0.2232227E+03_r8,     0.2218298E+03_r8,     0.2204432E+03_r8,     0.2190630E+03_r8,&  
         0.2176895E+03_r8,     0.2163228E+03_r8,     0.2149631E+03_r8,     0.2136105E+03_r8,     0.2122653E+03_r8,&  
         0.2109274E+03_r8,     0.2095970E+03_r8,     0.2082742E+03_r8,     0.2069589E+03_r8,     0.2056513E+03_r8,&  
         0.2043515E+03_r8,     0.2030593E+03_r8,     0.2017749E+03_r8,     0.2004982E+03_r8,     0.1992293E+03_r8,&  
         0.1979682E+03_r8,     0.1967148E+03_r8,     0.1954691E+03_r8,     0.1942311E+03_r8,     0.1930008E+03_r8,&  
         0.1917781E+03_r8,     0.1905631E+03_r8,     0.1893557E+03_r8,     0.1881558E+03_r8,     0.1869635E+03_r8,&  
         0.1857787E+03_r8,     0.1846013E+03_r8,     0.1834314E+03_r8,     0.1822688E+03_r8,     0.1811136E+03_r8,&  
         0.1799657E+03_r8,     0.1788250E+03_r8,     0.1776915E+03_r8,     0.1765652E+03_r8,     0.1754460E+03_r8,&  
         0.1743339E+03_r8,     0.1732289E+03_r8,     0.1721308E+03_r8,     0.1710397E+03_r8,     0.1699555E+03_r8,&  
         0.1688782E+03_r8,     0.1678077E+03_r8,     0.1667440E+03_r8,     0.1656870E+03_r8,     0.1646367E+03_r8,&  
         0.1635931E+03_r8,     0.1625561E+03_r8,     0.1615256E+03_r8,     0.1605017E+03_r8,     0.1594843E+03_r8,&  
         0.1584733E+03_r8,     0.1574688E+03_r8,     0.1564706E+03_r8,     0.1554787E+03_r8,     0.1544932E+03_r8,&  
         0.1535138E+03_r8,     0.1525407E+03_r8,     0.1515738E+03_r8,     0.1506129E+03_r8,     0.1496582E+03_r8,&  
         0.1487095E+03_r8,     0.1477669E+03_r8,     0.1468302E+03_r8,     0.1458994E+03_r8,     0.1449746E+03_r8,&  
         0.1440556E+03_r8,     0.1431424E+03_r8,     0.1422351E+03_r8,     0.1413334E+03_r8,     0.1404375E+03_r8,&  
         0.1395473E+03_r8,     0.1386627E+03_r8,     0.1377837E+03_r8,     0.1369103E+03_r8,     0.1360424E+03_r8,&  
         0.1351801E+03_r8,     0.1343232E+03_r8,     0.1334717E+03_r8,     0.1326256E+03_r8,     0.1317849E+03_r8,&  
         0.1309495E+03_r8,     0.1301194E+03_r8,     0.1292946E+03_r8,     0.1284750E+03_r8,     0.1276606E+03_r8,&  
         0.1268514E+03_r8,     0.1260473E+03_r8,     0.1252483E+03_r8,     0.1244543E+03_r8,     0.1236654E+03_r8,&  
         0.1228815E+03_r8,     0.1221026E+03_r8,     0.1213285E+03_r8,     0.1205594E+03_r8,     0.1197952E+03_r8,&  
         0.1190358E+03_r8,     0.1182813E+03_r8,     0.1175315E+03_r8,     0.1167865E+03_r8,     0.1160462E+03_r8,&  
         0.1153106E+03_r8,     0.1145796E+03_r8,     0.1138533E+03_r8,     0.1131316E+03_r8,     0.1124144E+03_r8,&  
         0.1117018E+03_r8,     0.1109938E+03_r8,     0.1102902E+03_r8,     0.1095910E+03_r8,     0.1088963E+03_r8,&  
         0.1082061E+03_r8,     0.1075201E+03_r8,     0.1068386E+03_r8,     0.1061613E+03_r8,     0.1054884E+03_r8,&  
         0.1048197E+03_r8,     0.1041552E+03_r8,     0.1034950E+03_r8,     0.1028389E+03_r8,     0.1021870E+03_r8/)  

    psaditmk(1:150,   32)= (/ &
         0.2617758E+03_r8,     0.2604211E+03_r8,     0.2590566E+03_r8,     0.2576828E+03_r8,     0.2563001E+03_r8,&  
         0.2549094E+03_r8,     0.2535110E+03_r8,     0.2521058E+03_r8,     0.2506944E+03_r8,     0.2492775E+03_r8,&  
         0.2478559E+03_r8,     0.2464304E+03_r8,     0.2450015E+03_r8,     0.2435701E+03_r8,     0.2421371E+03_r8,&  
         0.2407029E+03_r8,     0.2392685E+03_r8,     0.2378344E+03_r8,     0.2364013E+03_r8,     0.2349698E+03_r8,&  
         0.2335407E+03_r8,     0.2321143E+03_r8,     0.2306912E+03_r8,     0.2292720E+03_r8,     0.2278570E+03_r8,&  
         0.2264467E+03_r8,     0.2250416E+03_r8,     0.2236418E+03_r8,     0.2222478E+03_r8,     0.2208599E+03_r8,&  
         0.2194782E+03_r8,     0.2181031E+03_r8,     0.2167347E+03_r8,     0.2153732E+03_r8,     0.2140188E+03_r8,&  
         0.2126716E+03_r8,     0.2113317E+03_r8,     0.2099993E+03_r8,     0.2086743E+03_r8,     0.2073569E+03_r8,&  
         0.2060471E+03_r8,     0.2047450E+03_r8,     0.2034506E+03_r8,     0.2021639E+03_r8,     0.2008850E+03_r8,&  
         0.1996138E+03_r8,     0.1983503E+03_r8,     0.1970946E+03_r8,     0.1958466E+03_r8,     0.1946063E+03_r8,&  
         0.1933737E+03_r8,     0.1921488E+03_r8,     0.1909314E+03_r8,     0.1897217E+03_r8,     0.1885196E+03_r8,&  
         0.1873250E+03_r8,     0.1861379E+03_r8,     0.1849583E+03_r8,     0.1837861E+03_r8,     0.1826213E+03_r8,&  
         0.1814639E+03_r8,     0.1803137E+03_r8,     0.1791709E+03_r8,     0.1780352E+03_r8,     0.1769067E+03_r8,&  
         0.1757854E+03_r8,     0.1746711E+03_r8,     0.1735639E+03_r8,     0.1724638E+03_r8,     0.1713706E+03_r8,&  
         0.1702843E+03_r8,     0.1692048E+03_r8,     0.1681323E+03_r8,     0.1670665E+03_r8,     0.1660075E+03_r8,&  
         0.1649552E+03_r8,     0.1639096E+03_r8,     0.1628705E+03_r8,     0.1618381E+03_r8,     0.1608122E+03_r8,&  
         0.1597928E+03_r8,     0.1587799E+03_r8,     0.1577734E+03_r8,     0.1567733E+03_r8,     0.1557795E+03_r8,&  
         0.1547920E+03_r8,     0.1538108E+03_r8,     0.1528358E+03_r8,     0.1518670E+03_r8,     0.1509043E+03_r8,&  
         0.1499477E+03_r8,     0.1489972E+03_r8,     0.1480527E+03_r8,     0.1471142E+03_r8,     0.1461817E+03_r8,&  
         0.1452550E+03_r8,     0.1443343E+03_r8,     0.1434193E+03_r8,     0.1425102E+03_r8,     0.1416068E+03_r8,&  
         0.1407092E+03_r8,     0.1398172E+03_r8,     0.1389309E+03_r8,     0.1380503E+03_r8,     0.1371752E+03_r8,&  
         0.1363056E+03_r8,     0.1354416E+03_r8,     0.1345830E+03_r8,     0.1337299E+03_r8,     0.1328822E+03_r8,&  
         0.1320399E+03_r8,     0.1312029E+03_r8,     0.1303712E+03_r8,     0.1295447E+03_r8,     0.1287236E+03_r8,&  
         0.1279076E+03_r8,     0.1270968E+03_r8,     0.1262911E+03_r8,     0.1254906E+03_r8,     0.1246951E+03_r8,&  
         0.1239046E+03_r8,     0.1231192E+03_r8,     0.1223388E+03_r8,     0.1215632E+03_r8,     0.1207927E+03_r8,&  
         0.1200270E+03_r8,     0.1192661E+03_r8,     0.1185101E+03_r8,     0.1177589E+03_r8,     0.1170124E+03_r8,&  
         0.1162707E+03_r8,     0.1155336E+03_r8,     0.1148012E+03_r8,     0.1140735E+03_r8,     0.1133504E+03_r8,&  
         0.1126319E+03_r8,     0.1119179E+03_r8,     0.1112085E+03_r8,     0.1105035E+03_r8,     0.1098030E+03_r8,&  
         0.1091070E+03_r8,     0.1084154E+03_r8,     0.1077281E+03_r8,     0.1070453E+03_r8,     0.1063667E+03_r8,&  
         0.1056924E+03_r8,     0.1050225E+03_r8,     0.1043567E+03_r8,     0.1036952E+03_r8,     0.1030379E+03_r8/)  

    psaditmk(1:150,   33)= (/ &
         0.2634300E+03_r8,     0.2620948E+03_r8,     0.2607492E+03_r8,     0.2593935E+03_r8,     0.2580279E+03_r8,&  
         0.2566534E+03_r8,     0.2552701E+03_r8,     0.2538790E+03_r8,     0.2524804E+03_r8,     0.2510752E+03_r8,&  
         0.2496641E+03_r8,     0.2482478E+03_r8,     0.2468271E+03_r8,     0.2454026E+03_r8,     0.2439751E+03_r8,&  
         0.2425454E+03_r8,     0.2411142E+03_r8,     0.2396824E+03_r8,     0.2382504E+03_r8,     0.2368190E+03_r8,&  
         0.2353889E+03_r8,     0.2339608E+03_r8,     0.2325350E+03_r8,     0.2311123E+03_r8,     0.2296931E+03_r8,&  
         0.2282780E+03_r8,     0.2268673E+03_r8,     0.2254614E+03_r8,     0.2240608E+03_r8,     0.2226657E+03_r8,&  
         0.2212766E+03_r8,     0.2198936E+03_r8,     0.2185169E+03_r8,     0.2171469E+03_r8,     0.2157837E+03_r8,&  
         0.2144275E+03_r8,     0.2130784E+03_r8,     0.2117365E+03_r8,     0.2104020E+03_r8,     0.2090749E+03_r8,&  
         0.2077554E+03_r8,     0.2064435E+03_r8,     0.2051392E+03_r8,     0.2038425E+03_r8,     0.2025536E+03_r8,&  
         0.2012724E+03_r8,     0.1999989E+03_r8,     0.1987331E+03_r8,     0.1974751E+03_r8,     0.1962248E+03_r8,&  
         0.1949822E+03_r8,     0.1937473E+03_r8,     0.1925201E+03_r8,     0.1913004E+03_r8,     0.1900885E+03_r8,&  
         0.1888840E+03_r8,     0.1876872E+03_r8,     0.1864978E+03_r8,     0.1853160E+03_r8,     0.1841415E+03_r8,&  
         0.1829745E+03_r8,     0.1818148E+03_r8,     0.1806625E+03_r8,     0.1795174E+03_r8,     0.1783795E+03_r8,&  
         0.1772489E+03_r8,     0.1761254E+03_r8,     0.1750090E+03_r8,     0.1738997E+03_r8,     0.1727974E+03_r8,&  
         0.1717020E+03_r8,     0.1706136E+03_r8,     0.1695322E+03_r8,     0.1684575E+03_r8,     0.1673897E+03_r8,&  
         0.1663286E+03_r8,     0.1652743E+03_r8,     0.1642266E+03_r8,     0.1631856E+03_r8,     0.1621512E+03_r8,&  
         0.1611233E+03_r8,     0.1601020E+03_r8,     0.1590871E+03_r8,     0.1580786E+03_r8,     0.1570766E+03_r8,&  
         0.1560809E+03_r8,     0.1550915E+03_r8,     0.1541084E+03_r8,     0.1531315E+03_r8,     0.1521608E+03_r8,&  
         0.1511962E+03_r8,     0.1502378E+03_r8,     0.1492854E+03_r8,     0.1483391E+03_r8,     0.1473988E+03_r8,&  
         0.1464644E+03_r8,     0.1455360E+03_r8,     0.1446135E+03_r8,     0.1436967E+03_r8,     0.1427859E+03_r8,&  
         0.1418808E+03_r8,     0.1409814E+03_r8,     0.1400877E+03_r8,     0.1391997E+03_r8,     0.1383173E+03_r8,&  
         0.1374405E+03_r8,     0.1365693E+03_r8,     0.1357036E+03_r8,     0.1348434E+03_r8,     0.1339886E+03_r8,&  
         0.1331392E+03_r8,     0.1322953E+03_r8,     0.1314566E+03_r8,     0.1306234E+03_r8,     0.1297953E+03_r8,&  
         0.1289726E+03_r8,     0.1281550E+03_r8,     0.1273426E+03_r8,     0.1265354E+03_r8,     0.1257333E+03_r8,&  
         0.1249363E+03_r8,     0.1241443E+03_r8,     0.1233574E+03_r8,     0.1225754E+03_r8,     0.1217984E+03_r8,&  
         0.1210263E+03_r8,     0.1202591E+03_r8,     0.1194968E+03_r8,     0.1187393E+03_r8,     0.1179866E+03_r8,&  
         0.1172388E+03_r8,     0.1164956E+03_r8,     0.1157571E+03_r8,     0.1150233E+03_r8,     0.1142942E+03_r8,&  
         0.1135697E+03_r8,     0.1128498E+03_r8,     0.1121344E+03_r8,     0.1114236E+03_r8,     0.1107173E+03_r8,&  
         0.1100155E+03_r8,     0.1093181E+03_r8,     0.1086251E+03_r8,     0.1079365E+03_r8,     0.1072523E+03_r8,&  
         0.1065725E+03_r8,     0.1058969E+03_r8,     0.1052256E+03_r8,     0.1045586E+03_r8,     0.1038958E+03_r8/)  

    psaditmk(1:150,   34)= (/ &
         0.2650474E+03_r8,     0.2637332E+03_r8,     0.2624078E+03_r8,     0.2610713E+03_r8,     0.2597245E+03_r8,&  
         0.2583677E+03_r8,     0.2570013E+03_r8,     0.2556258E+03_r8,     0.2542421E+03_r8,     0.2528506E+03_r8,&  
         0.2514519E+03_r8,     0.2500468E+03_r8,     0.2486361E+03_r8,     0.2472204E+03_r8,     0.2458005E+03_r8,&  
         0.2443772E+03_r8,     0.2429512E+03_r8,     0.2415233E+03_r8,     0.2400942E+03_r8,     0.2386646E+03_r8,&  
         0.2372352E+03_r8,     0.2358066E+03_r8,     0.2343797E+03_r8,     0.2329548E+03_r8,     0.2315326E+03_r8,&  
         0.2301137E+03_r8,     0.2286985E+03_r8,     0.2272874E+03_r8,     0.2258811E+03_r8,     0.2244797E+03_r8,&  
         0.2230837E+03_r8,     0.2216934E+03_r8,     0.2203091E+03_r8,     0.2189310E+03_r8,     0.2175594E+03_r8,&  
         0.2161945E+03_r8,     0.2148365E+03_r8,     0.2134856E+03_r8,     0.2121418E+03_r8,     0.2108053E+03_r8,&  
         0.2094761E+03_r8,     0.2081545E+03_r8,     0.2068404E+03_r8,     0.2055339E+03_r8,     0.2042351E+03_r8,&  
         0.2029439E+03_r8,     0.2016604E+03_r8,     0.2003846E+03_r8,     0.1991166E+03_r8,     0.1978563E+03_r8,&  
         0.1966037E+03_r8,     0.1953588E+03_r8,     0.1941216E+03_r8,     0.1928920E+03_r8,     0.1916702E+03_r8,&  
         0.1904559E+03_r8,     0.1892492E+03_r8,     0.1880500E+03_r8,     0.1868584E+03_r8,     0.1856743E+03_r8,&  
         0.1844976E+03_r8,     0.1833283E+03_r8,     0.1821664E+03_r8,     0.1810119E+03_r8,     0.1798646E+03_r8,&  
         0.1787245E+03_r8,     0.1775917E+03_r8,     0.1764660E+03_r8,     0.1753475E+03_r8,     0.1742360E+03_r8,&  
         0.1731316E+03_r8,     0.1720342E+03_r8,     0.1709437E+03_r8,     0.1698601E+03_r8,     0.1687834E+03_r8,&  
         0.1677135E+03_r8,     0.1666504E+03_r8,     0.1655940E+03_r8,     0.1645443E+03_r8,     0.1635013E+03_r8,&  
         0.1624648E+03_r8,     0.1614350E+03_r8,     0.1604117E+03_r8,     0.1593948E+03_r8,     0.1583844E+03_r8,&  
         0.1573804E+03_r8,     0.1563828E+03_r8,     0.1553915E+03_r8,     0.1544064E+03_r8,     0.1534277E+03_r8,&  
         0.1524551E+03_r8,     0.1514887E+03_r8,     0.1505284E+03_r8,     0.1495742E+03_r8,     0.1486261E+03_r8,&  
         0.1476839E+03_r8,     0.1467478E+03_r8,     0.1458176E+03_r8,     0.1448932E+03_r8,     0.1439747E+03_r8,&  
         0.1430621E+03_r8,     0.1421552E+03_r8,     0.1412541E+03_r8,     0.1403587E+03_r8,     0.1394690E+03_r8,&  
         0.1385849E+03_r8,     0.1377064E+03_r8,     0.1368335E+03_r8,     0.1359661E+03_r8,     0.1351042E+03_r8,&  
         0.1342478E+03_r8,     0.1333968E+03_r8,     0.1325512E+03_r8,     0.1317109E+03_r8,     0.1308760E+03_r8,&  
         0.1300464E+03_r8,     0.1292220E+03_r8,     0.1284029E+03_r8,     0.1275890E+03_r8,     0.1267802E+03_r8,&  
         0.1259765E+03_r8,     0.1251780E+03_r8,     0.1243845E+03_r8,     0.1235960E+03_r8,     0.1228125E+03_r8,&  
         0.1220340E+03_r8,     0.1212604E+03_r8,     0.1204918E+03_r8,     0.1197280E+03_r8,     0.1189690E+03_r8,&  
         0.1182149E+03_r8,     0.1174655E+03_r8,     0.1167209E+03_r8,     0.1159810E+03_r8,     0.1152458E+03_r8,&  
         0.1145153E+03_r8,     0.1137894E+03_r8,     0.1130681E+03_r8,     0.1123513E+03_r8,     0.1116391E+03_r8,&  
         0.1109315E+03_r8,     0.1102283E+03_r8,     0.1095295E+03_r8,     0.1088352E+03_r8,     0.1081453E+03_r8,&  
         0.1074598E+03_r8,     0.1067786E+03_r8,     0.1061017E+03_r8,     0.1054292E+03_r8,     0.1047608E+03_r8/)  

    psaditmk(1:150,   35)= (/ &
         0.2666269E+03_r8,     0.2653343E+03_r8,     0.2640300E+03_r8,     0.2627142E+03_r8,     0.2613873E+03_r8,&  
         0.2600497E+03_r8,     0.2587016E+03_r8,     0.2573437E+03_r8,     0.2559765E+03_r8,     0.2546003E+03_r8,&  
         0.2532160E+03_r8,     0.2518242E+03_r8,     0.2504254E+03_r8,     0.2490205E+03_r8,     0.2476102E+03_r8,&  
         0.2461953E+03_r8,     0.2447764E+03_r8,     0.2433543E+03_r8,     0.2419299E+03_r8,     0.2405038E+03_r8,&  
         0.2390768E+03_r8,     0.2376496E+03_r8,     0.2362229E+03_r8,     0.2347973E+03_r8,     0.2333735E+03_r8,&  
         0.2319520E+03_r8,     0.2305335E+03_r8,     0.2291184E+03_r8,     0.2277072E+03_r8,     0.2263004E+03_r8,&  
         0.2248984E+03_r8,     0.2235015E+03_r8,     0.2221102E+03_r8,     0.2207247E+03_r8,     0.2193453E+03_r8,&  
         0.2179722E+03_r8,     0.2166057E+03_r8,     0.2152460E+03_r8,     0.2138932E+03_r8,     0.2125475E+03_r8,&  
         0.2112090E+03_r8,     0.2098779E+03_r8,     0.2085542E+03_r8,     0.2072379E+03_r8,     0.2059293E+03_r8,&  
         0.2046282E+03_r8,     0.2033348E+03_r8,     0.2020491E+03_r8,     0.2007711E+03_r8,     0.1995008E+03_r8,&  
         0.1982381E+03_r8,     0.1969833E+03_r8,     0.1957361E+03_r8,     0.1944966E+03_r8,     0.1932647E+03_r8,&  
         0.1920406E+03_r8,     0.1908240E+03_r8,     0.1896150E+03_r8,     0.1884136E+03_r8,     0.1872197E+03_r8,&  
         0.1860333E+03_r8,     0.1848544E+03_r8,     0.1836828E+03_r8,     0.1825187E+03_r8,     0.1813619E+03_r8,&  
         0.1802124E+03_r8,     0.1790702E+03_r8,     0.1779352E+03_r8,     0.1768074E+03_r8,     0.1756867E+03_r8,&  
         0.1745730E+03_r8,     0.1734665E+03_r8,     0.1723669E+03_r8,     0.1712743E+03_r8,     0.1701887E+03_r8,&  
         0.1691099E+03_r8,     0.1680379E+03_r8,     0.1669727E+03_r8,     0.1659143E+03_r8,     0.1648626E+03_r8,&  
         0.1638175E+03_r8,     0.1627791E+03_r8,     0.1617473E+03_r8,     0.1607219E+03_r8,     0.1597032E+03_r8,&  
         0.1586908E+03_r8,     0.1576848E+03_r8,     0.1566853E+03_r8,     0.1556920E+03_r8,     0.1547051E+03_r8,&  
         0.1537245E+03_r8,     0.1527500E+03_r8,     0.1517817E+03_r8,     0.1508196E+03_r8,     0.1498636E+03_r8,&  
         0.1489136E+03_r8,     0.1479696E+03_r8,     0.1470316E+03_r8,     0.1460996E+03_r8,     0.1451735E+03_r8,&  
         0.1442533E+03_r8,     0.1433388E+03_r8,     0.1424302E+03_r8,     0.1415274E+03_r8,     0.1406302E+03_r8,&  
         0.1397387E+03_r8,     0.1388530E+03_r8,     0.1379728E+03_r8,     0.1370982E+03_r8,     0.1362291E+03_r8,&  
         0.1353655E+03_r8,     0.1345075E+03_r8,     0.1336548E+03_r8,     0.1328076E+03_r8,     0.1319657E+03_r8,&  
         0.1311292E+03_r8,     0.1302980E+03_r8,     0.1294720E+03_r8,     0.1286513E+03_r8,     0.1278358E+03_r8,&  
         0.1270254E+03_r8,     0.1262202E+03_r8,     0.1254201E+03_r8,     0.1246251E+03_r8,     0.1238351E+03_r8,&  
         0.1230501E+03_r8,     0.1222701E+03_r8,     0.1214950E+03_r8,     0.1207249E+03_r8,     0.1199596E+03_r8,&  
         0.1191992E+03_r8,     0.1184436E+03_r8,     0.1176927E+03_r8,     0.1169467E+03_r8,     0.1162054E+03_r8,&  
         0.1154688E+03_r8,     0.1147368E+03_r8,     0.1140095E+03_r8,     0.1132868E+03_r8,     0.1125687E+03_r8,&  
         0.1118551E+03_r8,     0.1111460E+03_r8,     0.1104415E+03_r8,     0.1097414E+03_r8,     0.1090458E+03_r8,&  
         0.1083545E+03_r8,     0.1076677E+03_r8,     0.1069852E+03_r8,     0.1063070E+03_r8,     0.1056331E+03_r8/)  

    psaditmk(1:150,   36)= (/ &
         0.2681671E+03_r8,     0.2668968E+03_r8,     0.2656145E+03_r8,     0.2643203E+03_r8,     0.2630144E+03_r8,&  
         0.2616972E+03_r8,     0.2603688E+03_r8,     0.2590299E+03_r8,     0.2576805E+03_r8,     0.2563216E+03_r8,&  
         0.2549535E+03_r8,     0.2535766E+03_r8,     0.2521918E+03_r8,     0.2507998E+03_r8,     0.2494010E+03_r8,&  
         0.2479964E+03_r8,     0.2465866E+03_r8,     0.2451724E+03_r8,     0.2437546E+03_r8,     0.2423339E+03_r8,&  
         0.2409112E+03_r8,     0.2394870E+03_r8,     0.2380623E+03_r8,     0.2366376E+03_r8,     0.2352136E+03_r8,&  
         0.2337910E+03_r8,     0.2323705E+03_r8,     0.2309525E+03_r8,     0.2295377E+03_r8,     0.2281265E+03_r8,&  
         0.2267194E+03_r8,     0.2253169E+03_r8,     0.2239193E+03_r8,     0.2225270E+03_r8,     0.2211404E+03_r8,&  
         0.2197597E+03_r8,     0.2183852E+03_r8,     0.2170172E+03_r8,     0.2156558E+03_r8,     0.2143012E+03_r8,&  
         0.2129537E+03_r8,     0.2116133E+03_r8,     0.2102802E+03_r8,     0.2089544E+03_r8,     0.2076361E+03_r8,&  
         0.2063253E+03_r8,     0.2050220E+03_r8,     0.2037264E+03_r8,     0.2024385E+03_r8,     0.2011582E+03_r8,&  
         0.1998856E+03_r8,     0.1986207E+03_r8,     0.1973635E+03_r8,     0.1961140E+03_r8,     0.1948723E+03_r8,&  
         0.1936381E+03_r8,     0.1924116E+03_r8,     0.1911928E+03_r8,     0.1899815E+03_r8,     0.1887778E+03_r8,&  
         0.1875817E+03_r8,     0.1863930E+03_r8,     0.1852118E+03_r8,     0.1840381E+03_r8,     0.1828717E+03_r8,&  
         0.1817127E+03_r8,     0.1805610E+03_r8,     0.1794166E+03_r8,     0.1782794E+03_r8,     0.1771494E+03_r8,&  
         0.1760265E+03_r8,     0.1749108E+03_r8,     0.1738021E+03_r8,     0.1727004E+03_r8,     0.1716057E+03_r8,&  
         0.1705179E+03_r8,     0.1694370E+03_r8,     0.1683630E+03_r8,     0.1672957E+03_r8,     0.1662353E+03_r8,&  
         0.1651815E+03_r8,     0.1641344E+03_r8,     0.1630940E+03_r8,     0.1620602E+03_r8,     0.1610329E+03_r8,&  
         0.1600121E+03_r8,     0.1589978E+03_r8,     0.1579899E+03_r8,     0.1569884E+03_r8,     0.1559932E+03_r8,&  
         0.1550044E+03_r8,     0.1540219E+03_r8,     0.1530455E+03_r8,     0.1520754E+03_r8,     0.1511114E+03_r8,&  
         0.1501535E+03_r8,     0.1492017E+03_r8,     0.1482559E+03_r8,     0.1473161E+03_r8,     0.1463822E+03_r8,&  
         0.1454543E+03_r8,     0.1445323E+03_r8,     0.1436161E+03_r8,     0.1427057E+03_r8,     0.1418011E+03_r8,&  
         0.1409023E+03_r8,     0.1400091E+03_r8,     0.1391216E+03_r8,     0.1382397E+03_r8,     0.1373634E+03_r8,&  
         0.1364926E+03_r8,     0.1356274E+03_r8,     0.1347677E+03_r8,     0.1339134E+03_r8,     0.1330645E+03_r8,&  
         0.1322210E+03_r8,     0.1313829E+03_r8,     0.1305500E+03_r8,     0.1297225E+03_r8,     0.1289002E+03_r8,&  
         0.1280831E+03_r8,     0.1272712E+03_r8,     0.1264644E+03_r8,     0.1256627E+03_r8,     0.1248662E+03_r8,&  
         0.1240746E+03_r8,     0.1232881E+03_r8,     0.1225066E+03_r8,     0.1217300E+03_r8,     0.1209584E+03_r8,&  
         0.1201917E+03_r8,     0.1194298E+03_r8,     0.1186727E+03_r8,     0.1179204E+03_r8,     0.1171729E+03_r8,&  
         0.1164302E+03_r8,     0.1156921E+03_r8,     0.1149588E+03_r8,     0.1142300E+03_r8,     0.1135059E+03_r8,&  
         0.1127864E+03_r8,     0.1120715E+03_r8,     0.1113611E+03_r8,     0.1106551E+03_r8,     0.1099537E+03_r8,&  
         0.1092567E+03_r8,     0.1085641E+03_r8,     0.1078759E+03_r8,     0.1071921E+03_r8,     0.1065126E+03_r8/)  

    psaditmk(1:150,   37)= (/ &
         0.2696669E+03_r8,     0.2684195E+03_r8,     0.2671598E+03_r8,     0.2658879E+03_r8,     0.2646039E+03_r8,&  
         0.2633082E+03_r8,     0.2620007E+03_r8,     0.2606819E+03_r8,     0.2593521E+03_r8,     0.2580118E+03_r8,&  
         0.2566614E+03_r8,     0.2553013E+03_r8,     0.2539324E+03_r8,     0.2525548E+03_r8,     0.2511696E+03_r8,&  
         0.2497773E+03_r8,     0.2483786E+03_r8,     0.2469744E+03_r8,     0.2455652E+03_r8,     0.2441518E+03_r8,&  
         0.2427352E+03_r8,     0.2413160E+03_r8,     0.2398950E+03_r8,     0.2384729E+03_r8,     0.2370505E+03_r8,&  
         0.2356284E+03_r8,     0.2342072E+03_r8,     0.2327878E+03_r8,     0.2313706E+03_r8,     0.2299562E+03_r8,&  
         0.2285451E+03_r8,     0.2271379E+03_r8,     0.2257350E+03_r8,     0.2243368E+03_r8,     0.2229437E+03_r8,&  
         0.2215561E+03_r8,     0.2201742E+03_r8,     0.2187984E+03_r8,     0.2174289E+03_r8,     0.2160659E+03_r8,&  
         0.2147097E+03_r8,     0.2133603E+03_r8,     0.2120180E+03_r8,     0.2106829E+03_r8,     0.2093551E+03_r8,&  
         0.2080347E+03_r8,     0.2067218E+03_r8,     0.2054164E+03_r8,     0.2041186E+03_r8,     0.2028284E+03_r8,&  
         0.2015459E+03_r8,     0.2002711E+03_r8,     0.1990039E+03_r8,     0.1977444E+03_r8,     0.1964927E+03_r8,&  
         0.1952486E+03_r8,     0.1940121E+03_r8,     0.1927834E+03_r8,     0.1915622E+03_r8,     0.1903487E+03_r8,&  
         0.1891427E+03_r8,     0.1879443E+03_r8,     0.1867534E+03_r8,     0.1855699E+03_r8,     0.1843939E+03_r8,&  
         0.1832253E+03_r8,     0.1820641E+03_r8,     0.1809102E+03_r8,     0.1797636E+03_r8,     0.1786242E+03_r8,&  
         0.1774920E+03_r8,     0.1763670E+03_r8,     0.1752491E+03_r8,     0.1741382E+03_r8,     0.1730344E+03_r8,&  
         0.1719376E+03_r8,     0.1708477E+03_r8,     0.1697647E+03_r8,     0.1686886E+03_r8,     0.1676193E+03_r8,&  
         0.1665568E+03_r8,     0.1655010E+03_r8,     0.1644519E+03_r8,     0.1634095E+03_r8,     0.1623736E+03_r8,&  
         0.1613444E+03_r8,     0.1603216E+03_r8,     0.1593053E+03_r8,     0.1582955E+03_r8,     0.1572921E+03_r8,&  
         0.1562950E+03_r8,     0.1553043E+03_r8,     0.1543198E+03_r8,     0.1533416E+03_r8,     0.1523696E+03_r8,&  
         0.1514037E+03_r8,     0.1504439E+03_r8,     0.1494903E+03_r8,     0.1485426E+03_r8,     0.1476010E+03_r8,&  
         0.1466654E+03_r8,     0.1457357E+03_r8,     0.1448119E+03_r8,     0.1438939E+03_r8,     0.1429818E+03_r8,&  
         0.1420754E+03_r8,     0.1411748E+03_r8,     0.1402799E+03_r8,     0.1393907E+03_r8,     0.1385071E+03_r8,&  
         0.1376291E+03_r8,     0.1367567E+03_r8,     0.1358898E+03_r8,     0.1350284E+03_r8,     0.1341724E+03_r8,&  
         0.1333219E+03_r8,     0.1324768E+03_r8,     0.1316370E+03_r8,     0.1308026E+03_r8,     0.1299734E+03_r8,&  
         0.1291495E+03_r8,     0.1283308E+03_r8,     0.1275173E+03_r8,     0.1267090E+03_r8,     0.1259058E+03_r8,&  
         0.1251077E+03_r8,     0.1243146E+03_r8,     0.1235266E+03_r8,     0.1227436E+03_r8,     0.1219655E+03_r8,&  
         0.1211924E+03_r8,     0.1204241E+03_r8,     0.1196608E+03_r8,     0.1189023E+03_r8,     0.1181485E+03_r8,&  
         0.1173996E+03_r8,     0.1166554E+03_r8,     0.1159159E+03_r8,     0.1151811E+03_r8,     0.1144510E+03_r8,&  
         0.1137255E+03_r8,     0.1130046E+03_r8,     0.1122883E+03_r8,     0.1115765E+03_r8,     0.1108692E+03_r8,&  
         0.1101664E+03_r8,     0.1094680E+03_r8,     0.1087741E+03_r8,     0.1080846E+03_r8,     0.1073995E+03_r8/)  

    psaditmk(1:150,   38)= (/ &
         0.2711260E+03_r8,     0.2699018E+03_r8,     0.2686649E+03_r8,     0.2674160E+03_r8,     0.2661545E+03_r8,&  
         0.2648809E+03_r8,     0.2635952E+03_r8,     0.2622977E+03_r8,     0.2609886E+03_r8,     0.2596683E+03_r8,&  
         0.2583370E+03_r8,     0.2569955E+03_r8,     0.2556439E+03_r8,     0.2542828E+03_r8,     0.2529130E+03_r8,&  
         0.2515349E+03_r8,     0.2501492E+03_r8,     0.2487568E+03_r8,     0.2473584E+03_r8,     0.2459544E+03_r8,&  
         0.2445459E+03_r8,     0.2431337E+03_r8,     0.2417183E+03_r8,     0.2403007E+03_r8,     0.2388815E+03_r8,&  
         0.2374616E+03_r8,     0.2360415E+03_r8,     0.2346221E+03_r8,     0.2332039E+03_r8,     0.2317877E+03_r8,&  
         0.2303739E+03_r8,     0.2289631E+03_r8,     0.2275559E+03_r8,     0.2261528E+03_r8,     0.2247541E+03_r8,&  
         0.2233603E+03_r8,     0.2219717E+03_r8,     0.2205888E+03_r8,     0.2192117E+03_r8,     0.2178408E+03_r8,&  
         0.2164763E+03_r8,     0.2151184E+03_r8,     0.2137673E+03_r8,     0.2124232E+03_r8,     0.2110861E+03_r8,&  
         0.2097564E+03_r8,     0.2084339E+03_r8,     0.2071189E+03_r8,     0.2058114E+03_r8,     0.2045114E+03_r8,&  
         0.2032190E+03_r8,     0.2019343E+03_r8,     0.2006572E+03_r8,     0.1993877E+03_r8,     0.1981260E+03_r8,&  
         0.1968720E+03_r8,     0.1956256E+03_r8,     0.1943869E+03_r8,     0.1931558E+03_r8,     0.1919324E+03_r8,&  
         0.1907166E+03_r8,     0.1895083E+03_r8,     0.1883076E+03_r8,     0.1871144E+03_r8,     0.1859287E+03_r8,&  
         0.1847505E+03_r8,     0.1835797E+03_r8,     0.1824162E+03_r8,     0.1812601E+03_r8,     0.1801113E+03_r8,&  
         0.1789697E+03_r8,     0.1778353E+03_r8,     0.1767081E+03_r8,     0.1755880E+03_r8,     0.1744750E+03_r8,&  
         0.1733691E+03_r8,     0.1722702E+03_r8,     0.1711782E+03_r8,     0.1700931E+03_r8,     0.1690149E+03_r8,&  
         0.1679436E+03_r8,     0.1668790E+03_r8,     0.1658212E+03_r8,     0.1647700E+03_r8,     0.1637256E+03_r8,&  
         0.1626877E+03_r8,     0.1616565E+03_r8,     0.1606317E+03_r8,     0.1596135E+03_r8,     0.1586017E+03_r8,&  
         0.1575963E+03_r8,     0.1565974E+03_r8,     0.1556047E+03_r8,     0.1546183E+03_r8,     0.1536382E+03_r8,&  
         0.1526643E+03_r8,     0.1516966E+03_r8,     0.1507350E+03_r8,     0.1497794E+03_r8,     0.1488300E+03_r8,&  
         0.1478866E+03_r8,     0.1469491E+03_r8,     0.1460176E+03_r8,     0.1450920E+03_r8,     0.1441723E+03_r8,&  
         0.1432584E+03_r8,     0.1423503E+03_r8,     0.1414479E+03_r8,     0.1405513E+03_r8,     0.1396603E+03_r8,&  
         0.1387750E+03_r8,     0.1378953E+03_r8,     0.1370212E+03_r8,     0.1361526E+03_r8,     0.1352896E+03_r8,&  
         0.1344320E+03_r8,     0.1335798E+03_r8,     0.1327330E+03_r8,     0.1318917E+03_r8,     0.1310556E+03_r8,&  
         0.1302248E+03_r8,     0.1293994E+03_r8,     0.1285791E+03_r8,     0.1277640E+03_r8,     0.1269541E+03_r8,&  
         0.1261494E+03_r8,     0.1253497E+03_r8,     0.1245551E+03_r8,     0.1237656E+03_r8,     0.1229810E+03_r8,&  
         0.1222015E+03_r8,     0.1214268E+03_r8,     0.1206571E+03_r8,     0.1198923E+03_r8,     0.1191323E+03_r8,&  
         0.1183771E+03_r8,     0.1176267E+03_r8,     0.1168811E+03_r8,     0.1161402E+03_r8,     0.1154039E+03_r8,&  
         0.1146724E+03_r8,     0.1139455E+03_r8,     0.1132232E+03_r8,     0.1125055E+03_r8,     0.1117923E+03_r8,&  
         0.1110837E+03_r8,     0.1103795E+03_r8,     0.1096798E+03_r8,     0.1089846E+03_r8,     0.1082937E+03_r8/)  

    psaditmk(1:150,   39)= (/ &
         0.2725441E+03_r8,     0.2713430E+03_r8,     0.2701295E+03_r8,     0.2689036E+03_r8,     0.2676651E+03_r8,&  
         0.2664143E+03_r8,     0.2651512E+03_r8,     0.2638758E+03_r8,     0.2625883E+03_r8,     0.2612892E+03_r8,&  
         0.2599785E+03_r8,     0.2586566E+03_r8,     0.2573239E+03_r8,     0.2559809E+03_r8,     0.2546283E+03_r8,&  
         0.2532661E+03_r8,     0.2518956E+03_r8,     0.2505170E+03_r8,     0.2491310E+03_r8,     0.2477385E+03_r8,&  
         0.2463402E+03_r8,     0.2449368E+03_r8,     0.2435291E+03_r8,     0.2421179E+03_r8,     0.2407039E+03_r8,&  
         0.2392879E+03_r8,     0.2378707E+03_r8,     0.2364530E+03_r8,     0.2350355E+03_r8,     0.2336189E+03_r8,&  
         0.2322037E+03_r8,     0.2307908E+03_r8,     0.2293805E+03_r8,     0.2279735E+03_r8,     0.2265702E+03_r8,&  
         0.2251712E+03_r8,     0.2237768E+03_r8,     0.2223875E+03_r8,     0.2210035E+03_r8,     0.2196252E+03_r8,&  
         0.2182530E+03_r8,     0.2168870E+03_r8,     0.2155276E+03_r8,     0.2141748E+03_r8,     0.2128288E+03_r8,&  
         0.2114900E+03_r8,     0.2101582E+03_r8,     0.2088338E+03_r8,     0.2075166E+03_r8,     0.2062070E+03_r8,&  
         0.2049049E+03_r8,     0.2036103E+03_r8,     0.2023234E+03_r8,     0.2010440E+03_r8,     0.1997724E+03_r8,&  
         0.1985084E+03_r8,     0.1972520E+03_r8,     0.1960034E+03_r8,     0.1947624E+03_r8,     0.1935290E+03_r8,&  
         0.1923033E+03_r8,     0.1910852E+03_r8,     0.1898747E+03_r8,     0.1886717E+03_r8,     0.1874763E+03_r8,&  
         0.1862883E+03_r8,     0.1851078E+03_r8,     0.1839348E+03_r8,     0.1827690E+03_r8,     0.1816107E+03_r8,&  
         0.1804597E+03_r8,     0.1793159E+03_r8,     0.1781794E+03_r8,     0.1770500E+03_r8,     0.1759278E+03_r8,&  
         0.1748126E+03_r8,     0.1737046E+03_r8,     0.1726035E+03_r8,     0.1715094E+03_r8,     0.1704222E+03_r8,&  
         0.1693419E+03_r8,     0.1682685E+03_r8,     0.1672019E+03_r8,     0.1661420E+03_r8,     0.1650889E+03_r8,&  
         0.1640424E+03_r8,     0.1630025E+03_r8,     0.1619693E+03_r8,     0.1609426E+03_r8,     0.1599223E+03_r8,&  
         0.1589086E+03_r8,     0.1579013E+03_r8,     0.1569003E+03_r8,     0.1559058E+03_r8,     0.1549175E+03_r8,&  
         0.1539355E+03_r8,     0.1529597E+03_r8,     0.1519901E+03_r8,     0.1510266E+03_r8,     0.1500693E+03_r8,&  
         0.1491180E+03_r8,     0.1481727E+03_r8,     0.1472334E+03_r8,     0.1463002E+03_r8,     0.1453728E+03_r8,&  
         0.1444512E+03_r8,     0.1435356E+03_r8,     0.1426257E+03_r8,     0.1417216E+03_r8,     0.1408232E+03_r8,&  
         0.1399306E+03_r8,     0.1390435E+03_r8,     0.1381621E+03_r8,     0.1372863E+03_r8,     0.1364161E+03_r8,&  
         0.1355513E+03_r8,     0.1346921E+03_r8,     0.1338383E+03_r8,     0.1329899E+03_r8,     0.1321469E+03_r8,&  
         0.1313092E+03_r8,     0.1304768E+03_r8,     0.1296497E+03_r8,     0.1288279E+03_r8,     0.1280112E+03_r8,&  
         0.1271998E+03_r8,     0.1263935E+03_r8,     0.1255923E+03_r8,     0.1247961E+03_r8,     0.1240050E+03_r8,&  
         0.1232190E+03_r8,     0.1224379E+03_r8,     0.1216618E+03_r8,     0.1208906E+03_r8,     0.1201242E+03_r8,&  
         0.1193628E+03_r8,     0.1186061E+03_r8,     0.1178543E+03_r8,     0.1171072E+03_r8,     0.1163649E+03_r8,&  
         0.1156273E+03_r8,     0.1148943E+03_r8,     0.1141660E+03_r8,     0.1134423E+03_r8,     0.1127232E+03_r8,&  
         0.1120086E+03_r8,     0.1112986E+03_r8,     0.1105931E+03_r8,     0.1098920E+03_r8,     0.1091954E+03_r8/)  

    psaditmk(1:150,   40)= (/ &
         0.2739214E+03_r8,     0.2727433E+03_r8,     0.2715530E+03_r8,     0.2703502E+03_r8,     0.2691350E+03_r8,&  
         0.2679073E+03_r8,     0.2666672E+03_r8,     0.2654145E+03_r8,     0.2641496E+03_r8,     0.2628723E+03_r8,&  
         0.2615832E+03_r8,     0.2602823E+03_r8,     0.2589699E+03_r8,     0.2576465E+03_r8,     0.2563123E+03_r8,&  
         0.2549680E+03_r8,     0.2536142E+03_r8,     0.2522512E+03_r8,     0.2508799E+03_r8,     0.2495008E+03_r8,&  
         0.2481146E+03_r8,     0.2467222E+03_r8,     0.2453241E+03_r8,     0.2439213E+03_r8,     0.2425145E+03_r8,&  
         0.2411044E+03_r8,     0.2396919E+03_r8,     0.2382777E+03_r8,     0.2368625E+03_r8,     0.2354472E+03_r8,&  
         0.2340322E+03_r8,     0.2326185E+03_r8,     0.2312065E+03_r8,     0.2297969E+03_r8,     0.2283902E+03_r8,&  
         0.2269870E+03_r8,     0.2255878E+03_r8,     0.2241929E+03_r8,     0.2228029E+03_r8,     0.2214180E+03_r8,&  
         0.2200387E+03_r8,     0.2186652E+03_r8,     0.2172979E+03_r8,     0.2159369E+03_r8,     0.2145825E+03_r8,&  
         0.2132348E+03_r8,     0.2118941E+03_r8,     0.2105604E+03_r8,     0.2092340E+03_r8,     0.2079149E+03_r8,&  
         0.2066031E+03_r8,     0.2052989E+03_r8,     0.2040022E+03_r8,     0.2027130E+03_r8,     0.2014315E+03_r8,&  
         0.2001575E+03_r8,     0.1988913E+03_r8,     0.1976327E+03_r8,     0.1963817E+03_r8,     0.1951385E+03_r8,&  
         0.1939028E+03_r8,     0.1926749E+03_r8,     0.1914544E+03_r8,     0.1902416E+03_r8,     0.1890364E+03_r8,&  
         0.1878387E+03_r8,     0.1866485E+03_r8,     0.1854657E+03_r8,     0.1842904E+03_r8,     0.1831225E+03_r8,&  
         0.1819619E+03_r8,     0.1808087E+03_r8,     0.1796627E+03_r8,     0.1785240E+03_r8,     0.1773924E+03_r8,&  
         0.1762681E+03_r8,     0.1751508E+03_r8,     0.1740405E+03_r8,     0.1729374E+03_r8,     0.1718412E+03_r8,&  
         0.1707519E+03_r8,     0.1696695E+03_r8,     0.1685940E+03_r8,     0.1675253E+03_r8,     0.1664634E+03_r8,&  
         0.1654082E+03_r8,     0.1643597E+03_r8,     0.1633178E+03_r8,     0.1622826E+03_r8,     0.1612539E+03_r8,&  
         0.1602317E+03_r8,     0.1592160E+03_r8,     0.1582067E+03_r8,     0.1572039E+03_r8,     0.1562074E+03_r8,&  
         0.1552171E+03_r8,     0.1542333E+03_r8,     0.1532556E+03_r8,     0.1522841E+03_r8,     0.1513188E+03_r8,&  
         0.1503596E+03_r8,     0.1494064E+03_r8,     0.1484594E+03_r8,     0.1475183E+03_r8,     0.1465832E+03_r8,&  
         0.1456540E+03_r8,     0.1447307E+03_r8,     0.1438132E+03_r8,     0.1429016E+03_r8,     0.1419958E+03_r8,&  
         0.1410956E+03_r8,     0.1402012E+03_r8,     0.1393125E+03_r8,     0.1384294E+03_r8,     0.1375519E+03_r8,&  
         0.1366800E+03_r8,     0.1358136E+03_r8,     0.1349526E+03_r8,     0.1340972E+03_r8,     0.1332471E+03_r8,&  
         0.1324025E+03_r8,     0.1315632E+03_r8,     0.1307292E+03_r8,     0.1299005E+03_r8,     0.1290771E+03_r8,&  
         0.1282589E+03_r8,     0.1274458E+03_r8,     0.1266380E+03_r8,     0.1258352E+03_r8,     0.1250375E+03_r8,&  
         0.1242449E+03_r8,     0.1234574E+03_r8,     0.1226748E+03_r8,     0.1218971E+03_r8,     0.1211244E+03_r8,&  
         0.1203566E+03_r8,     0.1195937E+03_r8,     0.1188356E+03_r8,     0.1180823E+03_r8,     0.1173338E+03_r8,&  
         0.1165900E+03_r8,     0.1158509E+03_r8,     0.1151166E+03_r8,     0.1143868E+03_r8,     0.1136617E+03_r8,&  
         0.1129412E+03_r8,     0.1122253E+03_r8,     0.1115139E+03_r8,     0.1108070E+03_r8,     0.1101046E+03_r8/)  

    psaditmk(1:150,   41)= (/ &
         0.2752579E+03_r8,     0.2741026E+03_r8,     0.2729354E+03_r8,     0.2717556E+03_r8,     0.2705638E+03_r8,&  
         0.2693593E+03_r8,     0.2681426E+03_r8,     0.2669130E+03_r8,     0.2656710E+03_r8,     0.2644165E+03_r8,&  
         0.2631498E+03_r8,     0.2618707E+03_r8,     0.2605797E+03_r8,     0.2592769E+03_r8,     0.2579629E+03_r8,&  
         0.2566379E+03_r8,     0.2553024E+03_r8,     0.2539568E+03_r8,     0.2526019E+03_r8,     0.2512380E+03_r8,&  
         0.2498661E+03_r8,     0.2484865E+03_r8,     0.2471001E+03_r8,     0.2457078E+03_r8,     0.2443101E+03_r8,&  
         0.2429080E+03_r8,     0.2415021E+03_r8,     0.2400933E+03_r8,     0.2386824E+03_r8,     0.2372700E+03_r8,&  
         0.2358570E+03_r8,     0.2344441E+03_r8,     0.2330319E+03_r8,     0.2316211E+03_r8,     0.2302123E+03_r8,&  
         0.2288062E+03_r8,     0.2274032E+03_r8,     0.2260039E+03_r8,     0.2246087E+03_r8,     0.2232181E+03_r8,&  
         0.2218325E+03_r8,     0.2204522E+03_r8,     0.2190776E+03_r8,     0.2177089E+03_r8,     0.2163465E+03_r8,&  
         0.2149905E+03_r8,     0.2136411E+03_r8,     0.2122986E+03_r8,     0.2109631E+03_r8,     0.2096347E+03_r8,&  
         0.2083136E+03_r8,     0.2069998E+03_r8,     0.2056935E+03_r8,     0.2043945E+03_r8,     0.2031032E+03_r8,&  
         0.2018195E+03_r8,     0.2005433E+03_r8,     0.1992748E+03_r8,     0.1980140E+03_r8,     0.1967608E+03_r8,&  
         0.1955152E+03_r8,     0.1942773E+03_r8,     0.1930470E+03_r8,     0.1918244E+03_r8,     0.1906093E+03_r8,&  
         0.1894017E+03_r8,     0.1882018E+03_r8,     0.1870093E+03_r8,     0.1858243E+03_r8,     0.1846468E+03_r8,&  
         0.1834766E+03_r8,     0.1823138E+03_r8,     0.1811584E+03_r8,     0.1800102E+03_r8,     0.1788692E+03_r8,&  
         0.1777355E+03_r8,     0.1766090E+03_r8,     0.1754895E+03_r8,     0.1743772E+03_r8,     0.1732719E+03_r8,&  
         0.1721735E+03_r8,     0.1710822E+03_r8,     0.1699977E+03_r8,     0.1689201E+03_r8,     0.1678494E+03_r8,&  
         0.1667854E+03_r8,     0.1657282E+03_r8,     0.1646776E+03_r8,     0.1636337E+03_r8,     0.1625965E+03_r8,&  
         0.1615658E+03_r8,     0.1605417E+03_r8,     0.1595240E+03_r8,     0.1585128E+03_r8,     0.1575080E+03_r8,&  
         0.1565095E+03_r8,     0.1555174E+03_r8,     0.1545316E+03_r8,     0.1535520E+03_r8,     0.1525787E+03_r8,&  
         0.1516115E+03_r8,     0.1506504E+03_r8,     0.1496954E+03_r8,     0.1487466E+03_r8,     0.1478036E+03_r8,&  
         0.1468667E+03_r8,     0.1459357E+03_r8,     0.1450106E+03_r8,     0.1440914E+03_r8,     0.1431780E+03_r8,&  
         0.1422704E+03_r8,     0.1413686E+03_r8,     0.1404725E+03_r8,     0.1395820E+03_r8,     0.1386972E+03_r8,&  
         0.1378180E+03_r8,     0.1369444E+03_r8,     0.1360763E+03_r8,     0.1352137E+03_r8,     0.1343566E+03_r8,&  
         0.1335049E+03_r8,     0.1326586E+03_r8,     0.1318177E+03_r8,     0.1309821E+03_r8,     0.1301518E+03_r8,&  
         0.1293268E+03_r8,     0.1285070E+03_r8,     0.1276924E+03_r8,     0.1268829E+03_r8,     0.1260786E+03_r8,&  
         0.1252794E+03_r8,     0.1244853E+03_r8,     0.1236962E+03_r8,     0.1229121E+03_r8,     0.1221329E+03_r8,&  
         0.1213587E+03_r8,     0.1205894E+03_r8,     0.1198250E+03_r8,     0.1190655E+03_r8,     0.1183107E+03_r8,&  
         0.1175607E+03_r8,     0.1168155E+03_r8,     0.1160750E+03_r8,     0.1153392E+03_r8,     0.1146081E+03_r8,&  
         0.1138816E+03_r8,     0.1131597E+03_r8,     0.1124424E+03_r8,     0.1117296E+03_r8,     0.1110214E+03_r8/)  

    psaditmk(1:150,   42)= (/ &
         0.2765544E+03_r8,     0.2754214E+03_r8,     0.2742769E+03_r8,     0.2731202E+03_r8,     0.2719514E+03_r8,&  
         0.2707703E+03_r8,     0.2695768E+03_r8,     0.2683706E+03_r8,     0.2671519E+03_r8,     0.2659205E+03_r8,&  
         0.2646767E+03_r8,     0.2634203E+03_r8,     0.2621515E+03_r8,     0.2608706E+03_r8,     0.2595779E+03_r8,&  
         0.2582733E+03_r8,     0.2569576E+03_r8,     0.2556310E+03_r8,     0.2542941E+03_r8,     0.2529473E+03_r8,&  
         0.2515912E+03_r8,     0.2502265E+03_r8,     0.2488539E+03_r8,     0.2474740E+03_r8,     0.2460876E+03_r8,&  
         0.2446954E+03_r8,     0.2432982E+03_r8,     0.2418968E+03_r8,     0.2404921E+03_r8,     0.2390846E+03_r8,&  
         0.2376754E+03_r8,     0.2362650E+03_r8,     0.2348542E+03_r8,     0.2334438E+03_r8,     0.2320345E+03_r8,&  
         0.2306267E+03_r8,     0.2292213E+03_r8,     0.2278187E+03_r8,     0.2264194E+03_r8,     0.2250240E+03_r8,&  
         0.2236331E+03_r8,     0.2222467E+03_r8,     0.2208656E+03_r8,     0.2194899E+03_r8,     0.2181200E+03_r8,&  
         0.2167562E+03_r8,     0.2153987E+03_r8,     0.2140477E+03_r8,     0.2127036E+03_r8,     0.2113662E+03_r8,&  
         0.2100359E+03_r8,     0.2087128E+03_r8,     0.2073970E+03_r8,     0.2060886E+03_r8,     0.2047876E+03_r8,&  
         0.2034941E+03_r8,     0.2022081E+03_r8,     0.2009298E+03_r8,     0.1996590E+03_r8,     0.1983959E+03_r8,&  
         0.1971405E+03_r8,     0.1958926E+03_r8,     0.1946524E+03_r8,     0.1934199E+03_r8,     0.1921950E+03_r8,&  
         0.1909776E+03_r8,     0.1897678E+03_r8,     0.1885656E+03_r8,     0.1873708E+03_r8,     0.1861836E+03_r8,&  
         0.1850038E+03_r8,     0.1838314E+03_r8,     0.1826664E+03_r8,     0.1815087E+03_r8,     0.1803583E+03_r8,&  
         0.1792152E+03_r8,     0.1780793E+03_r8,     0.1769506E+03_r8,     0.1758290E+03_r8,     0.1747145E+03_r8,&  
         0.1736070E+03_r8,     0.1725066E+03_r8,     0.1714131E+03_r8,     0.1703265E+03_r8,     0.1692469E+03_r8,&  
         0.1681741E+03_r8,     0.1671080E+03_r8,     0.1660488E+03_r8,     0.1649962E+03_r8,     0.1639503E+03_r8,&  
         0.1629110E+03_r8,     0.1618784E+03_r8,     0.1608522E+03_r8,     0.1598326E+03_r8,     0.1588194E+03_r8,&  
         0.1578126E+03_r8,     0.1568123E+03_r8,     0.1558183E+03_r8,     0.1548305E+03_r8,     0.1538491E+03_r8,&  
         0.1528738E+03_r8,     0.1519048E+03_r8,     0.1509418E+03_r8,     0.1499850E+03_r8,     0.1490343E+03_r8,&  
         0.1480895E+03_r8,     0.1471508E+03_r8,     0.1462180E+03_r8,     0.1452912E+03_r8,     0.1443702E+03_r8,&  
         0.1434550E+03_r8,     0.1425456E+03_r8,     0.1416420E+03_r8,     0.1407442E+03_r8,     0.1398520E+03_r8,&  
         0.1389655E+03_r8,     0.1380846E+03_r8,     0.1372093E+03_r8,     0.1363395E+03_r8,     0.1354753E+03_r8,&  
         0.1346165E+03_r8,     0.1337632E+03_r8,     0.1329152E+03_r8,     0.1320727E+03_r8,     0.1312355E+03_r8,&  
         0.1304036E+03_r8,     0.1295770E+03_r8,     0.1287556E+03_r8,     0.1279394E+03_r8,     0.1271284E+03_r8,&  
         0.1263225E+03_r8,     0.1255218E+03_r8,     0.1247261E+03_r8,     0.1239354E+03_r8,     0.1231498E+03_r8,&  
         0.1223692E+03_r8,     0.1215935E+03_r8,     0.1208227E+03_r8,     0.1200568E+03_r8,     0.1192958E+03_r8,&  
         0.1185396E+03_r8,     0.1177881E+03_r8,     0.1170415E+03_r8,     0.1162996E+03_r8,     0.1155623E+03_r8,&  
         0.1148298E+03_r8,     0.1141019E+03_r8,     0.1133786E+03_r8,     0.1126599E+03_r8,     0.1119457E+03_r8/)  

    psaditmk(1:150,   43)= (/ &
         0.2778113E+03_r8,     0.2767004E+03_r8,     0.2755780E+03_r8,     0.2744440E+03_r8,     0.2732981E+03_r8,&  
         0.2721401E+03_r8,     0.2709697E+03_r8,     0.2697870E+03_r8,     0.2685917E+03_r8,     0.2673839E+03_r8,&  
         0.2661632E+03_r8,     0.2649300E+03_r8,     0.2636841E+03_r8,     0.2624258E+03_r8,     0.2611551E+03_r8,&  
         0.2598722E+03_r8,     0.2585775E+03_r8,     0.2572713E+03_r8,     0.2559538E+03_r8,     0.2546257E+03_r8,&  
         0.2532874E+03_r8,     0.2519393E+03_r8,     0.2505823E+03_r8,     0.2492168E+03_r8,     0.2478436E+03_r8,&  
         0.2464635E+03_r8,     0.2450770E+03_r8,     0.2436851E+03_r8,     0.2422885E+03_r8,     0.2408880E+03_r8,&  
         0.2394844E+03_r8,     0.2380784E+03_r8,     0.2366709E+03_r8,     0.2352626E+03_r8,     0.2338542E+03_r8,&  
         0.2324465E+03_r8,     0.2310400E+03_r8,     0.2296354E+03_r8,     0.2282334E+03_r8,     0.2268344E+03_r8,&  
         0.2254390E+03_r8,     0.2240477E+03_r8,     0.2226609E+03_r8,     0.2212789E+03_r8,     0.2199023E+03_r8,&  
         0.2185313E+03_r8,     0.2171662E+03_r8,     0.2158072E+03_r8,     0.2144548E+03_r8,     0.2131089E+03_r8,&  
         0.2117697E+03_r8,     0.2104376E+03_r8,     0.2091126E+03_r8,     0.2077947E+03_r8,     0.2064843E+03_r8,&  
         0.2051812E+03_r8,     0.2038855E+03_r8,     0.2025974E+03_r8,     0.2013169E+03_r8,     0.2000439E+03_r8,&  
         0.1987785E+03_r8,     0.1975208E+03_r8,     0.1962707E+03_r8,     0.1950283E+03_r8,     0.1937935E+03_r8,&  
         0.1925662E+03_r8,     0.1913466E+03_r8,     0.1901346E+03_r8,     0.1889300E+03_r8,     0.1877330E+03_r8,&  
         0.1865435E+03_r8,     0.1853615E+03_r8,     0.1841869E+03_r8,     0.1830196E+03_r8,     0.1818597E+03_r8,&  
         0.1807071E+03_r8,     0.1795618E+03_r8,     0.1784237E+03_r8,     0.1772928E+03_r8,     0.1761691E+03_r8,&  
         0.1750524E+03_r8,     0.1739429E+03_r8,     0.1728403E+03_r8,     0.1717447E+03_r8,     0.1706561E+03_r8,&  
         0.1695743E+03_r8,     0.1684994E+03_r8,     0.1674313E+03_r8,     0.1663700E+03_r8,     0.1653154E+03_r8,&  
         0.1642675E+03_r8,     0.1632262E+03_r8,     0.1621915E+03_r8,     0.1611634E+03_r8,     0.1601418E+03_r8,&  
         0.1591266E+03_r8,     0.1581180E+03_r8,     0.1571156E+03_r8,     0.1561197E+03_r8,     0.1551301E+03_r8,&  
         0.1541467E+03_r8,     0.1531696E+03_r8,     0.1521986E+03_r8,     0.1512339E+03_r8,     0.1502752E+03_r8,&  
         0.1493226E+03_r8,     0.1483760E+03_r8,     0.1474355E+03_r8,     0.1465009E+03_r8,     0.1455722E+03_r8,&  
         0.1446494E+03_r8,     0.1437325E+03_r8,     0.1428214E+03_r8,     0.1419160E+03_r8,     0.1410164E+03_r8,&  
         0.1401226E+03_r8,     0.1392343E+03_r8,     0.1383517E+03_r8,     0.1374747E+03_r8,     0.1366033E+03_r8,&  
         0.1357373E+03_r8,     0.1348769E+03_r8,     0.1340219E+03_r8,     0.1331723E+03_r8,     0.1323282E+03_r8,&  
         0.1314893E+03_r8,     0.1306558E+03_r8,     0.1298276E+03_r8,     0.1290046E+03_r8,     0.1281869E+03_r8,&  
         0.1273743E+03_r8,     0.1265669E+03_r8,     0.1257646E+03_r8,     0.1249674E+03_r8,     0.1241752E+03_r8,&  
         0.1233880E+03_r8,     0.1226059E+03_r8,     0.1218287E+03_r8,     0.1210564E+03_r8,     0.1202891E+03_r8,&  
         0.1195266E+03_r8,     0.1187689E+03_r8,     0.1180160E+03_r8,     0.1172679E+03_r8,     0.1165245E+03_r8,&  
         0.1157859E+03_r8,     0.1150519E+03_r8,     0.1143226E+03_r8,     0.1135979E+03_r8,     0.1128778E+03_r8/)  

    psaditmk(1:150,   44)= (/ &
         0.2790297E+03_r8,     0.2779402E+03_r8,     0.2768397E+03_r8,     0.2757278E+03_r8,     0.2746043E+03_r8,&  
         0.2734689E+03_r8,     0.2723217E+03_r8,     0.2711623E+03_r8,     0.2699902E+03_r8,     0.2688058E+03_r8,&  
         0.2676086E+03_r8,     0.2663989E+03_r8,     0.2651764E+03_r8,     0.2639411E+03_r8,     0.2626932E+03_r8,&  
         0.2614329E+03_r8,     0.2601602E+03_r8,     0.2588754E+03_r8,     0.2575788E+03_r8,     0.2562707E+03_r8,&  
         0.2549516E+03_r8,     0.2536219E+03_r8,     0.2522822E+03_r8,     0.2509330E+03_r8,     0.2495750E+03_r8,&  
         0.2482088E+03_r8,     0.2468351E+03_r8,     0.2454548E+03_r8,     0.2440684E+03_r8,     0.2426768E+03_r8,&  
         0.2412809E+03_r8,     0.2398814E+03_r8,     0.2384790E+03_r8,     0.2370746E+03_r8,     0.2356691E+03_r8,&  
         0.2342629E+03_r8,     0.2328570E+03_r8,     0.2314520E+03_r8,     0.2300485E+03_r8,     0.2286472E+03_r8,&  
         0.2272486E+03_r8,     0.2258533E+03_r8,     0.2244619E+03_r8,     0.2230746E+03_r8,     0.2216921E+03_r8,&  
         0.2203146E+03_r8,     0.2189426E+03_r8,     0.2175762E+03_r8,     0.2162160E+03_r8,     0.2148620E+03_r8,&  
         0.2135144E+03_r8,     0.2121736E+03_r8,     0.2108397E+03_r8,     0.2095128E+03_r8,     0.2081930E+03_r8,&  
         0.2068805E+03_r8,     0.2055753E+03_r8,     0.2042775E+03_r8,     0.2029873E+03_r8,     0.2017045E+03_r8,&  
         0.2004294E+03_r8,     0.1991618E+03_r8,     0.1979018E+03_r8,     0.1966495E+03_r8,     0.1954048E+03_r8,&  
         0.1941677E+03_r8,     0.1929382E+03_r8,     0.1917163E+03_r8,     0.1905019E+03_r8,     0.1892952E+03_r8,&  
         0.1880959E+03_r8,     0.1869041E+03_r8,     0.1857198E+03_r8,     0.1845430E+03_r8,     0.1833735E+03_r8,&  
         0.1822114E+03_r8,     0.1810566E+03_r8,     0.1799091E+03_r8,     0.1787688E+03_r8,     0.1776357E+03_r8,&  
         0.1765098E+03_r8,     0.1753910E+03_r8,     0.1742793E+03_r8,     0.1731746E+03_r8,     0.1720769E+03_r8,&  
         0.1709861E+03_r8,     0.1699023E+03_r8,     0.1688253E+03_r8,     0.1677552E+03_r8,     0.1666918E+03_r8,&  
         0.1656351E+03_r8,     0.1645852E+03_r8,     0.1635419E+03_r8,     0.1625052E+03_r8,     0.1614751E+03_r8,&  
         0.1604515E+03_r8,     0.1594345E+03_r8,     0.1584238E+03_r8,     0.1574195E+03_r8,     0.1564217E+03_r8,&  
         0.1554301E+03_r8,     0.1544449E+03_r8,     0.1534658E+03_r8,     0.1524931E+03_r8,     0.1515264E+03_r8,&  
         0.1505659E+03_r8,     0.1496114E+03_r8,     0.1486630E+03_r8,     0.1477207E+03_r8,     0.1467843E+03_r8,&  
         0.1458538E+03_r8,     0.1449292E+03_r8,     0.1440105E+03_r8,     0.1430977E+03_r8,     0.1421906E+03_r8,&  
         0.1412892E+03_r8,     0.1403936E+03_r8,     0.1395036E+03_r8,     0.1386193E+03_r8,     0.1377406E+03_r8,&  
         0.1368675E+03_r8,     0.1359999E+03_r8,     0.1351378E+03_r8,     0.1342812E+03_r8,     0.1334300E+03_r8,&  
         0.1325842E+03_r8,     0.1317437E+03_r8,     0.1309086E+03_r8,     0.1300788E+03_r8,     0.1292542E+03_r8,&  
         0.1284349E+03_r8,     0.1276207E+03_r8,     0.1268117E+03_r8,     0.1260079E+03_r8,     0.1252091E+03_r8,&  
         0.1244154E+03_r8,     0.1236267E+03_r8,     0.1228431E+03_r8,     0.1220644E+03_r8,     0.1212906E+03_r8,&  
         0.1205218E+03_r8,     0.1197578E+03_r8,     0.1189986E+03_r8,     0.1182443E+03_r8,     0.1174947E+03_r8,&  
         0.1167500E+03_r8,     0.1160099E+03_r8,     0.1152745E+03_r8,     0.1145438E+03_r8,     0.1138177E+03_r8/)  

    psaditmk(1:150,   45)= (/ &
         0.2802106E+03_r8,     0.2791419E+03_r8,     0.2780625E+03_r8,     0.2769722E+03_r8,     0.2758707E+03_r8,&  
         0.2747577E+03_r8,     0.2736330E+03_r8,     0.2724964E+03_r8,     0.2713476E+03_r8,     0.2701866E+03_r8,&  
         0.2690128E+03_r8,     0.2678265E+03_r8,     0.2666275E+03_r8,     0.2654157E+03_r8,     0.2641913E+03_r8,&  
         0.2629539E+03_r8,     0.2617040E+03_r8,     0.2604416E+03_r8,     0.2591667E+03_r8,     0.2578800E+03_r8,&  
         0.2565815E+03_r8,     0.2552716E+03_r8,     0.2539510E+03_r8,     0.2526197E+03_r8,     0.2512786E+03_r8,&  
         0.2499284E+03_r8,     0.2485694E+03_r8,     0.2472025E+03_r8,     0.2458285E+03_r8,     0.2444479E+03_r8,&  
         0.2430617E+03_r8,     0.2416707E+03_r8,     0.2402755E+03_r8,     0.2388770E+03_r8,     0.2374761E+03_r8,&  
         0.2360734E+03_r8,     0.2346697E+03_r8,     0.2332659E+03_r8,     0.2318626E+03_r8,     0.2304604E+03_r8,&  
         0.2290600E+03_r8,     0.2276620E+03_r8,     0.2262671E+03_r8,     0.2248756E+03_r8,     0.2234881E+03_r8,&  
         0.2221050E+03_r8,     0.2207268E+03_r8,     0.2193538E+03_r8,     0.2179864E+03_r8,     0.2166248E+03_r8,&  
         0.2152694E+03_r8,     0.2139203E+03_r8,     0.2125779E+03_r8,     0.2112422E+03_r8,     0.2099134E+03_r8,&  
         0.2085917E+03_r8,     0.2072772E+03_r8,     0.2059700E+03_r8,     0.2046701E+03_r8,     0.2033777E+03_r8,&  
         0.2020928E+03_r8,     0.2008155E+03_r8,     0.1995457E+03_r8,     0.1982835E+03_r8,     0.1970289E+03_r8,&  
         0.1957819E+03_r8,     0.1945425E+03_r8,     0.1933108E+03_r8,     0.1920866E+03_r8,     0.1908700E+03_r8,&  
         0.1896610E+03_r8,     0.1884594E+03_r8,     0.1872655E+03_r8,     0.1860789E+03_r8,     0.1848998E+03_r8,&  
         0.1837281E+03_r8,     0.1825637E+03_r8,     0.1814067E+03_r8,     0.1802570E+03_r8,     0.1791145E+03_r8,&  
         0.1779793E+03_r8,     0.1768512E+03_r8,     0.1757302E+03_r8,     0.1746164E+03_r8,     0.1735096E+03_r8,&  
         0.1724097E+03_r8,     0.1713169E+03_r8,     0.1702310E+03_r8,     0.1691519E+03_r8,     0.1680797E+03_r8,&  
         0.1670143E+03_r8,     0.1659556E+03_r8,     0.1649036E+03_r8,     0.1638583E+03_r8,     0.1628196E+03_r8,&  
         0.1617875E+03_r8,     0.1607619E+03_r8,     0.1597429E+03_r8,     0.1587302E+03_r8,     0.1577241E+03_r8,&  
         0.1567243E+03_r8,     0.1557308E+03_r8,     0.1547436E+03_r8,     0.1537627E+03_r8,     0.1527880E+03_r8,&  
         0.1518195E+03_r8,     0.1508571E+03_r8,     0.1499008E+03_r8,     0.1489506E+03_r8,     0.1480064E+03_r8,&  
         0.1470682E+03_r8,     0.1461360E+03_r8,     0.1452096E+03_r8,     0.1442891E+03_r8,     0.1433745E+03_r8,&  
         0.1424656E+03_r8,     0.1415625E+03_r8,     0.1406652E+03_r8,     0.1397735E+03_r8,     0.1388875E+03_r8,&  
         0.1380071E+03_r8,     0.1371323E+03_r8,     0.1362630E+03_r8,     0.1353992E+03_r8,     0.1345409E+03_r8,&  
         0.1336881E+03_r8,     0.1328406E+03_r8,     0.1319986E+03_r8,     0.1311618E+03_r8,     0.1303304E+03_r8,&  
         0.1295042E+03_r8,     0.1286833E+03_r8,     0.1278676E+03_r8,     0.1270570E+03_r8,     0.1262516E+03_r8,&  
         0.1254513E+03_r8,     0.1246561E+03_r8,     0.1238659E+03_r8,     0.1230807E+03_r8,     0.1223005E+03_r8,&  
         0.1215252E+03_r8,     0.1207549E+03_r8,     0.1199894E+03_r8,     0.1192288E+03_r8,     0.1184730E+03_r8,&  
         0.1177220E+03_r8,     0.1169758E+03_r8,     0.1162343E+03_r8,     0.1154975E+03_r8,     0.1147653E+03_r8/)  

    psaditmk(1:150,   46)= (/ &
         0.2813550E+03_r8,     0.2803063E+03_r8,     0.2792475E+03_r8,     0.2781781E+03_r8,     0.2770980E+03_r8,&  
         0.2760068E+03_r8,     0.2749042E+03_r8,     0.2737901E+03_r8,     0.2726641E+03_r8,     0.2715260E+03_r8,&  
         0.2703756E+03_r8,     0.2692127E+03_r8,     0.2680374E+03_r8,     0.2668492E+03_r8,     0.2656482E+03_r8,&  
         0.2644343E+03_r8,     0.2632078E+03_r8,     0.2619685E+03_r8,     0.2607164E+03_r8,     0.2594519E+03_r8,&  
         0.2581750E+03_r8,     0.2568862E+03_r8,     0.2555858E+03_r8,     0.2542742E+03_r8,     0.2529517E+03_r8,&  
         0.2516191E+03_r8,     0.2502767E+03_r8,     0.2489253E+03_r8,     0.2475656E+03_r8,     0.2461980E+03_r8,&  
         0.2448236E+03_r8,     0.2434431E+03_r8,     0.2420571E+03_r8,     0.2406666E+03_r8,     0.2392723E+03_r8,&  
         0.2378750E+03_r8,     0.2364756E+03_r8,     0.2350747E+03_r8,     0.2336732E+03_r8,     0.2322717E+03_r8,&  
         0.2308710E+03_r8,     0.2294717E+03_r8,     0.2280746E+03_r8,     0.2266800E+03_r8,     0.2252887E+03_r8,&  
         0.2239011E+03_r8,     0.2225177E+03_r8,     0.2211388E+03_r8,     0.2197650E+03_r8,     0.2183966E+03_r8,&  
         0.2170338E+03_r8,     0.2156771E+03_r8,     0.2143265E+03_r8,     0.2129825E+03_r8,     0.2116451E+03_r8,&  
         0.2103145E+03_r8,     0.2089909E+03_r8,     0.2076744E+03_r8,     0.2063652E+03_r8,     0.2050633E+03_r8,&  
         0.2037688E+03_r8,     0.2024817E+03_r8,     0.2012022E+03_r8,     0.1999302E+03_r8,     0.1986658E+03_r8,&  
         0.1974090E+03_r8,     0.1961597E+03_r8,     0.1949181E+03_r8,     0.1936841E+03_r8,     0.1924576E+03_r8,&  
         0.1912388E+03_r8,     0.1900275E+03_r8,     0.1888237E+03_r8,     0.1876274E+03_r8,     0.1864386E+03_r8,&  
         0.1852573E+03_r8,     0.1840833E+03_r8,     0.1829167E+03_r8,     0.1817575E+03_r8,     0.1806056E+03_r8,&  
         0.1794610E+03_r8,     0.1783235E+03_r8,     0.1771933E+03_r8,     0.1760701E+03_r8,     0.1749541E+03_r8,&  
         0.1738452E+03_r8,     0.1727432E+03_r8,     0.1716483E+03_r8,     0.1705602E+03_r8,     0.1694791E+03_r8,&  
         0.1684048E+03_r8,     0.1673373E+03_r8,     0.1662766E+03_r8,     0.1652226E+03_r8,     0.1641753E+03_r8,&  
         0.1631345E+03_r8,     0.1621005E+03_r8,     0.1610729E+03_r8,     0.1600519E+03_r8,     0.1590373E+03_r8,&  
         0.1580292E+03_r8,     0.1570274E+03_r8,     0.1560320E+03_r8,     0.1550430E+03_r8,     0.1540602E+03_r8,&  
         0.1530836E+03_r8,     0.1521132E+03_r8,     0.1511489E+03_r8,     0.1501908E+03_r8,     0.1492388E+03_r8,&  
         0.1482927E+03_r8,     0.1473527E+03_r8,     0.1464187E+03_r8,     0.1454905E+03_r8,     0.1445682E+03_r8,&  
         0.1436518E+03_r8,     0.1427412E+03_r8,     0.1418364E+03_r8,     0.1409373E+03_r8,     0.1400439E+03_r8,&  
         0.1391561E+03_r8,     0.1382740E+03_r8,     0.1373975E+03_r8,     0.1365266E+03_r8,     0.1356611E+03_r8,&  
         0.1348012E+03_r8,     0.1339467E+03_r8,     0.1330976E+03_r8,     0.1322539E+03_r8,     0.1314155E+03_r8,&  
         0.1305825E+03_r8,     0.1297547E+03_r8,     0.1289322E+03_r8,     0.1281149E+03_r8,     0.1273028E+03_r8,&  
         0.1264958E+03_r8,     0.1256940E+03_r8,     0.1248972E+03_r8,     0.1241055E+03_r8,     0.1233188E+03_r8,&  
         0.1225371E+03_r8,     0.1217603E+03_r8,     0.1209885E+03_r8,     0.1202215E+03_r8,     0.1194595E+03_r8,&  
         0.1187022E+03_r8,     0.1179498E+03_r8,     0.1172021E+03_r8,     0.1164591E+03_r8,     0.1157209E+03_r8/)  

    psaditmk(1:150,   47)= (/ &
         0.2824640E+03_r8,     0.2814348E+03_r8,     0.2803958E+03_r8,     0.2793467E+03_r8,     0.2782873E+03_r8,&  
         0.2772172E+03_r8,     0.2761363E+03_r8,     0.2750441E+03_r8,     0.2739404E+03_r8,     0.2728250E+03_r8,&  
         0.2716977E+03_r8,     0.2705579E+03_r8,     0.2694058E+03_r8,     0.2682412E+03_r8,     0.2670639E+03_r8,&  
         0.2658737E+03_r8,     0.2646708E+03_r8,     0.2634548E+03_r8,     0.2622260E+03_r8,     0.2609844E+03_r8,&  
         0.2597303E+03_r8,     0.2584636E+03_r8,     0.2571847E+03_r8,     0.2558939E+03_r8,     0.2545915E+03_r8,&  
         0.2532781E+03_r8,     0.2519540E+03_r8,     0.2506199E+03_r8,     0.2492763E+03_r8,     0.2479238E+03_r8,&  
         0.2465633E+03_r8,     0.2451953E+03_r8,     0.2438206E+03_r8,     0.2424401E+03_r8,     0.2410546E+03_r8,&  
         0.2396647E+03_r8,     0.2382714E+03_r8,     0.2368754E+03_r8,     0.2354775E+03_r8,     0.2340785E+03_r8,&  
         0.2326792E+03_r8,     0.2312803E+03_r8,     0.2298824E+03_r8,     0.2284861E+03_r8,     0.2270922E+03_r8,&  
         0.2257012E+03_r8,     0.2243136E+03_r8,     0.2229299E+03_r8,     0.2215506E+03_r8,     0.2201761E+03_r8,&  
         0.2188067E+03_r8,     0.2174429E+03_r8,     0.2160849E+03_r8,     0.2147329E+03_r8,     0.2133874E+03_r8,&  
         0.2120483E+03_r8,     0.2107159E+03_r8,     0.2093905E+03_r8,     0.2080722E+03_r8,     0.2067609E+03_r8,&  
         0.2054570E+03_r8,     0.2041604E+03_r8,     0.2028712E+03_r8,     0.2015896E+03_r8,     0.2003154E+03_r8,&  
         0.1990487E+03_r8,     0.1977897E+03_r8,     0.1965382E+03_r8,     0.1952943E+03_r8,     0.1940580E+03_r8,&  
         0.1928293E+03_r8,     0.1916082E+03_r8,     0.1903946E+03_r8,     0.1891886E+03_r8,     0.1879901E+03_r8,&  
         0.1867990E+03_r8,     0.1856154E+03_r8,     0.1844392E+03_r8,     0.1832705E+03_r8,     0.1821090E+03_r8,&  
         0.1809549E+03_r8,     0.1798080E+03_r8,     0.1786684E+03_r8,     0.1775360E+03_r8,     0.1764107E+03_r8,&  
         0.1752925E+03_r8,     0.1741814E+03_r8,     0.1730774E+03_r8,     0.1719803E+03_r8,     0.1708902E+03_r8,&  
         0.1698069E+03_r8,     0.1687306E+03_r8,     0.1676610E+03_r8,     0.1665982E+03_r8,     0.1655422E+03_r8,&  
         0.1644928E+03_r8,     0.1634501E+03_r8,     0.1624140E+03_r8,     0.1613845E+03_r8,     0.1603615E+03_r8,&  
         0.1593450E+03_r8,     0.1583349E+03_r8,     0.1573312E+03_r8,     0.1563339E+03_r8,     0.1553429E+03_r8,&  
         0.1543582E+03_r8,     0.1533797E+03_r8,     0.1524074E+03_r8,     0.1514413E+03_r8,     0.1504814E+03_r8,&  
         0.1495275E+03_r8,     0.1485796E+03_r8,     0.1476378E+03_r8,     0.1467019E+03_r8,     0.1457719E+03_r8,&  
         0.1448479E+03_r8,     0.1439297E+03_r8,     0.1430173E+03_r8,     0.1421108E+03_r8,     0.1412099E+03_r8,&  
         0.1403148E+03_r8,     0.1394254E+03_r8,     0.1385415E+03_r8,     0.1376633E+03_r8,     0.1367907E+03_r8,&  
         0.1359236E+03_r8,     0.1350620E+03_r8,     0.1342058E+03_r8,     0.1333551E+03_r8,     0.1325097E+03_r8,&  
         0.1316698E+03_r8,     0.1308351E+03_r8,     0.1300057E+03_r8,     0.1291817E+03_r8,     0.1283628E+03_r8,&  
         0.1275491E+03_r8,     0.1267405E+03_r8,     0.1259371E+03_r8,     0.1251388E+03_r8,     0.1243456E+03_r8,&  
         0.1235574E+03_r8,     0.1227741E+03_r8,     0.1219959E+03_r8,     0.1212225E+03_r8,     0.1204541E+03_r8,&  
         0.1196906E+03_r8,     0.1189318E+03_r8,     0.1181779E+03_r8,     0.1174288E+03_r8,     0.1166844E+03_r8/)  

    psaditmk(1:150,   48)= (/ &
         0.2835388E+03_r8,     0.2825283E+03_r8,     0.2815084E+03_r8,     0.2804789E+03_r8,     0.2794395E+03_r8,&  
         0.2783900E+03_r8,     0.2773300E+03_r8,     0.2762591E+03_r8,     0.2751772E+03_r8,     0.2740839E+03_r8,&  
         0.2729789E+03_r8,     0.2718621E+03_r8,     0.2707331E+03_r8,     0.2695919E+03_r8,     0.2684380E+03_r8,&  
         0.2672717E+03_r8,     0.2660922E+03_r8,     0.2649001E+03_r8,     0.2636949E+03_r8,     0.2624767E+03_r8,&  
         0.2612458E+03_r8,     0.2600020E+03_r8,     0.2587455E+03_r8,     0.2574767E+03_r8,     0.2561956E+03_r8,&  
         0.2549029E+03_r8,     0.2535987E+03_r8,     0.2522834E+03_r8,     0.2509577E+03_r8,     0.2496222E+03_r8,&  
         0.2482774E+03_r8,     0.2469240E+03_r8,     0.2455628E+03_r8,     0.2441943E+03_r8,     0.2428195E+03_r8,&  
         0.2414392E+03_r8,     0.2400540E+03_r8,     0.2386649E+03_r8,     0.2372727E+03_r8,     0.2358781E+03_r8,&  
         0.2344819E+03_r8,     0.2330850E+03_r8,     0.2316880E+03_r8,     0.2302916E+03_r8,     0.2288965E+03_r8,&  
         0.2275035E+03_r8,     0.2261130E+03_r8,     0.2247255E+03_r8,     0.2233418E+03_r8,     0.2219621E+03_r8,&  
         0.2205870E+03_r8,     0.2192168E+03_r8,     0.2178520E+03_r8,     0.2164928E+03_r8,     0.2151396E+03_r8,&  
         0.2137925E+03_r8,     0.2124518E+03_r8,     0.2111178E+03_r8,     0.2097906E+03_r8,     0.2084703E+03_r8,&  
         0.2071572E+03_r8,     0.2058512E+03_r8,     0.2045526E+03_r8,     0.2032613E+03_r8,     0.2019775E+03_r8,&  
         0.2007011E+03_r8,     0.1994323E+03_r8,     0.1981710E+03_r8,     0.1969173E+03_r8,     0.1956712E+03_r8,&  
         0.1944327E+03_r8,     0.1932017E+03_r8,     0.1919783E+03_r8,     0.1907625E+03_r8,     0.1895542E+03_r8,&  
         0.1883534E+03_r8,     0.1871601E+03_r8,     0.1859742E+03_r8,     0.1847958E+03_r8,     0.1836248E+03_r8,&  
         0.1824612E+03_r8,     0.1813048E+03_r8,     0.1801557E+03_r8,     0.1790139E+03_r8,     0.1778793E+03_r8,&  
         0.1767519E+03_r8,     0.1756316E+03_r8,     0.1745184E+03_r8,     0.1734122E+03_r8,     0.1723129E+03_r8,&  
         0.1712207E+03_r8,     0.1701354E+03_r8,     0.1690569E+03_r8,     0.1679853E+03_r8,     0.1669205E+03_r8,&  
         0.1658624E+03_r8,     0.1648110E+03_r8,     0.1637663E+03_r8,     0.1627282E+03_r8,     0.1616967E+03_r8,&  
         0.1606717E+03_r8,     0.1596532E+03_r8,     0.1586412E+03_r8,     0.1576356E+03_r8,     0.1566363E+03_r8,&  
         0.1556434E+03_r8,     0.1546568E+03_r8,     0.1536764E+03_r8,     0.1527023E+03_r8,     0.1517343E+03_r8,&  
         0.1507724E+03_r8,     0.1498167E+03_r8,     0.1488670E+03_r8,     0.1479234E+03_r8,     0.1469857E+03_r8,&  
         0.1460539E+03_r8,     0.1451281E+03_r8,     0.1442081E+03_r8,     0.1432940E+03_r8,     0.1423857E+03_r8,&  
         0.1414831E+03_r8,     0.1405862E+03_r8,     0.1396951E+03_r8,     0.1388095E+03_r8,     0.1379296E+03_r8,&  
         0.1370553E+03_r8,     0.1361865E+03_r8,     0.1353232E+03_r8,     0.1344654E+03_r8,     0.1336130E+03_r8,&  
         0.1327661E+03_r8,     0.1319245E+03_r8,     0.1310882E+03_r8,     0.1302572E+03_r8,     0.1294315E+03_r8,&  
         0.1286111E+03_r8,     0.1277958E+03_r8,     0.1269857E+03_r8,     0.1261808E+03_r8,     0.1253809E+03_r8,&  
         0.1245861E+03_r8,     0.1237964E+03_r8,     0.1230116E+03_r8,     0.1222318E+03_r8,     0.1214570E+03_r8,&  
         0.1206871E+03_r8,     0.1199221E+03_r8,     0.1191619E+03_r8,     0.1184065E+03_r8,     0.1176560E+03_r8/)  

    psaditmk(1:150,   49)= (/ &
         0.2845808E+03_r8,     0.2835882E+03_r8,     0.2825866E+03_r8,     0.2815757E+03_r8,     0.2805559E+03_r8,&  
         0.2795262E+03_r8,     0.2784862E+03_r8,     0.2774361E+03_r8,     0.2763753E+03_r8,     0.2753036E+03_r8,&  
         0.2742206E+03_r8,     0.2731261E+03_r8,     0.2720198E+03_r8,     0.2709015E+03_r8,     0.2697710E+03_r8,&  
         0.2686279E+03_r8,     0.2674722E+03_r8,     0.2663037E+03_r8,     0.2651224E+03_r8,     0.2639280E+03_r8,&  
         0.2627207E+03_r8,     0.2615004E+03_r8,     0.2602671E+03_r8,     0.2590210E+03_r8,     0.2577623E+03_r8,&  
         0.2564913E+03_r8,     0.2552082E+03_r8,     0.2539133E+03_r8,     0.2526071E+03_r8,     0.2512902E+03_r8,&  
         0.2499629E+03_r8,     0.2486261E+03_r8,     0.2472800E+03_r8,     0.2459257E+03_r8,     0.2445639E+03_r8,&  
         0.2431951E+03_r8,     0.2418202E+03_r8,     0.2404401E+03_r8,     0.2390555E+03_r8,     0.2376673E+03_r8,&  
         0.2362763E+03_r8,     0.2348831E+03_r8,     0.2334889E+03_r8,     0.2320940E+03_r8,     0.2306994E+03_r8,&  
         0.2293058E+03_r8,     0.2279137E+03_r8,     0.2265239E+03_r8,     0.2251368E+03_r8,     0.2237531E+03_r8,&  
         0.2223732E+03_r8,     0.2209977E+03_r8,     0.2196268E+03_r8,     0.2182611E+03_r8,     0.2169008E+03_r8,&  
         0.2155463E+03_r8,     0.2141978E+03_r8,     0.2128557E+03_r8,     0.2115200E+03_r8,     0.2101911E+03_r8,&  
         0.2088690E+03_r8,     0.2075539E+03_r8,     0.2062460E+03_r8,     0.2049453E+03_r8,     0.2036519E+03_r8,&  
         0.2023660E+03_r8,     0.2010875E+03_r8,     0.1998165E+03_r8,     0.1985530E+03_r8,     0.1972971E+03_r8,&  
         0.1960487E+03_r8,     0.1948079E+03_r8,     0.1935747E+03_r8,     0.1923491E+03_r8,     0.1911310E+03_r8,&  
         0.1899204E+03_r8,     0.1887174E+03_r8,     0.1875218E+03_r8,     0.1863337E+03_r8,     0.1851531E+03_r8,&  
         0.1839798E+03_r8,     0.1828140E+03_r8,     0.1816554E+03_r8,     0.1805041E+03_r8,     0.1793601E+03_r8,&  
         0.1782234E+03_r8,     0.1770937E+03_r8,     0.1759713E+03_r8,     0.1748559E+03_r8,     0.1737476E+03_r8,&  
         0.1726463E+03_r8,     0.1715519E+03_r8,     0.1704645E+03_r8,     0.1693839E+03_r8,     0.1683103E+03_r8,&  
         0.1672434E+03_r8,     0.1661833E+03_r8,     0.1651298E+03_r8,     0.1640831E+03_r8,     0.1630430E+03_r8,&  
         0.1620095E+03_r8,     0.1609825E+03_r8,     0.1599620E+03_r8,     0.1589481E+03_r8,     0.1579405E+03_r8,&  
         0.1569393E+03_r8,     0.1559445E+03_r8,     0.1549559E+03_r8,     0.1539737E+03_r8,     0.1529977E+03_r8,&  
         0.1520278E+03_r8,     0.1510641E+03_r8,     0.1501065E+03_r8,     0.1491550E+03_r8,     0.1482095E+03_r8,&  
         0.1472700E+03_r8,     0.1463365E+03_r8,     0.1454088E+03_r8,     0.1444871E+03_r8,     0.1435712E+03_r8,&  
         0.1426611E+03_r8,     0.1417568E+03_r8,     0.1408582E+03_r8,     0.1399653E+03_r8,     0.1390780E+03_r8,&  
         0.1381964E+03_r8,     0.1373204E+03_r8,     0.1364499E+03_r8,     0.1355850E+03_r8,     0.1347255E+03_r8,&  
         0.1338715E+03_r8,     0.1330229E+03_r8,     0.1321796E+03_r8,     0.1313418E+03_r8,     0.1305092E+03_r8,&  
         0.1296819E+03_r8,     0.1288599E+03_r8,     0.1280430E+03_r8,     0.1272314E+03_r8,     0.1264248E+03_r8,&  
         0.1256234E+03_r8,     0.1248271E+03_r8,     0.1240358E+03_r8,     0.1232496E+03_r8,     0.1224683E+03_r8,&  
         0.1216920E+03_r8,     0.1209206E+03_r8,     0.1201541E+03_r8,     0.1193924E+03_r8,     0.1186356E+03_r8/)  

    psaditmk(1:150,   50)= (/ &
         0.2855910E+03_r8,     0.2846154E+03_r8,     0.2836314E+03_r8,     0.2826389E+03_r8,     0.2816375E+03_r8,&  
         0.2806268E+03_r8,     0.2796065E+03_r8,     0.2785764E+03_r8,     0.2775360E+03_r8,     0.2764851E+03_r8,&  
         0.2754234E+03_r8,     0.2743507E+03_r8,     0.2732665E+03_r8,     0.2721705E+03_r8,     0.2710629E+03_r8,&  
         0.2699431E+03_r8,     0.2688108E+03_r8,     0.2676659E+03_r8,     0.2665083E+03_r8,     0.2653378E+03_r8,&  
         0.2641542E+03_r8,     0.2629576E+03_r8,     0.2617482E+03_r8,     0.2605254E+03_r8,     0.2592898E+03_r8,&  
         0.2580415E+03_r8,     0.2567805E+03_r8,     0.2555072E+03_r8,     0.2542220E+03_r8,     0.2529251E+03_r8,&  
         0.2516170E+03_r8,     0.2502983E+03_r8,     0.2489696E+03_r8,     0.2476314E+03_r8,     0.2462843E+03_r8,&  
         0.2449291E+03_r8,     0.2435666E+03_r8,     0.2421976E+03_r8,     0.2408228E+03_r8,     0.2394430E+03_r8,&  
         0.2380591E+03_r8,     0.2366719E+03_r8,     0.2352822E+03_r8,     0.2338908E+03_r8,     0.2324984E+03_r8,&  
         0.2311058E+03_r8,     0.2297137E+03_r8,     0.2283229E+03_r8,     0.2269339E+03_r8,     0.2255473E+03_r8,&  
         0.2241638E+03_r8,     0.2227839E+03_r8,     0.2214080E+03_r8,     0.2200366E+03_r8,     0.2186701E+03_r8,&  
         0.2173088E+03_r8,     0.2159532E+03_r8,     0.2146034E+03_r8,     0.2132598E+03_r8,     0.2119225E+03_r8,&  
         0.2105919E+03_r8,     0.2092681E+03_r8,     0.2079511E+03_r8,     0.2066413E+03_r8,     0.2053386E+03_r8,&  
         0.2040432E+03_r8,     0.2027551E+03_r8,     0.2014745E+03_r8,     0.2002013E+03_r8,     0.1989356E+03_r8,&  
         0.1976775E+03_r8,     0.1964269E+03_r8,     0.1951839E+03_r8,     0.1939484E+03_r8,     0.1927205E+03_r8,&  
         0.1915002E+03_r8,     0.1902874E+03_r8,     0.1890821E+03_r8,     0.1878843E+03_r8,     0.1866939E+03_r8,&  
         0.1855110E+03_r8,     0.1843356E+03_r8,     0.1831674E+03_r8,     0.1820067E+03_r8,     0.1808532E+03_r8,&  
         0.1797070E+03_r8,     0.1785680E+03_r8,     0.1774363E+03_r8,     0.1763116E+03_r8,     0.1751941E+03_r8,&  
         0.1740836E+03_r8,     0.1729802E+03_r8,     0.1718837E+03_r8,     0.1707942E+03_r8,     0.1697116E+03_r8,&  
         0.1686358E+03_r8,     0.1675669E+03_r8,     0.1665047E+03_r8,     0.1654492E+03_r8,     0.1644005E+03_r8,&  
         0.1633584E+03_r8,     0.1623228E+03_r8,     0.1612939E+03_r8,     0.1602715E+03_r8,     0.1592555E+03_r8,&  
         0.1582460E+03_r8,     0.1572429E+03_r8,     0.1562461E+03_r8,     0.1552557E+03_r8,     0.1542715E+03_r8,&  
         0.1532936E+03_r8,     0.1523219E+03_r8,     0.1513563E+03_r8,     0.1503969E+03_r8,     0.1494435E+03_r8,&  
         0.1484962E+03_r8,     0.1475549E+03_r8,     0.1466195E+03_r8,     0.1456901E+03_r8,     0.1447666E+03_r8,&  
         0.1438489E+03_r8,     0.1429371E+03_r8,     0.1420310E+03_r8,     0.1411307E+03_r8,     0.1402361E+03_r8,&  
         0.1393471E+03_r8,     0.1384638E+03_r8,     0.1375861E+03_r8,     0.1367139E+03_r8,     0.1358473E+03_r8,&  
         0.1349861E+03_r8,     0.1341305E+03_r8,     0.1332802E+03_r8,     0.1324354E+03_r8,     0.1315959E+03_r8,&  
         0.1307617E+03_r8,     0.1299328E+03_r8,     0.1291091E+03_r8,     0.1282907E+03_r8,     0.1274775E+03_r8,&  
         0.1266694E+03_r8,     0.1258665E+03_r8,     0.1250686E+03_r8,     0.1242758E+03_r8,     0.1234880E+03_r8,&  
         0.1227052E+03_r8,     0.1219274E+03_r8,     0.1211545E+03_r8,     0.1203865E+03_r8,     0.1196234E+03_r8/)  

    psaditmk(1:150,   51)= (/ &
         0.2865706E+03_r8,     0.2856115E+03_r8,     0.2846445E+03_r8,     0.2836693E+03_r8,     0.2826856E+03_r8,&  
         0.2816933E+03_r8,     0.2806917E+03_r8,     0.2796808E+03_r8,     0.2786601E+03_r8,     0.2776295E+03_r8,&  
         0.2765884E+03_r8,     0.2755368E+03_r8,     0.2744742E+03_r8,     0.2734003E+03_r8,     0.2723148E+03_r8,&  
         0.2712178E+03_r8,     0.2701082E+03_r8,     0.2689868E+03_r8,     0.2678527E+03_r8,     0.2667058E+03_r8,&  
         0.2655461E+03_r8,     0.2643735E+03_r8,     0.2631878E+03_r8,     0.2619889E+03_r8,     0.2607769E+03_r8,&  
         0.2595520E+03_r8,     0.2583141E+03_r8,     0.2570633E+03_r8,     0.2558000E+03_r8,     0.2545246E+03_r8,&  
         0.2532371E+03_r8,     0.2519382E+03_r8,     0.2506283E+03_r8,     0.2493079E+03_r8,     0.2479776E+03_r8,&  
         0.2466381E+03_r8,     0.2452900E+03_r8,     0.2439341E+03_r8,     0.2425712E+03_r8,     0.2412020E+03_r8,&  
         0.2398273E+03_r8,     0.2384480E+03_r8,     0.2370649E+03_r8,     0.2356788E+03_r8,     0.2342905E+03_r8,&  
         0.2329008E+03_r8,     0.2315104E+03_r8,     0.2301202E+03_r8,     0.2287308E+03_r8,     0.2273429E+03_r8,&  
         0.2259570E+03_r8,     0.2245739E+03_r8,     0.2231941E+03_r8,     0.2218180E+03_r8,     0.2204461E+03_r8,&  
         0.2190789E+03_r8,     0.2177168E+03_r8,     0.2163601E+03_r8,     0.2150091E+03_r8,     0.2136641E+03_r8,&  
         0.2123254E+03_r8,     0.2109931E+03_r8,     0.2096675E+03_r8,     0.2083488E+03_r8,     0.2070370E+03_r8,&  
         0.2057323E+03_r8,     0.2044349E+03_r8,     0.2031448E+03_r8,     0.2018620E+03_r8,     0.2005867E+03_r8,&  
         0.1993189E+03_r8,     0.1980585E+03_r8,     0.1968057E+03_r8,     0.1955605E+03_r8,     0.1943228E+03_r8,&  
         0.1930926E+03_r8,     0.1918700E+03_r8,     0.1906549E+03_r8,     0.1894474E+03_r8,     0.1882474E+03_r8,&  
         0.1870548E+03_r8,     0.1858697E+03_r8,     0.1846919E+03_r8,     0.1835216E+03_r8,     0.1823586E+03_r8,&  
         0.1812029E+03_r8,     0.1800546E+03_r8,     0.1789134E+03_r8,     0.1777795E+03_r8,     0.1766526E+03_r8,&  
         0.1755330E+03_r8,     0.1744204E+03_r8,     0.1733148E+03_r8,     0.1722162E+03_r8,     0.1711246E+03_r8,&  
         0.1700399E+03_r8,     0.1689621E+03_r8,     0.1678910E+03_r8,     0.1668268E+03_r8,     0.1657693E+03_r8,&  
         0.1647185E+03_r8,     0.1636744E+03_r8,     0.1626369E+03_r8,     0.1616059E+03_r8,     0.1605815E+03_r8,&  
         0.1595636E+03_r8,     0.1585521E+03_r8,     0.1575471E+03_r8,     0.1565484E+03_r8,     0.1555560E+03_r8,&  
         0.1545700E+03_r8,     0.1535901E+03_r8,     0.1526166E+03_r8,     0.1516491E+03_r8,     0.1506878E+03_r8,&  
         0.1497326E+03_r8,     0.1487835E+03_r8,     0.1478403E+03_r8,     0.1469032E+03_r8,     0.1459720E+03_r8,&  
         0.1450466E+03_r8,     0.1441272E+03_r8,     0.1432136E+03_r8,     0.1423057E+03_r8,     0.1414037E+03_r8,&  
         0.1405073E+03_r8,     0.1396167E+03_r8,     0.1387316E+03_r8,     0.1378522E+03_r8,     0.1369784E+03_r8,&  
         0.1361101E+03_r8,     0.1352473E+03_r8,     0.1343899E+03_r8,     0.1335380E+03_r8,     0.1326915E+03_r8,&  
         0.1318504E+03_r8,     0.1310146E+03_r8,     0.1301841E+03_r8,     0.1293589E+03_r8,     0.1285389E+03_r8,&  
         0.1277241E+03_r8,     0.1269144E+03_r8,     0.1261099E+03_r8,     0.1253105E+03_r8,     0.1245162E+03_r8,&  
         0.1237269E+03_r8,     0.1229426E+03_r8,     0.1221632E+03_r8,     0.1213889E+03_r8,     0.1206194E+03_r8/)  

    psaditmk(1:150,   52)= (/ &
         0.2875207E+03_r8,     0.2865776E+03_r8,     0.2856266E+03_r8,     0.2846679E+03_r8,     0.2837015E+03_r8,&  
         0.2827266E+03_r8,     0.2817432E+03_r8,     0.2807507E+03_r8,     0.2797491E+03_r8,     0.2787379E+03_r8,&  
         0.2777167E+03_r8,     0.2766854E+03_r8,     0.2756436E+03_r8,     0.2745911E+03_r8,     0.2735274E+03_r8,&  
         0.2724522E+03_r8,     0.2713654E+03_r8,     0.2702670E+03_r8,     0.2691559E+03_r8,     0.2680325E+03_r8,&  
         0.2668964E+03_r8,     0.2657476E+03_r8,     0.2645857E+03_r8,     0.2634109E+03_r8,     0.2622229E+03_r8,&  
         0.2610216E+03_r8,     0.2598073E+03_r8,     0.2585798E+03_r8,     0.2573395E+03_r8,     0.2560864E+03_r8,&  
         0.2548208E+03_r8,     0.2535430E+03_r8,     0.2522535E+03_r8,     0.2509526E+03_r8,     0.2496408E+03_r8,&  
         0.2483188E+03_r8,     0.2469870E+03_r8,     0.2456464E+03_r8,     0.2442973E+03_r8,     0.2429408E+03_r8,&  
         0.2415775E+03_r8,     0.2402082E+03_r8,     0.2388338E+03_r8,     0.2374550E+03_r8,     0.2360728E+03_r8,&  
         0.2346879E+03_r8,     0.2333012E+03_r8,     0.2319133E+03_r8,     0.2305252E+03_r8,     0.2291374E+03_r8,&  
         0.2277507E+03_r8,     0.2263658E+03_r8,     0.2249832E+03_r8,     0.2236036E+03_r8,     0.2222275E+03_r8,&  
         0.2208553E+03_r8,     0.2194876E+03_r8,     0.2181247E+03_r8,     0.2167670E+03_r8,     0.2154149E+03_r8,&  
         0.2140686E+03_r8,     0.2127284E+03_r8,     0.2113946E+03_r8,     0.2100673E+03_r8,     0.2087468E+03_r8,&  
         0.2074332E+03_r8,     0.2061266E+03_r8,     0.2048272E+03_r8,     0.2035350E+03_r8,     0.2022502E+03_r8,&  
         0.2009727E+03_r8,     0.1997027E+03_r8,     0.1984402E+03_r8,     0.1971852E+03_r8,     0.1959377E+03_r8,&  
         0.1946978E+03_r8,     0.1934654E+03_r8,     0.1922405E+03_r8,     0.1910232E+03_r8,     0.1898134E+03_r8,&  
         0.1886111E+03_r8,     0.1874163E+03_r8,     0.1862289E+03_r8,     0.1850490E+03_r8,     0.1838764E+03_r8,&  
         0.1827112E+03_r8,     0.1815533E+03_r8,     0.1804028E+03_r8,     0.1792594E+03_r8,     0.1781233E+03_r8,&  
         0.1769943E+03_r8,     0.1758725E+03_r8,     0.1747577E+03_r8,     0.1736500E+03_r8,     0.1725493E+03_r8,&  
         0.1714556E+03_r8,     0.1703688E+03_r8,     0.1692889E+03_r8,     0.1682158E+03_r8,     0.1671495E+03_r8,&  
         0.1660900E+03_r8,     0.1650371E+03_r8,     0.1639910E+03_r8,     0.1629515E+03_r8,     0.1619185E+03_r8,&  
         0.1608921E+03_r8,     0.1598723E+03_r8,     0.1588588E+03_r8,     0.1578519E+03_r8,     0.1568512E+03_r8,&  
         0.1558569E+03_r8,     0.1548690E+03_r8,     0.1538873E+03_r8,     0.1529118E+03_r8,     0.1519425E+03_r8,&  
         0.1509793E+03_r8,     0.1500223E+03_r8,     0.1490713E+03_r8,     0.1481263E+03_r8,     0.1471874E+03_r8,&  
         0.1462543E+03_r8,     0.1453272E+03_r8,     0.1444060E+03_r8,     0.1434906E+03_r8,     0.1425810E+03_r8,&  
         0.1416772E+03_r8,     0.1407791E+03_r8,     0.1398867E+03_r8,     0.1390000E+03_r8,     0.1381189E+03_r8,&  
         0.1372433E+03_r8,     0.1363734E+03_r8,     0.1355089E+03_r8,     0.1346499E+03_r8,     0.1337964E+03_r8,&  
         0.1329482E+03_r8,     0.1321055E+03_r8,     0.1312681E+03_r8,     0.1304360E+03_r8,     0.1296091E+03_r8,&  
         0.1287875E+03_r8,     0.1279712E+03_r8,     0.1271600E+03_r8,     0.1263539E+03_r8,     0.1255529E+03_r8,&  
         0.1247571E+03_r8,     0.1239662E+03_r8,     0.1231804E+03_r8,     0.1223996E+03_r8,     0.1216237E+03_r8/)  

    psaditmk(1:150,   53)= (/ &
         0.2884430E+03_r8,     0.2875147E+03_r8,     0.2865791E+03_r8,     0.2856365E+03_r8,     0.2846863E+03_r8,&  
         0.2837283E+03_r8,     0.2827620E+03_r8,     0.2817874E+03_r8,     0.2808039E+03_r8,     0.2798114E+03_r8,&  
         0.2788095E+03_r8,     0.2777979E+03_r8,     0.2767762E+03_r8,     0.2757442E+03_r8,     0.2747016E+03_r8,&  
         0.2736480E+03_r8,     0.2725831E+03_r8,     0.2715067E+03_r8,     0.2704183E+03_r8,     0.2693181E+03_r8,&  
         0.2682053E+03_r8,     0.2670801E+03_r8,     0.2659421E+03_r8,     0.2647911E+03_r8,     0.2636271E+03_r8,&  
         0.2624499E+03_r8,     0.2612595E+03_r8,     0.2600558E+03_r8,     0.2588390E+03_r8,     0.2576091E+03_r8,&  
         0.2563663E+03_r8,     0.2551108E+03_r8,     0.2538429E+03_r8,     0.2525628E+03_r8,     0.2512711E+03_r8,&  
         0.2499683E+03_r8,     0.2486547E+03_r8,     0.2473310E+03_r8,     0.2459979E+03_r8,     0.2446561E+03_r8,&  
         0.2433062E+03_r8,     0.2419491E+03_r8,     0.2405855E+03_r8,     0.2392162E+03_r8,     0.2378422E+03_r8,&  
         0.2364641E+03_r8,     0.2350829E+03_r8,     0.2336994E+03_r8,     0.2323143E+03_r8,     0.2309284E+03_r8,&  
         0.2295425E+03_r8,     0.2281573E+03_r8,     0.2267735E+03_r8,     0.2253917E+03_r8,     0.2240125E+03_r8,&  
         0.2226365E+03_r8,     0.2212641E+03_r8,     0.2198960E+03_r8,     0.2185325E+03_r8,     0.2171739E+03_r8,&  
         0.2158207E+03_r8,     0.2144732E+03_r8,     0.2131317E+03_r8,     0.2117964E+03_r8,     0.2104675E+03_r8,&  
         0.2091452E+03_r8,     0.2078298E+03_r8,     0.2065214E+03_r8,     0.2052200E+03_r8,     0.2039258E+03_r8,&  
         0.2026389E+03_r8,     0.2013593E+03_r8,     0.2000872E+03_r8,     0.1988225E+03_r8,     0.1975653E+03_r8,&  
         0.1963156E+03_r8,     0.1950734E+03_r8,     0.1938388E+03_r8,     0.1926117E+03_r8,     0.1913922E+03_r8,&  
         0.1901801E+03_r8,     0.1889756E+03_r8,     0.1877785E+03_r8,     0.1865889E+03_r8,     0.1854067E+03_r8,&  
         0.1842319E+03_r8,     0.1830645E+03_r8,     0.1819044E+03_r8,     0.1807516E+03_r8,     0.1796061E+03_r8,&  
         0.1784678E+03_r8,     0.1773366E+03_r8,     0.1762126E+03_r8,     0.1750957E+03_r8,     0.1739859E+03_r8,&  
         0.1728831E+03_r8,     0.1717872E+03_r8,     0.1706984E+03_r8,     0.1696163E+03_r8,     0.1685412E+03_r8,&  
         0.1674728E+03_r8,     0.1664112E+03_r8,     0.1653564E+03_r8,     0.1643082E+03_r8,     0.1632667E+03_r8,&  
         0.1622317E+03_r8,     0.1612034E+03_r8,     0.1601815E+03_r8,     0.1591661E+03_r8,     0.1581572E+03_r8,&  
         0.1571546E+03_r8,     0.1561584E+03_r8,     0.1551685E+03_r8,     0.1541849E+03_r8,     0.1532076E+03_r8,&  
         0.1522364E+03_r8,     0.1512714E+03_r8,     0.1503125E+03_r8,     0.1493596E+03_r8,     0.1484129E+03_r8,&  
         0.1474721E+03_r8,     0.1465372E+03_r8,     0.1456083E+03_r8,     0.1446853E+03_r8,     0.1437682E+03_r8,&  
         0.1428568E+03_r8,     0.1419513E+03_r8,     0.1410515E+03_r8,     0.1401573E+03_r8,     0.1392689E+03_r8,&  
         0.1383861E+03_r8,     0.1375088E+03_r8,     0.1366371E+03_r8,     0.1357710E+03_r8,     0.1349104E+03_r8,&  
         0.1340552E+03_r8,     0.1332054E+03_r8,     0.1323610E+03_r8,     0.1315220E+03_r8,     0.1306883E+03_r8,&  
         0.1298598E+03_r8,     0.1290367E+03_r8,     0.1282187E+03_r8,     0.1274059E+03_r8,     0.1265983E+03_r8,&  
         0.1257958E+03_r8,     0.1249984E+03_r8,     0.1242060E+03_r8,     0.1234187E+03_r8,     0.1226363E+03_r8/)  

    psaditmk(1:150,   54)= (/ &
         0.2893380E+03_r8,     0.2884240E+03_r8,     0.2875036E+03_r8,     0.2865759E+03_r8,     0.2856413E+03_r8,&  
         0.2846994E+03_r8,     0.2837497E+03_r8,     0.2827921E+03_r8,     0.2818260E+03_r8,     0.2808514E+03_r8,&  
         0.2798679E+03_r8,     0.2788752E+03_r8,     0.2778730E+03_r8,     0.2768608E+03_r8,     0.2758386E+03_r8,&  
         0.2748057E+03_r8,     0.2737621E+03_r8,     0.2727074E+03_r8,     0.2716412E+03_r8,     0.2705635E+03_r8,&  
         0.2694736E+03_r8,     0.2683716E+03_r8,     0.2672570E+03_r8,     0.2661297E+03_r8,     0.2649897E+03_r8,&  
         0.2638364E+03_r8,     0.2626701E+03_r8,     0.2614905E+03_r8,     0.2602976E+03_r8,     0.2590915E+03_r8,&  
         0.2578722E+03_r8,     0.2566398E+03_r8,     0.2553945E+03_r8,     0.2541365E+03_r8,     0.2528662E+03_r8,&  
         0.2515839E+03_r8,     0.2502901E+03_r8,     0.2489853E+03_r8,     0.2476699E+03_r8,     0.2463447E+03_r8,&  
         0.2450103E+03_r8,     0.2436675E+03_r8,     0.2423168E+03_r8,     0.2409592E+03_r8,     0.2395954E+03_r8,&  
         0.2382263E+03_r8,     0.2368527E+03_r8,     0.2354754E+03_r8,     0.2340954E+03_r8,     0.2327132E+03_r8,&  
         0.2313299E+03_r8,     0.2299461E+03_r8,     0.2285626E+03_r8,     0.2271801E+03_r8,     0.2257993E+03_r8,&  
         0.2244207E+03_r8,     0.2230449E+03_r8,     0.2216726E+03_r8,     0.2203042E+03_r8,     0.2189401E+03_r8,&  
         0.2175808E+03_r8,     0.2162267E+03_r8,     0.2148781E+03_r8,     0.2135352E+03_r8,     0.2121985E+03_r8,&  
         0.2108680E+03_r8,     0.2095442E+03_r8,     0.2082270E+03_r8,     0.2069166E+03_r8,     0.2056134E+03_r8,&  
         0.2043172E+03_r8,     0.2030282E+03_r8,     0.2017466E+03_r8,     0.2004723E+03_r8,     0.1992054E+03_r8,&  
         0.1979461E+03_r8,     0.1966942E+03_r8,     0.1954498E+03_r8,     0.1942130E+03_r8,     0.1929837E+03_r8,&  
         0.1917618E+03_r8,     0.1905475E+03_r8,     0.1893407E+03_r8,     0.1881414E+03_r8,     0.1869496E+03_r8,&  
         0.1857652E+03_r8,     0.1845882E+03_r8,     0.1834185E+03_r8,     0.1822562E+03_r8,     0.1811012E+03_r8,&  
         0.1799535E+03_r8,     0.1788130E+03_r8,     0.1776797E+03_r8,     0.1765535E+03_r8,     0.1754345E+03_r8,&  
         0.1743225E+03_r8,     0.1732176E+03_r8,     0.1721196E+03_r8,     0.1710286E+03_r8,     0.1699445E+03_r8,&  
         0.1688673E+03_r8,     0.1677969E+03_r8,     0.1667332E+03_r8,     0.1656763E+03_r8,     0.1646261E+03_r8,&  
         0.1635826E+03_r8,     0.1625457E+03_r8,     0.1615153E+03_r8,     0.1604915E+03_r8,     0.1594741E+03_r8,&  
         0.1584632E+03_r8,     0.1574587E+03_r8,     0.1564606E+03_r8,     0.1554688E+03_r8,     0.1544833E+03_r8,&  
         0.1535040E+03_r8,     0.1525310E+03_r8,     0.1515641E+03_r8,     0.1506033E+03_r8,     0.1496486E+03_r8,&  
         0.1487000E+03_r8,     0.1477574E+03_r8,     0.1468208E+03_r8,     0.1458901E+03_r8,     0.1449653E+03_r8,&  
         0.1440464E+03_r8,     0.1431333E+03_r8,     0.1422259E+03_r8,     0.1413244E+03_r8,     0.1404285E+03_r8,&  
         0.1395384E+03_r8,     0.1386538E+03_r8,     0.1377749E+03_r8,     0.1369016E+03_r8,     0.1360337E+03_r8,&  
         0.1351714E+03_r8,     0.1343146E+03_r8,     0.1334631E+03_r8,     0.1326171E+03_r8,     0.1317765E+03_r8,&  
         0.1309411E+03_r8,     0.1301111E+03_r8,     0.1292863E+03_r8,     0.1284668E+03_r8,     0.1276524E+03_r8,&  
         0.1268433E+03_r8,     0.1260392E+03_r8,     0.1252402E+03_r8,     0.1244464E+03_r8,     0.1236575E+03_r8/)  

    psaditmk(1:150,   55)= (/ &
         0.2902071E+03_r8,     0.2893069E+03_r8,     0.2884004E+03_r8,     0.2874875E+03_r8,     0.2865676E+03_r8,&  
         0.2856411E+03_r8,     0.2847072E+03_r8,     0.2837658E+03_r8,     0.2828165E+03_r8,     0.2818591E+03_r8,&  
         0.2808933E+03_r8,     0.2799188E+03_r8,     0.2789351E+03_r8,     0.2779421E+03_r8,     0.2769394E+03_r8,&  
         0.2759267E+03_r8,     0.2749036E+03_r8,     0.2738698E+03_r8,     0.2728251E+03_r8,     0.2717692E+03_r8,&  
         0.2707017E+03_r8,     0.2696223E+03_r8,     0.2685309E+03_r8,     0.2674269E+03_r8,     0.2663104E+03_r8,&  
         0.2651812E+03_r8,     0.2640388E+03_r8,     0.2628832E+03_r8,     0.2617146E+03_r8,     0.2605325E+03_r8,&  
         0.2593372E+03_r8,     0.2581286E+03_r8,     0.2569066E+03_r8,     0.2556715E+03_r8,     0.2544237E+03_r8,&  
         0.2531633E+03_r8,     0.2518907E+03_r8,     0.2506062E+03_r8,     0.2493103E+03_r8,     0.2480034E+03_r8,&  
         0.2466864E+03_r8,     0.2453597E+03_r8,     0.2440240E+03_r8,     0.2426801E+03_r8,     0.2413288E+03_r8,&  
         0.2399708E+03_r8,     0.2386070E+03_r8,     0.2372381E+03_r8,     0.2358651E+03_r8,     0.2344888E+03_r8,&  
         0.2331100E+03_r8,     0.2317294E+03_r8,     0.2303480E+03_r8,     0.2289664E+03_r8,     0.2275854E+03_r8,&  
         0.2262057E+03_r8,     0.2248279E+03_r8,     0.2234526E+03_r8,     0.2220804E+03_r8,     0.2207119E+03_r8,&  
         0.2193474E+03_r8,     0.2179875E+03_r8,     0.2166325E+03_r8,     0.2152829E+03_r8,     0.2139388E+03_r8,&  
         0.2126007E+03_r8,     0.2112688E+03_r8,     0.2099433E+03_r8,     0.2086244E+03_r8,     0.2073123E+03_r8,&  
         0.2060071E+03_r8,     0.2047090E+03_r8,     0.2034181E+03_r8,     0.2021344E+03_r8,     0.2008580E+03_r8,&  
         0.1995890E+03_r8,     0.1983274E+03_r8,     0.1970734E+03_r8,     0.1958268E+03_r8,     0.1945877E+03_r8,&  
         0.1933561E+03_r8,     0.1921321E+03_r8,     0.1909156E+03_r8,     0.1897065E+03_r8,     0.1885050E+03_r8,&  
         0.1873109E+03_r8,     0.1861242E+03_r8,     0.1849450E+03_r8,     0.1837731E+03_r8,     0.1826086E+03_r8,&  
         0.1814514E+03_r8,     0.1803015E+03_r8,     0.1791588E+03_r8,     0.1780233E+03_r8,     0.1768950E+03_r8,&  
         0.1757738E+03_r8,     0.1746597E+03_r8,     0.1735526E+03_r8,     0.1724525E+03_r8,     0.1713594E+03_r8,&  
         0.1702733E+03_r8,     0.1691939E+03_r8,     0.1681214E+03_r8,     0.1670558E+03_r8,     0.1659968E+03_r8,&  
         0.1649446E+03_r8,     0.1638990E+03_r8,     0.1628601E+03_r8,     0.1618277E+03_r8,     0.1608019E+03_r8,&  
         0.1597826E+03_r8,     0.1587697E+03_r8,     0.1577633E+03_r8,     0.1567632E+03_r8,     0.1557695E+03_r8,&  
         0.1547821E+03_r8,     0.1538010E+03_r8,     0.1528260E+03_r8,     0.1518573E+03_r8,     0.1508947E+03_r8,&  
         0.1499381E+03_r8,     0.1489877E+03_r8,     0.1480432E+03_r8,     0.1471048E+03_r8,     0.1461723E+03_r8,&  
         0.1452457E+03_r8,     0.1443250E+03_r8,     0.1434101E+03_r8,     0.1425011E+03_r8,     0.1415977E+03_r8,&  
         0.1407002E+03_r8,     0.1398083E+03_r8,     0.1389220E+03_r8,     0.1380414E+03_r8,     0.1371664E+03_r8,&  
         0.1362969E+03_r8,     0.1354329E+03_r8,     0.1345744E+03_r8,     0.1337213E+03_r8,     0.1328737E+03_r8,&  
         0.1320314E+03_r8,     0.1311944E+03_r8,     0.1303628E+03_r8,     0.1295364E+03_r8,     0.1287153E+03_r8,&  
         0.1278994E+03_r8,     0.1270886E+03_r8,     0.1262830E+03_r8,     0.1254825E+03_r8,     0.1246871E+03_r8/)  

    psaditmk(1:150,   56)= (/ &
         0.2910515E+03_r8,     0.2901642E+03_r8,     0.2892712E+03_r8,     0.2883721E+03_r8,     0.2874668E+03_r8,&  
         0.2865545E+03_r8,     0.2856360E+03_r8,     0.2847100E+03_r8,     0.2837767E+03_r8,     0.2828357E+03_r8,&  
         0.2818868E+03_r8,     0.2809296E+03_r8,     0.2799638E+03_r8,     0.2789891E+03_r8,     0.2780052E+03_r8,&  
         0.2770118E+03_r8,     0.2760085E+03_r8,     0.2749952E+03_r8,     0.2739711E+03_r8,     0.2729363E+03_r8,&  
         0.2718905E+03_r8,     0.2708333E+03_r8,     0.2697644E+03_r8,     0.2686833E+03_r8,     0.2675901E+03_r8,&  
         0.2664843E+03_r8,     0.2653657E+03_r8,     0.2642342E+03_r8,     0.2630896E+03_r8,     0.2619317E+03_r8,&  
         0.2607605E+03_r8,     0.2595760E+03_r8,     0.2583779E+03_r8,     0.2571666E+03_r8,     0.2559421E+03_r8,&  
         0.2547046E+03_r8,     0.2534542E+03_r8,     0.2521914E+03_r8,     0.2509162E+03_r8,     0.2496296E+03_r8,&  
         0.2483316E+03_r8,     0.2470228E+03_r8,     0.2457041E+03_r8,     0.2443760E+03_r8,     0.2430391E+03_r8,&  
         0.2416943E+03_r8,     0.2403424E+03_r8,     0.2389841E+03_r8,     0.2376204E+03_r8,     0.2362519E+03_r8,&  
         0.2348796E+03_r8,     0.2335043E+03_r8,     0.2321268E+03_r8,     0.2307480E+03_r8,     0.2293686E+03_r8,&  
         0.2279893E+03_r8,     0.2266110E+03_r8,     0.2252341E+03_r8,     0.2238595E+03_r8,     0.2224876E+03_r8,&  
         0.2211191E+03_r8,     0.2197543E+03_r8,     0.2183939E+03_r8,     0.2170383E+03_r8,     0.2156877E+03_r8,&  
         0.2143425E+03_r8,     0.2130032E+03_r8,     0.2116698E+03_r8,     0.2103428E+03_r8,     0.2090222E+03_r8,&  
         0.2077083E+03_r8,     0.2064013E+03_r8,     0.2051013E+03_r8,     0.2038083E+03_r8,     0.2025226E+03_r8,&  
         0.2012442E+03_r8,     0.1999731E+03_r8,     0.1987094E+03_r8,     0.1974531E+03_r8,     0.1962044E+03_r8,&  
         0.1949631E+03_r8,     0.1937293E+03_r8,     0.1925030E+03_r8,     0.1912842E+03_r8,     0.1900730E+03_r8,&  
         0.1888692E+03_r8,     0.1876729E+03_r8,     0.1864840E+03_r8,     0.1853025E+03_r8,     0.1841284E+03_r8,&  
         0.1829617E+03_r8,     0.1818023E+03_r8,     0.1806501E+03_r8,     0.1795053E+03_r8,     0.1783676E+03_r8,&  
         0.1772371E+03_r8,     0.1761138E+03_r8,     0.1749975E+03_r8,     0.1738883E+03_r8,     0.1727861E+03_r8,&  
         0.1716909E+03_r8,     0.1706026E+03_r8,     0.1695212E+03_r8,     0.1684467E+03_r8,     0.1673789E+03_r8,&  
         0.1663179E+03_r8,     0.1652637E+03_r8,     0.1642161E+03_r8,     0.1631751E+03_r8,     0.1621408E+03_r8,&  
         0.1611130E+03_r8,     0.1600917E+03_r8,     0.1590769E+03_r8,     0.1580685E+03_r8,     0.1570665E+03_r8,&  
         0.1560709E+03_r8,     0.1550815E+03_r8,     0.1540985E+03_r8,     0.1531216E+03_r8,     0.1521510E+03_r8,&  
         0.1511865E+03_r8,     0.1502282E+03_r8,     0.1492759E+03_r8,     0.1483296E+03_r8,     0.1473894E+03_r8,&  
         0.1464551E+03_r8,     0.1455267E+03_r8,     0.1446042E+03_r8,     0.1436876E+03_r8,     0.1427767E+03_r8,&  
         0.1418717E+03_r8,     0.1409724E+03_r8,     0.1400787E+03_r8,     0.1391908E+03_r8,     0.1383084E+03_r8,&  
         0.1374317E+03_r8,     0.1365605E+03_r8,     0.1356949E+03_r8,     0.1348347E+03_r8,     0.1339800E+03_r8,&  
         0.1331307E+03_r8,     0.1322868E+03_r8,     0.1314482E+03_r8,     0.1306150E+03_r8,     0.1297870E+03_r8,&  
         0.1289643E+03_r8,     0.1281468E+03_r8,     0.1273345E+03_r8,     0.1265273E+03_r8,     0.1257253E+03_r8/)  

    psaditmk(1:150,   57)= (/ &
         0.2918719E+03_r8,     0.2909971E+03_r8,     0.2901168E+03_r8,     0.2892310E+03_r8,     0.2883392E+03_r8,&  
         0.2874413E+03_r8,     0.2865366E+03_r8,     0.2856259E+03_r8,     0.2847077E+03_r8,     0.2837824E+03_r8,&  
         0.2828498E+03_r8,     0.2819091E+03_r8,     0.2809604E+03_r8,     0.2800033E+03_r8,     0.2790374E+03_r8,&  
         0.2780626E+03_r8,     0.2770783E+03_r8,     0.2760844E+03_r8,     0.2750805E+03_r8,     0.2740663E+03_r8,&  
         0.2730414E+03_r8,     0.2720056E+03_r8,     0.2709586E+03_r8,     0.2698997E+03_r8,     0.2688292E+03_r8,&  
         0.2677464E+03_r8,     0.2666512E+03_r8,     0.2655434E+03_r8,     0.2644226E+03_r8,     0.2632889E+03_r8,&  
         0.2621419E+03_r8,     0.2609815E+03_r8,     0.2598078E+03_r8,     0.2586206E+03_r8,     0.2574200E+03_r8,&  
         0.2562061E+03_r8,     0.2549789E+03_r8,     0.2537387E+03_r8,     0.2524858E+03_r8,     0.2512204E+03_r8,&  
         0.2499431E+03_r8,     0.2486541E+03_r8,     0.2473540E+03_r8,     0.2460435E+03_r8,     0.2447232E+03_r8,&  
         0.2433937E+03_r8,     0.2420557E+03_r8,     0.2407101E+03_r8,     0.2393577E+03_r8,     0.2379992E+03_r8,&  
         0.2366356E+03_r8,     0.2352676E+03_r8,     0.2338961E+03_r8,     0.2325220E+03_r8,     0.2311460E+03_r8,&  
         0.2297690E+03_r8,     0.2283918E+03_r8,     0.2270150E+03_r8,     0.2256393E+03_r8,     0.2242655E+03_r8,&  
         0.2228941E+03_r8,     0.2215257E+03_r8,     0.2201609E+03_r8,     0.2188001E+03_r8,     0.2174438E+03_r8,&  
         0.2160924E+03_r8,     0.2147463E+03_r8,     0.2134057E+03_r8,     0.2120710E+03_r8,     0.2107425E+03_r8,&  
         0.2094204E+03_r8,     0.2081048E+03_r8,     0.2067960E+03_r8,     0.2054941E+03_r8,     0.2041992E+03_r8,&  
         0.2029115E+03_r8,     0.2016309E+03_r8,     0.2003578E+03_r8,     0.1990919E+03_r8,     0.1978335E+03_r8,&  
         0.1965826E+03_r8,     0.1953391E+03_r8,     0.1941031E+03_r8,     0.1928746E+03_r8,     0.1916536E+03_r8,&  
         0.1904401E+03_r8,     0.1892340E+03_r8,     0.1880355E+03_r8,     0.1868444E+03_r8,     0.1856607E+03_r8,&  
         0.1844843E+03_r8,     0.1833154E+03_r8,     0.1821538E+03_r8,     0.1809995E+03_r8,     0.1798524E+03_r8,&  
         0.1787125E+03_r8,     0.1775799E+03_r8,     0.1764544E+03_r8,     0.1753360E+03_r8,     0.1742246E+03_r8,&  
         0.1731203E+03_r8,     0.1720230E+03_r8,     0.1709326E+03_r8,     0.1698491E+03_r8,     0.1687725E+03_r8,&  
         0.1677027E+03_r8,     0.1666396E+03_r8,     0.1655833E+03_r8,     0.1645337E+03_r8,     0.1634908E+03_r8,&  
         0.1624544E+03_r8,     0.1614246E+03_r8,     0.1604014E+03_r8,     0.1593846E+03_r8,     0.1583743E+03_r8,&  
         0.1573703E+03_r8,     0.1563728E+03_r8,     0.1553815E+03_r8,     0.1543966E+03_r8,     0.1534178E+03_r8,&  
         0.1524453E+03_r8,     0.1514790E+03_r8,     0.1505188E+03_r8,     0.1495646E+03_r8,     0.1486166E+03_r8,&  
         0.1476745E+03_r8,     0.1467384E+03_r8,     0.1458082E+03_r8,     0.1448839E+03_r8,     0.1439655E+03_r8,&  
         0.1430529E+03_r8,     0.1421461E+03_r8,     0.1412450E+03_r8,     0.1403497E+03_r8,     0.1394600E+03_r8,&  
         0.1385760E+03_r8,     0.1376976E+03_r8,     0.1368247E+03_r8,     0.1359574E+03_r8,     0.1350955E+03_r8,&  
         0.1342392E+03_r8,     0.1333882E+03_r8,     0.1325427E+03_r8,     0.1317025E+03_r8,     0.1308676E+03_r8,&  
         0.1300381E+03_r8,     0.1292138E+03_r8,     0.1283947E+03_r8,     0.1275808E+03_r8,     0.1267721E+03_r8/)  

    psaditmk(1:150,   58)= (/ &
         0.2926696E+03_r8,     0.2918066E+03_r8,     0.2909386E+03_r8,     0.2900654E+03_r8,     0.2891866E+03_r8,&  
         0.2883020E+03_r8,     0.2874114E+03_r8,     0.2865143E+03_r8,     0.2856110E+03_r8,     0.2847006E+03_r8,&  
         0.2837831E+03_r8,     0.2828585E+03_r8,     0.2819261E+03_r8,     0.2809858E+03_r8,     0.2800374E+03_r8,&  
         0.2790803E+03_r8,     0.2781143E+03_r8,     0.2771392E+03_r8,     0.2761545E+03_r8,     0.2751600E+03_r8,&  
         0.2741553E+03_r8,     0.2731402E+03_r8,     0.2721142E+03_r8,     0.2710771E+03_r8,     0.2700288E+03_r8,&  
         0.2689684E+03_r8,     0.2678961E+03_r8,     0.2668115E+03_r8,     0.2657143E+03_r8,     0.2646044E+03_r8,&  
         0.2634814E+03_r8,     0.2623452E+03_r8,     0.2611958E+03_r8,     0.2600329E+03_r8,     0.2588565E+03_r8,&  
         0.2576667E+03_r8,     0.2564633E+03_r8,     0.2552468E+03_r8,     0.2540169E+03_r8,     0.2527740E+03_r8,&  
         0.2515185E+03_r8,     0.2502506E+03_r8,     0.2489709E+03_r8,     0.2476797E+03_r8,     0.2463777E+03_r8,&  
         0.2450654E+03_r8,     0.2437435E+03_r8,     0.2424127E+03_r8,     0.2410738E+03_r8,     0.2397275E+03_r8,&  
         0.2383746E+03_r8,     0.2370161E+03_r8,     0.2356527E+03_r8,     0.2342853E+03_r8,     0.2329148E+03_r8,&  
         0.2315420E+03_r8,     0.2301677E+03_r8,     0.2287926E+03_r8,     0.2274176E+03_r8,     0.2260433E+03_r8,&  
         0.2246705E+03_r8,     0.2232998E+03_r8,     0.2219318E+03_r8,     0.2205670E+03_r8,     0.2192060E+03_r8,&  
         0.2178492E+03_r8,     0.2164971E+03_r8,     0.2151500E+03_r8,     0.2138084E+03_r8,     0.2124724E+03_r8,&  
         0.2111425E+03_r8,     0.2098188E+03_r8,     0.2085016E+03_r8,     0.2071911E+03_r8,     0.2058873E+03_r8,&  
         0.2045905E+03_r8,     0.2033008E+03_r8,     0.2020183E+03_r8,     0.2007430E+03_r8,     0.1994751E+03_r8,&  
         0.1982145E+03_r8,     0.1969614E+03_r8,     0.1957157E+03_r8,     0.1944775E+03_r8,     0.1932468E+03_r8,&  
         0.1920236E+03_r8,     0.1908079E+03_r8,     0.1895996E+03_r8,     0.1883988E+03_r8,     0.1872054E+03_r8,&  
         0.1860195E+03_r8,     0.1848410E+03_r8,     0.1836698E+03_r8,     0.1825060E+03_r8,     0.1813495E+03_r8,&  
         0.1802002E+03_r8,     0.1790581E+03_r8,     0.1779234E+03_r8,     0.1767957E+03_r8,     0.1756751E+03_r8,&  
         0.1745616E+03_r8,     0.1734552E+03_r8,     0.1723557E+03_r8,     0.1712632E+03_r8,     0.1701777E+03_r8,&  
         0.1690989E+03_r8,     0.1680271E+03_r8,     0.1669620E+03_r8,     0.1659037E+03_r8,     0.1648520E+03_r8,&  
         0.1638070E+03_r8,     0.1627687E+03_r8,     0.1617369E+03_r8,     0.1607117E+03_r8,     0.1596929E+03_r8,&  
         0.1586806E+03_r8,     0.1576748E+03_r8,     0.1566753E+03_r8,     0.1556821E+03_r8,     0.1546952E+03_r8,&  
         0.1537146E+03_r8,     0.1527402E+03_r8,     0.1517720E+03_r8,     0.1508099E+03_r8,     0.1498540E+03_r8,&  
         0.1489041E+03_r8,     0.1479601E+03_r8,     0.1470222E+03_r8,     0.1460903E+03_r8,     0.1451642E+03_r8,&  
         0.1442440E+03_r8,     0.1433297E+03_r8,     0.1424211E+03_r8,     0.1415183E+03_r8,     0.1406212E+03_r8,&  
         0.1397298E+03_r8,     0.1388441E+03_r8,     0.1379639E+03_r8,     0.1370894E+03_r8,     0.1362204E+03_r8,&  
         0.1353569E+03_r8,     0.1344989E+03_r8,     0.1336463E+03_r8,     0.1327991E+03_r8,     0.1319573E+03_r8,&  
         0.1311208E+03_r8,     0.1302896E+03_r8,     0.1294637E+03_r8,     0.1286431E+03_r8,     0.1278276E+03_r8/)  

    psaditmk(1:150,   59)= (/ &
         0.2934454E+03_r8,     0.2925936E+03_r8,     0.2917375E+03_r8,     0.2908762E+03_r8,     0.2900096E+03_r8,&  
         0.2891378E+03_r8,     0.2882603E+03_r8,     0.2873769E+03_r8,     0.2864872E+03_r8,     0.2855913E+03_r8,&  
         0.2846886E+03_r8,     0.2837791E+03_r8,     0.2828623E+03_r8,     0.2819380E+03_r8,     0.2810060E+03_r8,&  
         0.2800659E+03_r8,     0.2791175E+03_r8,     0.2781603E+03_r8,     0.2771943E+03_r8,     0.2762185E+03_r8,&  
         0.2752334E+03_r8,     0.2742381E+03_r8,     0.2732326E+03_r8,     0.2722165E+03_r8,     0.2711893E+03_r8,&  
         0.2701510E+03_r8,     0.2691010E+03_r8,     0.2680390E+03_r8,     0.2669650E+03_r8,     0.2658784E+03_r8,&  
         0.2647792E+03_r8,     0.2636670E+03_r8,     0.2625417E+03_r8,     0.2614032E+03_r8,     0.2602510E+03_r8,&  
         0.2590856E+03_r8,     0.2579065E+03_r8,     0.2567139E+03_r8,     0.2555078E+03_r8,     0.2542884E+03_r8,&  
         0.2530558E+03_r8,     0.2518103E+03_r8,     0.2505522E+03_r8,     0.2492819E+03_r8,     0.2479998E+03_r8,&  
         0.2467065E+03_r8,     0.2454025E+03_r8,     0.2440884E+03_r8,     0.2427651E+03_r8,     0.2414331E+03_r8,&  
         0.2400932E+03_r8,     0.2387463E+03_r8,     0.2373932E+03_r8,     0.2360348E+03_r8,     0.2346718E+03_r8,&  
         0.2333051E+03_r8,     0.2319357E+03_r8,     0.2305643E+03_r8,     0.2291917E+03_r8,     0.2278187E+03_r8,&  
         0.2264460E+03_r8,     0.2250744E+03_r8,     0.2237046E+03_r8,     0.2223371E+03_r8,     0.2209725E+03_r8,&  
         0.2196114E+03_r8,     0.2182543E+03_r8,     0.2169016E+03_r8,     0.2155537E+03_r8,     0.2142111E+03_r8,&  
         0.2128740E+03_r8,     0.2115427E+03_r8,     0.2102176E+03_r8,     0.2088987E+03_r8,     0.2075865E+03_r8,&  
         0.2062810E+03_r8,     0.2049823E+03_r8,     0.2036907E+03_r8,     0.2024062E+03_r8,     0.2011289E+03_r8,&  
         0.1998588E+03_r8,     0.1985962E+03_r8,     0.1973409E+03_r8,     0.1960930E+03_r8,     0.1948526E+03_r8,&  
         0.1936197E+03_r8,     0.1923943E+03_r8,     0.1911763E+03_r8,     0.1899658E+03_r8,     0.1887628E+03_r8,&  
         0.1875672E+03_r8,     0.1863790E+03_r8,     0.1851982E+03_r8,     0.1840249E+03_r8,     0.1828588E+03_r8,&  
         0.1817001E+03_r8,     0.1805486E+03_r8,     0.1794044E+03_r8,     0.1782674E+03_r8,     0.1771376E+03_r8,&  
         0.1760149E+03_r8,     0.1748992E+03_r8,     0.1737907E+03_r8,     0.1726891E+03_r8,     0.1715945E+03_r8,&  
         0.1705068E+03_r8,     0.1694260E+03_r8,     0.1683521E+03_r8,     0.1672850E+03_r8,     0.1662246E+03_r8,&  
         0.1651709E+03_r8,     0.1641239E+03_r8,     0.1630835E+03_r8,     0.1620497E+03_r8,     0.1610225E+03_r8,&  
         0.1600018E+03_r8,     0.1589876E+03_r8,     0.1579798E+03_r8,     0.1569783E+03_r8,     0.1559833E+03_r8,&  
         0.1549945E+03_r8,     0.1540120E+03_r8,     0.1530357E+03_r8,     0.1520656E+03_r8,     0.1511017E+03_r8,&  
         0.1501439E+03_r8,     0.1491921E+03_r8,     0.1482464E+03_r8,     0.1473066E+03_r8,     0.1463729E+03_r8,&  
         0.1454450E+03_r8,     0.1445230E+03_r8,     0.1436069E+03_r8,     0.1426966E+03_r8,     0.1417920E+03_r8,&  
         0.1408932E+03_r8,     0.1400001E+03_r8,     0.1391126E+03_r8,     0.1382308E+03_r8,     0.1373546E+03_r8,&  
         0.1364839E+03_r8,     0.1356187E+03_r8,     0.1347590E+03_r8,     0.1339048E+03_r8,     0.1330560E+03_r8,&  
         0.1322125E+03_r8,     0.1313745E+03_r8,     0.1305417E+03_r8,     0.1297142E+03_r8,     0.1288919E+03_r8/)  

    psaditmk(1:150,   60)= (/ &
         0.2942003E+03_r8,     0.2933594E+03_r8,     0.2925143E+03_r8,     0.2916643E+03_r8,     0.2908097E+03_r8,&  
         0.2899499E+03_r8,     0.2890849E+03_r8,     0.2882143E+03_r8,     0.2873381E+03_r8,     0.2864555E+03_r8,&  
         0.2855671E+03_r8,     0.2846718E+03_r8,     0.2837701E+03_r8,     0.2828612E+03_r8,     0.2819449E+03_r8,&  
         0.2810211E+03_r8,     0.2800893E+03_r8,     0.2791493E+03_r8,     0.2782008E+03_r8,     0.2772433E+03_r8,&  
         0.2762768E+03_r8,     0.2753009E+03_r8,     0.2743150E+03_r8,     0.2733190E+03_r8,     0.2723125E+03_r8,&  
         0.2712951E+03_r8,     0.2702671E+03_r8,     0.2692270E+03_r8,     0.2681754E+03_r8,     0.2671118E+03_r8,&  
         0.2660358E+03_r8,     0.2649473E+03_r8,     0.2638457E+03_r8,     0.2627313E+03_r8,     0.2616035E+03_r8,&  
         0.2604624E+03_r8,     0.2593079E+03_r8,     0.2581395E+03_r8,     0.2569577E+03_r8,     0.2557623E+03_r8,&  
         0.2545534E+03_r8,     0.2533312E+03_r8,     0.2520959E+03_r8,     0.2508477E+03_r8,     0.2495870E+03_r8,&  
         0.2483142E+03_r8,     0.2470298E+03_r8,     0.2457344E+03_r8,     0.2444285E+03_r8,     0.2431128E+03_r8,&  
         0.2417879E+03_r8,     0.2404548E+03_r8,     0.2391142E+03_r8,     0.2377668E+03_r8,     0.2364135E+03_r8,&  
         0.2350552E+03_r8,     0.2336927E+03_r8,     0.2323270E+03_r8,     0.2309588E+03_r8,     0.2295889E+03_r8,&  
         0.2282181E+03_r8,     0.2268473E+03_r8,     0.2254772E+03_r8,     0.2241084E+03_r8,     0.2227416E+03_r8,&  
         0.2213773E+03_r8,     0.2200163E+03_r8,     0.2186590E+03_r8,     0.2173058E+03_r8,     0.2159573E+03_r8,&  
         0.2146137E+03_r8,     0.2132755E+03_r8,     0.2119431E+03_r8,     0.2106165E+03_r8,     0.2092962E+03_r8,&  
         0.2079823E+03_r8,     0.2066750E+03_r8,     0.2053746E+03_r8,     0.2040810E+03_r8,     0.2027946E+03_r8,&  
         0.2015152E+03_r8,     0.2002432E+03_r8,     0.1989784E+03_r8,     0.1977209E+03_r8,     0.1964709E+03_r8,&  
         0.1952283E+03_r8,     0.1939932E+03_r8,     0.1927656E+03_r8,     0.1915454E+03_r8,     0.1903326E+03_r8,&  
         0.1891274E+03_r8,     0.1879296E+03_r8,     0.1867392E+03_r8,     0.1855562E+03_r8,     0.1843806E+03_r8,&  
         0.1832123E+03_r8,     0.1820514E+03_r8,     0.1808978E+03_r8,     0.1797513E+03_r8,     0.1786122E+03_r8,&  
         0.1774802E+03_r8,     0.1763553E+03_r8,     0.1752375E+03_r8,     0.1741268E+03_r8,     0.1730231E+03_r8,&  
         0.1719264E+03_r8,     0.1708367E+03_r8,     0.1697538E+03_r8,     0.1686777E+03_r8,     0.1676085E+03_r8,&  
         0.1665461E+03_r8,     0.1654904E+03_r8,     0.1644414E+03_r8,     0.1633990E+03_r8,     0.1623632E+03_r8,&  
         0.1613340E+03_r8,     0.1603113E+03_r8,     0.1592951E+03_r8,     0.1582854E+03_r8,     0.1572820E+03_r8,&  
         0.1562850E+03_r8,     0.1552943E+03_r8,     0.1543099E+03_r8,     0.1533317E+03_r8,     0.1523598E+03_r8,&  
         0.1513940E+03_r8,     0.1504343E+03_r8,     0.1494807E+03_r8,     0.1485331E+03_r8,     0.1475916E+03_r8,&  
         0.1466560E+03_r8,     0.1457264E+03_r8,     0.1448026E+03_r8,     0.1438847E+03_r8,     0.1429726E+03_r8,&  
         0.1420663E+03_r8,     0.1411658E+03_r8,     0.1402709E+03_r8,     0.1393818E+03_r8,     0.1384982E+03_r8,&  
         0.1376203E+03_r8,     0.1367479E+03_r8,     0.1358811E+03_r8,     0.1350197E+03_r8,     0.1341638E+03_r8,&  
         0.1333134E+03_r8,     0.1324683E+03_r8,     0.1316286E+03_r8,     0.1307942E+03_r8,     0.1299651E+03_r8/)  

    psaditmk(1:150,   61)= (/ &
         0.2949350E+03_r8,     0.2941045E+03_r8,     0.2932699E+03_r8,     0.2924309E+03_r8,     0.2915874E+03_r8,&  
         0.2907392E+03_r8,     0.2898861E+03_r8,     0.2890279E+03_r8,     0.2881642E+03_r8,     0.2872947E+03_r8,&  
         0.2864192E+03_r8,     0.2855383E+03_r8,     0.2846506E+03_r8,     0.2837563E+03_r8,     0.2828551E+03_r8,&  
         0.2819468E+03_r8,     0.2810310E+03_r8,     0.2801074E+03_r8,     0.2791758E+03_r8,     0.2782357E+03_r8,&  
         0.2772872E+03_r8,     0.2763295E+03_r8,     0.2753626E+03_r8,     0.2743860E+03_r8,     0.2733994E+03_r8,&  
         0.2724025E+03_r8,     0.2713949E+03_r8,     0.2703764E+03_r8,     0.2693467E+03_r8,     0.2683053E+03_r8,&  
         0.2672520E+03_r8,     0.2661865E+03_r8,     0.2651084E+03_r8,     0.2640177E+03_r8,     0.2629140E+03_r8,&  
         0.2617971E+03_r8,     0.2606668E+03_r8,     0.2595231E+03_r8,     0.2583657E+03_r8,     0.2571947E+03_r8,&  
         0.2560100E+03_r8,     0.2548118E+03_r8,     0.2536001E+03_r8,     0.2523750E+03_r8,     0.2511369E+03_r8,&  
         0.2498860E+03_r8,     0.2486228E+03_r8,     0.2473474E+03_r8,     0.2460609E+03_r8,     0.2447634E+03_r8,&  
         0.2434557E+03_r8,     0.2421383E+03_r8,     0.2408122E+03_r8,     0.2394780E+03_r8,     0.2381367E+03_r8,&  
         0.2367889E+03_r8,     0.2354355E+03_r8,     0.2340775E+03_r8,     0.2327157E+03_r8,     0.2313509E+03_r8,&  
         0.2299840E+03_r8,     0.2286158E+03_r8,     0.2272471E+03_r8,     0.2258786E+03_r8,     0.2245110E+03_r8,&  
         0.2231451E+03_r8,     0.2217815E+03_r8,     0.2204207E+03_r8,     0.2190633E+03_r8,     0.2177098E+03_r8,&  
         0.2163607E+03_r8,     0.2150163E+03_r8,     0.2136772E+03_r8,     0.2123435E+03_r8,     0.2110157E+03_r8,&  
         0.2096939E+03_r8,     0.2083784E+03_r8,     0.2070695E+03_r8,     0.2057673E+03_r8,     0.2044718E+03_r8,&  
         0.2031834E+03_r8,     0.2019021E+03_r8,     0.2006280E+03_r8,     0.1993611E+03_r8,     0.1981016E+03_r8,&  
         0.1968494E+03_r8,     0.1956047E+03_r8,     0.1943674E+03_r8,     0.1931375E+03_r8,     0.1919151E+03_r8,&  
         0.1907002E+03_r8,     0.1894927E+03_r8,     0.1882926E+03_r8,     0.1871000E+03_r8,     0.1859148E+03_r8,&  
         0.1847370E+03_r8,     0.1835665E+03_r8,     0.1824034E+03_r8,     0.1812475E+03_r8,     0.1800989E+03_r8,&  
         0.1789576E+03_r8,     0.1778234E+03_r8,     0.1766964E+03_r8,     0.1755764E+03_r8,     0.1744636E+03_r8,&  
         0.1733578E+03_r8,     0.1722590E+03_r8,     0.1711671E+03_r8,     0.1700821E+03_r8,     0.1690040E+03_r8,&  
         0.1679327E+03_r8,     0.1668682E+03_r8,     0.1658105E+03_r8,     0.1647594E+03_r8,     0.1637151E+03_r8,&  
         0.1626773E+03_r8,     0.1616461E+03_r8,     0.1606214E+03_r8,     0.1596033E+03_r8,     0.1585916E+03_r8,&  
         0.1575862E+03_r8,     0.1565873E+03_r8,     0.1555947E+03_r8,     0.1546084E+03_r8,     0.1536283E+03_r8,&  
         0.1526545E+03_r8,     0.1516868E+03_r8,     0.1507253E+03_r8,     0.1497699E+03_r8,     0.1488204E+03_r8,&  
         0.1478771E+03_r8,     0.1469397E+03_r8,     0.1460082E+03_r8,     0.1450827E+03_r8,     0.1441630E+03_r8,&  
         0.1432492E+03_r8,     0.1423411E+03_r8,     0.1414388E+03_r8,     0.1405423E+03_r8,     0.1396514E+03_r8,&  
         0.1387661E+03_r8,     0.1378865E+03_r8,     0.1370124E+03_r8,     0.1361439E+03_r8,     0.1352809E+03_r8,&  
         0.1344234E+03_r8,     0.1335712E+03_r8,     0.1327245E+03_r8,     0.1318832E+03_r8,     0.1310472E+03_r8/)  

    psaditmk(1:150,   62)= (/ &
         0.2956505E+03_r8,     0.2948298E+03_r8,     0.2940053E+03_r8,     0.2931767E+03_r8,     0.2923439E+03_r8,&  
         0.2915068E+03_r8,     0.2906650E+03_r8,     0.2898184E+03_r8,     0.2889669E+03_r8,     0.2881097E+03_r8,&  
         0.2872472E+03_r8,     0.2863791E+03_r8,     0.2855051E+03_r8,     0.2846247E+03_r8,     0.2837379E+03_r8,&  
         0.2828443E+03_r8,     0.2819438E+03_r8,     0.2810359E+03_r8,     0.2801204E+03_r8,     0.2791969E+03_r8,&  
         0.2782656E+03_r8,     0.2773255E+03_r8,     0.2763766E+03_r8,     0.2754185E+03_r8,     0.2744510E+03_r8,&  
         0.2734738E+03_r8,     0.2724863E+03_r8,     0.2714885E+03_r8,     0.2704799E+03_r8,     0.2694599E+03_r8,&  
         0.2684288E+03_r8,     0.2673856E+03_r8,     0.2663304E+03_r8,     0.2652630E+03_r8,     0.2641829E+03_r8,&  
         0.2630899E+03_r8,     0.2619838E+03_r8,     0.2608644E+03_r8,     0.2597314E+03_r8,     0.2585851E+03_r8,&  
         0.2574249E+03_r8,     0.2562510E+03_r8,     0.2550636E+03_r8,     0.2538623E+03_r8,     0.2526476E+03_r8,&  
         0.2514198E+03_r8,     0.2501789E+03_r8,     0.2489253E+03_r8,     0.2476595E+03_r8,     0.2463819E+03_r8,&  
         0.2450930E+03_r8,     0.2437934E+03_r8,     0.2424840E+03_r8,     0.2411651E+03_r8,     0.2398377E+03_r8,&  
         0.2385027E+03_r8,     0.2371606E+03_r8,     0.2358126E+03_r8,     0.2344593E+03_r8,     0.2331017E+03_r8,&  
         0.2317407E+03_r8,     0.2303771E+03_r8,     0.2290116E+03_r8,     0.2276452E+03_r8,     0.2262786E+03_r8,&  
         0.2249125E+03_r8,     0.2235477E+03_r8,     0.2221848E+03_r8,     0.2208244E+03_r8,     0.2194671E+03_r8,&  
         0.2181134E+03_r8,     0.2167639E+03_r8,     0.2154189E+03_r8,     0.2140788E+03_r8,     0.2127441E+03_r8,&  
         0.2114151E+03_r8,     0.2100919E+03_r8,     0.2087749E+03_r8,     0.2074643E+03_r8,     0.2061603E+03_r8,&  
         0.2048631E+03_r8,     0.2035728E+03_r8,     0.2022896E+03_r8,     0.2010134E+03_r8,     0.1997445E+03_r8,&  
         0.1984828E+03_r8,     0.1972285E+03_r8,     0.1959816E+03_r8,     0.1947422E+03_r8,     0.1935101E+03_r8,&  
         0.1922855E+03_r8,     0.1910684E+03_r8,     0.1898586E+03_r8,     0.1886564E+03_r8,     0.1874615E+03_r8,&  
         0.1862741E+03_r8,     0.1850940E+03_r8,     0.1839214E+03_r8,     0.1827560E+03_r8,     0.1815979E+03_r8,&  
         0.1804472E+03_r8,     0.1793036E+03_r8,     0.1781673E+03_r8,     0.1770381E+03_r8,     0.1759160E+03_r8,&  
         0.1748010E+03_r8,     0.1736931E+03_r8,     0.1725921E+03_r8,     0.1714982E+03_r8,     0.1704111E+03_r8,&  
         0.1693309E+03_r8,     0.1682576E+03_r8,     0.1671910E+03_r8,     0.1661312E+03_r8,     0.1650782E+03_r8,&  
         0.1640318E+03_r8,     0.1629920E+03_r8,     0.1619588E+03_r8,     0.1609321E+03_r8,     0.1599120E+03_r8,&  
         0.1588983E+03_r8,     0.1578911E+03_r8,     0.1568902E+03_r8,     0.1558957E+03_r8,     0.1549075E+03_r8,&  
         0.1539255E+03_r8,     0.1529498E+03_r8,     0.1519803E+03_r8,     0.1510169E+03_r8,     0.1500596E+03_r8,&  
         0.1491084E+03_r8,     0.1481631E+03_r8,     0.1472239E+03_r8,     0.1462907E+03_r8,     0.1453634E+03_r8,&  
         0.1444419E+03_r8,     0.1435263E+03_r8,     0.1426165E+03_r8,     0.1417124E+03_r8,     0.1408141E+03_r8,&  
         0.1399215E+03_r8,     0.1390346E+03_r8,     0.1381532E+03_r8,     0.1372775E+03_r8,     0.1364073E+03_r8,&  
         0.1355426E+03_r8,     0.1346834E+03_r8,     0.1338296E+03_r8,     0.1329813E+03_r8,     0.1321383E+03_r8/)  

    psaditmk(1:150,   63)= (/ &
         0.2963476E+03_r8,     0.2955363E+03_r8,     0.2947213E+03_r8,     0.2939027E+03_r8,     0.2930801E+03_r8,&  
         0.2922534E+03_r8,     0.2914224E+03_r8,     0.2905869E+03_r8,     0.2897467E+03_r8,     0.2889015E+03_r8,&  
         0.2880512E+03_r8,     0.2871956E+03_r8,     0.2863342E+03_r8,     0.2854676E+03_r8,     0.2845943E+03_r8,&  
         0.2837149E+03_r8,     0.2828288E+03_r8,     0.2819359E+03_r8,     0.2810358E+03_r8,     0.2801283E+03_r8,&  
         0.2792132E+03_r8,     0.2782900E+03_r8,     0.2773584E+03_r8,     0.2764181E+03_r8,     0.2754689E+03_r8,&  
         0.2745105E+03_r8,     0.2735424E+03_r8,     0.2725642E+03_r8,     0.2715759E+03_r8,     0.2705767E+03_r8,&  
         0.2695669E+03_r8,     0.2685456E+03_r8,     0.2675128E+03_r8,     0.2664680E+03_r8,     0.2654109E+03_r8,&  
         0.2643414E+03_r8,     0.2632590E+03_r8,     0.2621637E+03_r8,     0.2610552E+03_r8,     0.2599331E+03_r8,&  
         0.2587976E+03_r8,     0.2576483E+03_r8,     0.2564853E+03_r8,     0.2553085E+03_r8,     0.2541180E+03_r8,&  
         0.2529139E+03_r8,     0.2516964E+03_r8,     0.2504655E+03_r8,     0.2492218E+03_r8,     0.2479656E+03_r8,&  
         0.2466972E+03_r8,     0.2454172E+03_r8,     0.2441262E+03_r8,     0.2428245E+03_r8,     0.2415134E+03_r8,&  
         0.2401931E+03_r8,     0.2388646E+03_r8,     0.2375287E+03_r8,     0.2361862E+03_r8,     0.2348380E+03_r8,&  
         0.2334850E+03_r8,     0.2321280E+03_r8,     0.2307678E+03_r8,     0.2294054E+03_r8,     0.2280415E+03_r8,&  
         0.2266770E+03_r8,     0.2253127E+03_r8,     0.2239491E+03_r8,     0.2225871E+03_r8,     0.2212273E+03_r8,&  
         0.2198703E+03_r8,     0.2185166E+03_r8,     0.2171668E+03_r8,     0.2158212E+03_r8,     0.2144804E+03_r8,&  
         0.2131447E+03_r8,     0.2118145E+03_r8,     0.2104900E+03_r8,     0.2091716E+03_r8,     0.2078595E+03_r8,&  
         0.2065538E+03_r8,     0.2052549E+03_r8,     0.2039627E+03_r8,     0.2026774E+03_r8,     0.2013994E+03_r8,&  
         0.2001284E+03_r8,     0.1988647E+03_r8,     0.1976083E+03_r8,     0.1963592E+03_r8,     0.1951176E+03_r8,&  
         0.1938833E+03_r8,     0.1926565E+03_r8,     0.1914371E+03_r8,     0.1902252E+03_r8,     0.1890208E+03_r8,&  
         0.1878237E+03_r8,     0.1866340E+03_r8,     0.1854518E+03_r8,     0.1842768E+03_r8,     0.1831093E+03_r8,&  
         0.1819491E+03_r8,     0.1807961E+03_r8,     0.1796503E+03_r8,     0.1785118E+03_r8,     0.1773805E+03_r8,&  
         0.1762562E+03_r8,     0.1751391E+03_r8,     0.1740290E+03_r8,     0.1729260E+03_r8,     0.1718299E+03_r8,&  
         0.1707407E+03_r8,     0.1696584E+03_r8,     0.1685830E+03_r8,     0.1675144E+03_r8,     0.1664526E+03_r8,&  
         0.1653975E+03_r8,     0.1643491E+03_r8,     0.1633073E+03_r8,     0.1622721E+03_r8,     0.1612435E+03_r8,&  
         0.1602213E+03_r8,     0.1592057E+03_r8,     0.1581965E+03_r8,     0.1571937E+03_r8,     0.1561973E+03_r8,&  
         0.1552071E+03_r8,     0.1542233E+03_r8,     0.1532457E+03_r8,     0.1522742E+03_r8,     0.1513090E+03_r8,&  
         0.1503499E+03_r8,     0.1493968E+03_r8,     0.1484498E+03_r8,     0.1475087E+03_r8,     0.1465737E+03_r8,&  
         0.1456445E+03_r8,     0.1447213E+03_r8,     0.1438039E+03_r8,     0.1428924E+03_r8,     0.1419866E+03_r8,&  
         0.1410865E+03_r8,     0.1401922E+03_r8,     0.1393035E+03_r8,     0.1384205E+03_r8,     0.1375430E+03_r8,&  
         0.1366711E+03_r8,     0.1358048E+03_r8,     0.1349439E+03_r8,     0.1340885E+03_r8,     0.1332385E+03_r8/)  

    psaditmk(1:150,   64)= (/ &
         0.2970271E+03_r8,     0.2962245E+03_r8,     0.2954188E+03_r8,     0.2946096E+03_r8,     0.2937967E+03_r8,&  
         0.2929800E+03_r8,     0.2921593E+03_r8,     0.2913344E+03_r8,     0.2905051E+03_r8,     0.2896712E+03_r8,&  
         0.2888324E+03_r8,     0.2879887E+03_r8,     0.2871400E+03_r8,     0.2862856E+03_r8,     0.2854256E+03_r8,&  
         0.2845596E+03_r8,     0.2836874E+03_r8,     0.2828086E+03_r8,     0.2819234E+03_r8,     0.2810310E+03_r8,&  
         0.2801314E+03_r8,     0.2792242E+03_r8,     0.2783092E+03_r8,     0.2773859E+03_r8,     0.2764542E+03_r8,&  
         0.2755138E+03_r8,     0.2745642E+03_r8,     0.2736051E+03_r8,     0.2726363E+03_r8,     0.2716573E+03_r8,&  
         0.2706678E+03_r8,     0.2696676E+03_r8,     0.2686563E+03_r8,     0.2676335E+03_r8,     0.2665989E+03_r8,&  
         0.2655522E+03_r8,     0.2644931E+03_r8,     0.2634215E+03_r8,     0.2623370E+03_r8,     0.2612390E+03_r8,&  
         0.2601280E+03_r8,     0.2590032E+03_r8,     0.2578648E+03_r8,     0.2567127E+03_r8,     0.2555467E+03_r8,&  
         0.2543670E+03_r8,     0.2531735E+03_r8,     0.2519664E+03_r8,     0.2507459E+03_r8,     0.2495122E+03_r8,&  
         0.2482657E+03_r8,     0.2470067E+03_r8,     0.2457358E+03_r8,     0.2444535E+03_r8,     0.2431603E+03_r8,&  
         0.2418568E+03_r8,     0.2405439E+03_r8,     0.2392224E+03_r8,     0.2378928E+03_r8,     0.2365562E+03_r8,&  
         0.2352133E+03_r8,     0.2338652E+03_r8,     0.2325125E+03_r8,     0.2311561E+03_r8,     0.2297970E+03_r8,&  
         0.2284359E+03_r8,     0.2270738E+03_r8,     0.2257114E+03_r8,     0.2243494E+03_r8,     0.2229885E+03_r8,&  
         0.2216295E+03_r8,     0.2202729E+03_r8,     0.2189193E+03_r8,     0.2175693E+03_r8,     0.2162233E+03_r8,&  
         0.2148819E+03_r8,     0.2135453E+03_r8,     0.2122140E+03_r8,     0.2108884E+03_r8,     0.2095686E+03_r8,&  
         0.2082549E+03_r8,     0.2069477E+03_r8,     0.2056469E+03_r8,     0.2043530E+03_r8,     0.2030659E+03_r8,&  
         0.2017858E+03_r8,     0.2005128E+03_r8,     0.1992471E+03_r8,     0.1979886E+03_r8,     0.1967374E+03_r8,&  
         0.1954936E+03_r8,     0.1942572E+03_r8,     0.1930282E+03_r8,     0.1918066E+03_r8,     0.1905925E+03_r8,&  
         0.1893858E+03_r8,     0.1881865E+03_r8,     0.1869947E+03_r8,     0.1858102E+03_r8,     0.1846330E+03_r8,&  
         0.1834633E+03_r8,     0.1823008E+03_r8,     0.1811457E+03_r8,     0.1799977E+03_r8,     0.1788570E+03_r8,&  
         0.1777235E+03_r8,     0.1765971E+03_r8,     0.1754778E+03_r8,     0.1743656E+03_r8,     0.1732605E+03_r8,&  
         0.1721623E+03_r8,     0.1710710E+03_r8,     0.1699866E+03_r8,     0.1689091E+03_r8,     0.1678385E+03_r8,&  
         0.1667746E+03_r8,     0.1657174E+03_r8,     0.1646670E+03_r8,     0.1636232E+03_r8,     0.1625860E+03_r8,&  
         0.1615554E+03_r8,     0.1605313E+03_r8,     0.1595137E+03_r8,     0.1585025E+03_r8,     0.1574978E+03_r8,&  
         0.1564994E+03_r8,     0.1555074E+03_r8,     0.1545216E+03_r8,     0.1535421E+03_r8,     0.1525688E+03_r8,&  
         0.1516017E+03_r8,     0.1506407E+03_r8,     0.1496858E+03_r8,     0.1487369E+03_r8,     0.1477941E+03_r8,&  
         0.1468572E+03_r8,     0.1459263E+03_r8,     0.1450013E+03_r8,     0.1440821E+03_r8,     0.1431688E+03_r8,&  
         0.1422612E+03_r8,     0.1413595E+03_r8,     0.1404634E+03_r8,     0.1395730E+03_r8,     0.1386882E+03_r8,&  
         0.1378091E+03_r8,     0.1369355E+03_r8,     0.1360675E+03_r8,     0.1352050E+03_r8,     0.1343479E+03_r8/)  

    psaditmk(1:150,   65)= (/ &
         0.2976893E+03_r8,     0.2968954E+03_r8,     0.2960984E+03_r8,     0.2952982E+03_r8,     0.2944946E+03_r8,&  
         0.2936875E+03_r8,     0.2928766E+03_r8,     0.2920619E+03_r8,     0.2912430E+03_r8,     0.2904198E+03_r8,&  
         0.2895921E+03_r8,     0.2887598E+03_r8,     0.2879224E+03_r8,     0.2870804E+03_r8,     0.2862325E+03_r8,&  
         0.2853795E+03_r8,     0.2845205E+03_r8,     0.2836555E+03_r8,     0.2827842E+03_r8,     0.2819062E+03_r8,&  
         0.2810214E+03_r8,     0.2801296E+03_r8,     0.2792303E+03_r8,     0.2783233E+03_r8,     0.2774082E+03_r8,&  
         0.2764851E+03_r8,     0.2755532E+03_r8,     0.2746123E+03_r8,     0.2736621E+03_r8,     0.2727025E+03_r8,&  
         0.2717328E+03_r8,     0.2707528E+03_r8,     0.2697622E+03_r8,     0.2687607E+03_r8,     0.2677479E+03_r8,&  
         0.2667234E+03_r8,     0.2656870E+03_r8,     0.2646384E+03_r8,     0.2635773E+03_r8,     0.2625032E+03_r8,&  
         0.2614162E+03_r8,     0.2603159E+03_r8,     0.2592021E+03_r8,     0.2580746E+03_r8,     0.2569333E+03_r8,&  
         0.2557782E+03_r8,     0.2546092E+03_r8,     0.2534264E+03_r8,     0.2522299E+03_r8,     0.2510198E+03_r8,&  
         0.2497963E+03_r8,     0.2485597E+03_r8,     0.2473104E+03_r8,     0.2460488E+03_r8,     0.2447754E+03_r8,&  
         0.2434907E+03_r8,     0.2421954E+03_r8,     0.2408902E+03_r8,     0.2395758E+03_r8,     0.2382529E+03_r8,&  
         0.2369225E+03_r8,     0.2355853E+03_r8,     0.2342422E+03_r8,     0.2328941E+03_r8,     0.2315418E+03_r8,&  
         0.2301862E+03_r8,     0.2288283E+03_r8,     0.2274688E+03_r8,     0.2261085E+03_r8,     0.2247482E+03_r8,&  
         0.2233887E+03_r8,     0.2220306E+03_r8,     0.2206747E+03_r8,     0.2193214E+03_r8,     0.2179714E+03_r8,&  
         0.2166251E+03_r8,     0.2152831E+03_r8,     0.2139458E+03_r8,     0.2126136E+03_r8,     0.2112868E+03_r8,&  
         0.2099657E+03_r8,     0.2086507E+03_r8,     0.2073418E+03_r8,     0.2060394E+03_r8,     0.2047437E+03_r8,&  
         0.2034548E+03_r8,     0.2021728E+03_r8,     0.2008978E+03_r8,     0.1996301E+03_r8,     0.1983695E+03_r8,&  
         0.1971162E+03_r8,     0.1958703E+03_r8,     0.1946317E+03_r8,     0.1934005E+03_r8,     0.1921768E+03_r8,&  
         0.1909604E+03_r8,     0.1897515E+03_r8,     0.1885500E+03_r8,     0.1873559E+03_r8,     0.1861692E+03_r8,&  
         0.1849899E+03_r8,     0.1838179E+03_r8,     0.1826532E+03_r8,     0.1814959E+03_r8,     0.1803458E+03_r8,&  
         0.1792029E+03_r8,     0.1780672E+03_r8,     0.1769386E+03_r8,     0.1758172E+03_r8,     0.1747029E+03_r8,&  
         0.1735956E+03_r8,     0.1724952E+03_r8,     0.1714019E+03_r8,     0.1703154E+03_r8,     0.1692359E+03_r8,&  
         0.1681631E+03_r8,     0.1670972E+03_r8,     0.1660380E+03_r8,     0.1649855E+03_r8,     0.1639397E+03_r8,&  
         0.1629005E+03_r8,     0.1618679E+03_r8,     0.1608418E+03_r8,     0.1598223E+03_r8,     0.1588092E+03_r8,&  
         0.1578025E+03_r8,     0.1568022E+03_r8,     0.1558082E+03_r8,     0.1548205E+03_r8,     0.1538391E+03_r8,&  
         0.1528640E+03_r8,     0.1518949E+03_r8,     0.1509321E+03_r8,     0.1499754E+03_r8,     0.1490247E+03_r8,&  
         0.1480800E+03_r8,     0.1471413E+03_r8,     0.1462086E+03_r8,     0.1452818E+03_r8,     0.1443608E+03_r8,&  
         0.1434457E+03_r8,     0.1425364E+03_r8,     0.1416329E+03_r8,     0.1407351E+03_r8,     0.1398430E+03_r8,&  
         0.1389565E+03_r8,     0.1380757E+03_r8,     0.1372004E+03_r8,     0.1363307E+03_r8,     0.1354665E+03_r8/)  

    psaditmk(1:150,   66)= (/ &
         0.2983355E+03_r8,     0.2975495E+03_r8,     0.2967610E+03_r8,     0.2959695E+03_r8,     0.2951747E+03_r8,&  
         0.2943766E+03_r8,     0.2935751E+03_r8,     0.2927698E+03_r8,     0.2919611E+03_r8,     0.2911481E+03_r8,&  
         0.2903310E+03_r8,     0.2895095E+03_r8,     0.2886833E+03_r8,     0.2878524E+03_r8,     0.2870168E+03_r8,&  
         0.2861757E+03_r8,     0.2853294E+03_r8,     0.2844773E+03_r8,     0.2836192E+03_r8,     0.2827552E+03_r8,&  
         0.2818845E+03_r8,     0.2810071E+03_r8,     0.2801230E+03_r8,     0.2792315E+03_r8,     0.2783324E+03_r8,&  
         0.2774257E+03_r8,     0.2765107E+03_r8,     0.2755872E+03_r8,     0.2746550E+03_r8,     0.2737137E+03_r8,&  
         0.2727630E+03_r8,     0.2718025E+03_r8,     0.2708316E+03_r8,     0.2698508E+03_r8,     0.2688589E+03_r8,&  
         0.2678559E+03_r8,     0.2668415E+03_r8,     0.2658152E+03_r8,     0.2647770E+03_r8,     0.2637263E+03_r8,&  
         0.2626629E+03_r8,     0.2615866E+03_r8,     0.2604970E+03_r8,     0.2593940E+03_r8,     0.2582774E+03_r8,&  
         0.2571471E+03_r8,     0.2560029E+03_r8,     0.2548449E+03_r8,     0.2536727E+03_r8,     0.2524869E+03_r8,&  
         0.2512872E+03_r8,     0.2500741E+03_r8,     0.2488475E+03_r8,     0.2476080E+03_r8,     0.2463560E+03_r8,&  
         0.2450917E+03_r8,     0.2438159E+03_r8,     0.2425289E+03_r8,     0.2412317E+03_r8,     0.2399246E+03_r8,&  
         0.2386088E+03_r8,     0.2372848E+03_r8,     0.2359535E+03_r8,     0.2346158E+03_r8,     0.2332726E+03_r8,&  
         0.2319247E+03_r8,     0.2305730E+03_r8,     0.2292185E+03_r8,     0.2278618E+03_r8,     0.2265040E+03_r8,&  
         0.2251457E+03_r8,     0.2237877E+03_r8,     0.2224308E+03_r8,     0.2210756E+03_r8,     0.2197228E+03_r8,&  
         0.2183729E+03_r8,     0.2170265E+03_r8,     0.2156841E+03_r8,     0.2143462E+03_r8,     0.2130131E+03_r8,&  
         0.2116853E+03_r8,     0.2103630E+03_r8,     0.2090466E+03_r8,     0.2077362E+03_r8,     0.2064323E+03_r8,&  
         0.2051348E+03_r8,     0.2038441E+03_r8,     0.2025602E+03_r8,     0.2012833E+03_r8,     0.2000135E+03_r8,&  
         0.1987510E+03_r8,     0.1974956E+03_r8,     0.1962475E+03_r8,     0.1950068E+03_r8,     0.1937735E+03_r8,&  
         0.1925475E+03_r8,     0.1913289E+03_r8,     0.1901178E+03_r8,     0.1889141E+03_r8,     0.1877178E+03_r8,&  
         0.1865289E+03_r8,     0.1853474E+03_r8,     0.1841732E+03_r8,     0.1830063E+03_r8,     0.1818468E+03_r8,&  
         0.1806944E+03_r8,     0.1795494E+03_r8,     0.1784115E+03_r8,     0.1772808E+03_r8,     0.1761573E+03_r8,&  
         0.1750408E+03_r8,     0.1739313E+03_r8,     0.1728289E+03_r8,     0.1717334E+03_r8,     0.1706449E+03_r8,&  
         0.1695632E+03_r8,     0.1684884E+03_r8,     0.1674204E+03_r8,     0.1663592E+03_r8,     0.1653046E+03_r8,&  
         0.1642568E+03_r8,     0.1632156E+03_r8,     0.1621810E+03_r8,     0.1611529E+03_r8,     0.1601314E+03_r8,&  
         0.1591163E+03_r8,     0.1581077E+03_r8,     0.1571055E+03_r8,     0.1561096E+03_r8,     0.1551200E+03_r8,&  
         0.1541367E+03_r8,     0.1531597E+03_r8,     0.1521888E+03_r8,     0.1512241E+03_r8,     0.1502655E+03_r8,&  
         0.1493129E+03_r8,     0.1483664E+03_r8,     0.1474259E+03_r8,     0.1464914E+03_r8,     0.1455628E+03_r8,&  
         0.1446401E+03_r8,     0.1437232E+03_r8,     0.1428122E+03_r8,     0.1419069E+03_r8,     0.1410073E+03_r8,&  
         0.1401135E+03_r8,     0.1392253E+03_r8,     0.1383428E+03_r8,     0.1374658E+03_r8,     0.1365944E+03_r8/)  

    psaditmk(1:150,   67)= (/ &
         0.2989657E+03_r8,     0.2981877E+03_r8,     0.2974071E+03_r8,     0.2966238E+03_r8,     0.2958374E+03_r8,&  
         0.2950481E+03_r8,     0.2942556E+03_r8,     0.2934595E+03_r8,     0.2926600E+03_r8,     0.2918569E+03_r8,&  
         0.2910499E+03_r8,     0.2902387E+03_r8,     0.2894233E+03_r8,     0.2886033E+03_r8,     0.2877786E+03_r8,&  
         0.2869494E+03_r8,     0.2861150E+03_r8,     0.2852750E+03_r8,     0.2844299E+03_r8,     0.2835788E+03_r8,&  
         0.2827219E+03_r8,     0.2818584E+03_r8,     0.2809884E+03_r8,     0.2801117E+03_r8,     0.2792279E+03_r8,&  
         0.2783366E+03_r8,     0.2774379E+03_r8,     0.2765310E+03_r8,     0.2756159E+03_r8,     0.2746924E+03_r8,&  
         0.2737597E+03_r8,     0.2728178E+03_r8,     0.2718663E+03_r8,     0.2709051E+03_r8,     0.2699334E+03_r8,&  
         0.2689511E+03_r8,     0.2679578E+03_r8,     0.2669534E+03_r8,     0.2659373E+03_r8,     0.2649093E+03_r8,&  
         0.2638688E+03_r8,     0.2628159E+03_r8,     0.2617501E+03_r8,     0.2606714E+03_r8,     0.2595792E+03_r8,&  
         0.2584735E+03_r8,     0.2573540E+03_r8,     0.2562207E+03_r8,     0.2550735E+03_r8,     0.2539123E+03_r8,&  
         0.2527372E+03_r8,     0.2515482E+03_r8,     0.2503454E+03_r8,     0.2491292E+03_r8,     0.2478996E+03_r8,&  
         0.2466572E+03_r8,     0.2454023E+03_r8,     0.2441355E+03_r8,     0.2428572E+03_r8,     0.2415680E+03_r8,&  
         0.2402689E+03_r8,     0.2389602E+03_r8,     0.2376430E+03_r8,     0.2363180E+03_r8,     0.2349860E+03_r8,&  
         0.2336479E+03_r8,     0.2323047E+03_r8,     0.2309572E+03_r8,     0.2296063E+03_r8,     0.2282528E+03_r8,&  
         0.2268976E+03_r8,     0.2255415E+03_r8,     0.2241854E+03_r8,     0.2228298E+03_r8,     0.2214756E+03_r8,&  
         0.2201234E+03_r8,     0.2187739E+03_r8,     0.2174275E+03_r8,     0.2160849E+03_r8,     0.2147464E+03_r8,&  
         0.2134126E+03_r8,     0.2120838E+03_r8,     0.2107604E+03_r8,     0.2094427E+03_r8,     0.2081310E+03_r8,&  
         0.2068254E+03_r8,     0.2055263E+03_r8,     0.2042339E+03_r8,     0.2029482E+03_r8,     0.2016693E+03_r8,&  
         0.2003976E+03_r8,     0.1991330E+03_r8,     0.1978755E+03_r8,     0.1966254E+03_r8,     0.1953825E+03_r8,&  
         0.1941470E+03_r8,     0.1929189E+03_r8,     0.1916982E+03_r8,     0.1904848E+03_r8,     0.1892789E+03_r8,&  
         0.1880804E+03_r8,     0.1868893E+03_r8,     0.1857055E+03_r8,     0.1845291E+03_r8,     0.1833601E+03_r8,&  
         0.1821983E+03_r8,     0.1810438E+03_r8,     0.1798966E+03_r8,     0.1787565E+03_r8,     0.1776237E+03_r8,&  
         0.1764979E+03_r8,     0.1753793E+03_r8,     0.1742677E+03_r8,     0.1731632E+03_r8,     0.1720656E+03_r8,&  
         0.1709750E+03_r8,     0.1698912E+03_r8,     0.1688143E+03_r8,     0.1677443E+03_r8,     0.1666810E+03_r8,&  
         0.1656244E+03_r8,     0.1645746E+03_r8,     0.1635313E+03_r8,     0.1624947E+03_r8,     0.1614647E+03_r8,&  
         0.1604412E+03_r8,     0.1594241E+03_r8,     0.1584136E+03_r8,     0.1574094E+03_r8,     0.1564116E+03_r8,&  
         0.1554201E+03_r8,     0.1544349E+03_r8,     0.1534559E+03_r8,     0.1524832E+03_r8,     0.1515166E+03_r8,&  
         0.1505561E+03_r8,     0.1496018E+03_r8,     0.1486534E+03_r8,     0.1477112E+03_r8,     0.1467748E+03_r8,&  
         0.1458444E+03_r8,     0.1449199E+03_r8,     0.1440013E+03_r8,     0.1430884E+03_r8,     0.1421814E+03_r8,&  
         0.1412801E+03_r8,     0.1403845E+03_r8,     0.1394946E+03_r8,     0.1386104E+03_r8,     0.1377318E+03_r8/)  

    psaditmk(1:150,   68)= (/ &
         0.2995810E+03_r8,     0.2988104E+03_r8,     0.2980375E+03_r8,     0.2972619E+03_r8,     0.2964836E+03_r8,&  
         0.2957026E+03_r8,     0.2949185E+03_r8,     0.2941314E+03_r8,     0.2933409E+03_r8,     0.2925470E+03_r8,&  
         0.2917495E+03_r8,     0.2909481E+03_r8,     0.2901428E+03_r8,     0.2893334E+03_r8,     0.2885196E+03_r8,&  
         0.2877014E+03_r8,     0.2868780E+03_r8,     0.2860502E+03_r8,     0.2852169E+03_r8,     0.2843785E+03_r8,&  
         0.2835342E+03_r8,     0.2826841E+03_r8,     0.2818278E+03_r8,     0.2809651E+03_r8,     0.2800959E+03_r8,&  
         0.2792196E+03_r8,     0.2783362E+03_r8,     0.2774452E+03_r8,     0.2765464E+03_r8,     0.2756396E+03_r8,&  
         0.2747243E+03_r8,     0.2738003E+03_r8,     0.2728671E+03_r8,     0.2719246E+03_r8,     0.2709724E+03_r8,&  
         0.2700100E+03_r8,     0.2690372E+03_r8,     0.2680537E+03_r8,     0.2670589E+03_r8,     0.2660527E+03_r8,&  
         0.2650348E+03_r8,     0.2640048E+03_r8,     0.2629622E+03_r8,     0.2619072E+03_r8,     0.2608391E+03_r8,&  
         0.2597577E+03_r8,     0.2586628E+03_r8,     0.2575541E+03_r8,     0.2564318E+03_r8,     0.2552955E+03_r8,&  
         0.2541451E+03_r8,     0.2529808E+03_r8,     0.2518024E+03_r8,     0.2506103E+03_r8,     0.2494044E+03_r8,&  
         0.2481850E+03_r8,     0.2469524E+03_r8,     0.2457071E+03_r8,     0.2444495E+03_r8,     0.2431801E+03_r8,&  
         0.2418995E+03_r8,     0.2406082E+03_r8,     0.2393071E+03_r8,     0.2379969E+03_r8,     0.2366784E+03_r8,&  
         0.2353524E+03_r8,     0.2340200E+03_r8,     0.2326817E+03_r8,     0.2313386E+03_r8,     0.2299916E+03_r8,&  
         0.2286415E+03_r8,     0.2272893E+03_r8,     0.2259357E+03_r8,     0.2245815E+03_r8,     0.2232275E+03_r8,&  
         0.2218746E+03_r8,     0.2205232E+03_r8,     0.2191741E+03_r8,     0.2178280E+03_r8,     0.2164852E+03_r8,&  
         0.2151464E+03_r8,     0.2138120E+03_r8,     0.2124823E+03_r8,     0.2111579E+03_r8,     0.2098390E+03_r8,&  
         0.2085259E+03_r8,     0.2072189E+03_r8,     0.2059182E+03_r8,     0.2046240E+03_r8,     0.2033365E+03_r8,&  
         0.2020559E+03_r8,     0.2007822E+03_r8,     0.1995155E+03_r8,     0.1982560E+03_r8,     0.1970038E+03_r8,&  
         0.1957588E+03_r8,     0.1945211E+03_r8,     0.1932909E+03_r8,     0.1920680E+03_r8,     0.1908525E+03_r8,&  
         0.1896444E+03_r8,     0.1884436E+03_r8,     0.1872503E+03_r8,     0.1860644E+03_r8,     0.1848857E+03_r8,&  
         0.1837145E+03_r8,     0.1825505E+03_r8,     0.1813938E+03_r8,     0.1802444E+03_r8,     0.1791022E+03_r8,&  
         0.1779671E+03_r8,     0.1768392E+03_r8,     0.1757185E+03_r8,     0.1746048E+03_r8,     0.1734981E+03_r8,&  
         0.1723984E+03_r8,     0.1713056E+03_r8,     0.1702198E+03_r8,     0.1691409E+03_r8,     0.1680687E+03_r8,&  
         0.1670034E+03_r8,     0.1659448E+03_r8,     0.1648929E+03_r8,     0.1638476E+03_r8,     0.1628091E+03_r8,&  
         0.1617770E+03_r8,     0.1607515E+03_r8,     0.1597325E+03_r8,     0.1587200E+03_r8,     0.1577139E+03_r8,&  
         0.1567142E+03_r8,     0.1557207E+03_r8,     0.1547336E+03_r8,     0.1537528E+03_r8,     0.1527781E+03_r8,&  
         0.1518097E+03_r8,     0.1508474E+03_r8,     0.1498912E+03_r8,     0.1489410E+03_r8,     0.1479969E+03_r8,&  
         0.1470587E+03_r8,     0.1461265E+03_r8,     0.1452002E+03_r8,     0.1442798E+03_r8,     0.1433652E+03_r8,&  
         0.1424564E+03_r8,     0.1415534E+03_r8,     0.1406561E+03_r8,     0.1397645E+03_r8,     0.1388785E+03_r8/)  

    psaditmk(1:150,   69)= (/ &
         0.3001818E+03_r8,     0.2994183E+03_r8,     0.2986526E+03_r8,     0.2978847E+03_r8,     0.2971142E+03_r8,&  
         0.2963410E+03_r8,     0.2955651E+03_r8,     0.2947862E+03_r8,     0.2940044E+03_r8,     0.2932194E+03_r8,&  
         0.2924310E+03_r8,     0.2916389E+03_r8,     0.2908434E+03_r8,     0.2900438E+03_r8,     0.2892401E+03_r8,&  
         0.2884325E+03_r8,     0.2876203E+03_r8,     0.2868036E+03_r8,     0.2859818E+03_r8,     0.2851549E+03_r8,&  
         0.2843232E+03_r8,     0.2834855E+03_r8,     0.2826422E+03_r8,     0.2817931E+03_r8,     0.2809374E+03_r8,&  
         0.2800756E+03_r8,     0.2792067E+03_r8,     0.2783310E+03_r8,     0.2774478E+03_r8,     0.2765568E+03_r8,&  
         0.2756582E+03_r8,     0.2747511E+03_r8,     0.2738355E+03_r8,     0.2729111E+03_r8,     0.2719774E+03_r8,&  
         0.2710341E+03_r8,     0.2700809E+03_r8,     0.2691176E+03_r8,     0.2681435E+03_r8,     0.2671586E+03_r8,&  
         0.2661623E+03_r8,     0.2651543E+03_r8,     0.2641345E+03_r8,     0.2631024E+03_r8,     0.2620577E+03_r8,&  
         0.2610002E+03_r8,     0.2599296E+03_r8,     0.2588453E+03_r8,     0.2577477E+03_r8,     0.2566361E+03_r8,&  
         0.2555107E+03_r8,     0.2543713E+03_r8,     0.2532178E+03_r8,     0.2520503E+03_r8,     0.2508687E+03_r8,&  
         0.2496733E+03_r8,     0.2484641E+03_r8,     0.2472417E+03_r8,     0.2460061E+03_r8,     0.2447579E+03_r8,&  
         0.2434977E+03_r8,     0.2422256E+03_r8,     0.2409427E+03_r8,     0.2396494E+03_r8,     0.2383465E+03_r8,&  
         0.2370349E+03_r8,     0.2357153E+03_r8,     0.2343885E+03_r8,     0.2330555E+03_r8,     0.2317171E+03_r8,&  
         0.2303744E+03_r8,     0.2290280E+03_r8,     0.2276790E+03_r8,     0.2263281E+03_r8,     0.2249761E+03_r8,&  
         0.2236241E+03_r8,     0.2222725E+03_r8,     0.2209221E+03_r8,     0.2195738E+03_r8,     0.2182279E+03_r8,&  
         0.2168852E+03_r8,     0.2155461E+03_r8,     0.2142112E+03_r8,     0.2128809E+03_r8,     0.2115555E+03_r8,&  
         0.2102355E+03_r8,     0.2089211E+03_r8,     0.2076127E+03_r8,     0.2063105E+03_r8,     0.2050146E+03_r8,&  
         0.2037254E+03_r8,     0.2024429E+03_r8,     0.2011673E+03_r8,     0.1998987E+03_r8,     0.1986372E+03_r8,&  
         0.1973829E+03_r8,     0.1961358E+03_r8,     0.1948961E+03_r8,     0.1936636E+03_r8,     0.1924385E+03_r8,&  
         0.1912208E+03_r8,     0.1900105E+03_r8,     0.1888076E+03_r8,     0.1876121E+03_r8,     0.1864239E+03_r8,&  
         0.1852431E+03_r8,     0.1840696E+03_r8,     0.1829035E+03_r8,     0.1817446E+03_r8,     0.1805930E+03_r8,&  
         0.1794486E+03_r8,     0.1783114E+03_r8,     0.1771813E+03_r8,     0.1760584E+03_r8,     0.1749425E+03_r8,&  
         0.1738337E+03_r8,     0.1727319E+03_r8,     0.1716371E+03_r8,     0.1705491E+03_r8,     0.1694681E+03_r8,&  
         0.1683939E+03_r8,     0.1673265E+03_r8,     0.1662659E+03_r8,     0.1652119E+03_r8,     0.1641647E+03_r8,&  
         0.1631241E+03_r8,     0.1620900E+03_r8,     0.1610626E+03_r8,     0.1600416E+03_r8,     0.1590271E+03_r8,&  
         0.1580191E+03_r8,     0.1570174E+03_r8,     0.1560220E+03_r8,     0.1550330E+03_r8,     0.1540503E+03_r8,&  
         0.1530737E+03_r8,     0.1521034E+03_r8,     0.1511393E+03_r8,     0.1501812E+03_r8,     0.1492292E+03_r8,&  
         0.1482832E+03_r8,     0.1473433E+03_r8,     0.1464093E+03_r8,     0.1454812E+03_r8,     0.1445590E+03_r8,&  
         0.1436426E+03_r8,     0.1427321E+03_r8,     0.1418273E+03_r8,     0.1409283E+03_r8,     0.1400349E+03_r8/)  

    psaditmk(1:150,   70)= (/ &
         0.3007687E+03_r8,     0.3000117E+03_r8,     0.2992534E+03_r8,     0.2984927E+03_r8,     0.2977292E+03_r8,&  
         0.2969637E+03_r8,     0.2961956E+03_r8,     0.2954248E+03_r8,     0.2946511E+03_r8,     0.2938745E+03_r8,&  
         0.2930948E+03_r8,     0.2923118E+03_r8,     0.2915255E+03_r8,     0.2907354E+03_r8,     0.2899415E+03_r8,&  
         0.2891436E+03_r8,     0.2883418E+03_r8,     0.2875357E+03_r8,     0.2867250E+03_r8,     0.2859096E+03_r8,&  
         0.2850892E+03_r8,     0.2842637E+03_r8,     0.2834328E+03_r8,     0.2825964E+03_r8,     0.2817540E+03_r8,&  
         0.2809055E+03_r8,     0.2800508E+03_r8,     0.2791892E+03_r8,     0.2783210E+03_r8,     0.2774454E+03_r8,&  
         0.2765625E+03_r8,     0.2756717E+03_r8,     0.2747728E+03_r8,     0.2738656E+03_r8,     0.2729496E+03_r8,&  
         0.2720246E+03_r8,     0.2710902E+03_r8,     0.2701460E+03_r8,     0.2691919E+03_r8,     0.2682273E+03_r8,&  
         0.2672519E+03_r8,     0.2662653E+03_r8,     0.2652675E+03_r8,     0.2642578E+03_r8,     0.2632358E+03_r8,&  
         0.2622015E+03_r8,     0.2611545E+03_r8,     0.2600944E+03_r8,     0.2590211E+03_r8,     0.2579342E+03_r8,&  
         0.2568336E+03_r8,     0.2557191E+03_r8,     0.2545906E+03_r8,     0.2534480E+03_r8,     0.2522913E+03_r8,&  
         0.2511204E+03_r8,     0.2499355E+03_r8,     0.2487369E+03_r8,     0.2475245E+03_r8,     0.2462990E+03_r8,&  
         0.2450605E+03_r8,     0.2438094E+03_r8,     0.2425464E+03_r8,     0.2412719E+03_r8,     0.2399868E+03_r8,&  
         0.2386915E+03_r8,     0.2373870E+03_r8,     0.2360740E+03_r8,     0.2347533E+03_r8,     0.2334259E+03_r8,&  
         0.2320925E+03_r8,     0.2307542E+03_r8,     0.2294119E+03_r8,     0.2280663E+03_r8,     0.2267184E+03_r8,&  
         0.2253690E+03_r8,     0.2240189E+03_r8,     0.2226690E+03_r8,     0.2213199E+03_r8,     0.2199724E+03_r8,&  
         0.2186270E+03_r8,     0.2172845E+03_r8,     0.2159454E+03_r8,     0.2146101E+03_r8,     0.2132792E+03_r8,&  
         0.2119530E+03_r8,     0.2106320E+03_r8,     0.2093164E+03_r8,     0.2080067E+03_r8,     0.2067029E+03_r8,&  
         0.2054055E+03_r8,     0.2041146E+03_r8,     0.2028303E+03_r8,     0.2015528E+03_r8,     0.2002823E+03_r8,&  
         0.1990188E+03_r8,     0.1977625E+03_r8,     0.1965134E+03_r8,     0.1952714E+03_r8,     0.1940369E+03_r8,&  
         0.1928097E+03_r8,     0.1915898E+03_r8,     0.1903773E+03_r8,     0.1891722E+03_r8,     0.1879744E+03_r8,&  
         0.1867840E+03_r8,     0.1856010E+03_r8,     0.1844254E+03_r8,     0.1832570E+03_r8,     0.1820959E+03_r8,&  
         0.1809421E+03_r8,     0.1797955E+03_r8,     0.1786562E+03_r8,     0.1775240E+03_r8,     0.1763989E+03_r8,&  
         0.1752809E+03_r8,     0.1741699E+03_r8,     0.1730660E+03_r8,     0.1719691E+03_r8,     0.1708790E+03_r8,&  
         0.1697959E+03_r8,     0.1687196E+03_r8,     0.1676501E+03_r8,     0.1665875E+03_r8,     0.1655315E+03_r8,&  
         0.1644822E+03_r8,     0.1634396E+03_r8,     0.1624036E+03_r8,     0.1613741E+03_r8,     0.1603512E+03_r8,&  
         0.1593347E+03_r8,     0.1583247E+03_r8,     0.1573211E+03_r8,     0.1563239E+03_r8,     0.1553329E+03_r8,&  
         0.1543483E+03_r8,     0.1533699E+03_r8,     0.1523977E+03_r8,     0.1514316E+03_r8,     0.1504717E+03_r8,&  
         0.1495179E+03_r8,     0.1485701E+03_r8,     0.1476283E+03_r8,     0.1466925E+03_r8,     0.1457626E+03_r8,&  
         0.1448386E+03_r8,     0.1439205E+03_r8,     0.1430082E+03_r8,     0.1421017E+03_r8,     0.1412009E+03_r8/)  

    psaditmk(1:150,   71)= (/ &
         0.3013422E+03_r8,     0.3005921E+03_r8,     0.2998403E+03_r8,     0.2990859E+03_r8,     0.2983299E+03_r8,&  
         0.2975716E+03_r8,     0.2968108E+03_r8,     0.2960478E+03_r8,     0.2952818E+03_r8,     0.2945133E+03_r8,&  
         0.2937419E+03_r8,     0.2929673E+03_r8,     0.2921896E+03_r8,     0.2914086E+03_r8,     0.2906241E+03_r8,&  
         0.2898362E+03_r8,     0.2890440E+03_r8,     0.2882477E+03_r8,     0.2874477E+03_r8,     0.2866430E+03_r8,&  
         0.2858337E+03_r8,     0.2850199E+03_r8,     0.2842006E+03_r8,     0.2833762E+03_r8,     0.2825464E+03_r8,&  
         0.2817108E+03_r8,     0.2808693E+03_r8,     0.2800217E+03_r8,     0.2791673E+03_r8,     0.2783065E+03_r8,&  
         0.2774384E+03_r8,     0.2765633E+03_r8,     0.2756803E+03_r8,     0.2747894E+03_r8,     0.2738904E+03_r8,&  
         0.2729828E+03_r8,     0.2720664E+03_r8,     0.2711408E+03_r8,     0.2702056E+03_r8,     0.2692605E+03_r8,&  
         0.2683052E+03_r8,     0.2673392E+03_r8,     0.2663626E+03_r8,     0.2653743E+03_r8,     0.2643745E+03_r8,&  
         0.2633629E+03_r8,     0.2623387E+03_r8,     0.2613023E+03_r8,     0.2602527E+03_r8,     0.2591902E+03_r8,&  
         0.2581141E+03_r8,     0.2570243E+03_r8,     0.2559206E+03_r8,     0.2548031E+03_r8,     0.2536714E+03_r8,&  
         0.2525256E+03_r8,     0.2513655E+03_r8,     0.2501913E+03_r8,     0.2490032E+03_r8,     0.2478012E+03_r8,&  
         0.2465857E+03_r8,     0.2453570E+03_r8,     0.2441154E+03_r8,     0.2428616E+03_r8,     0.2415959E+03_r8,&  
         0.2403190E+03_r8,     0.2390317E+03_r8,     0.2377345E+03_r8,     0.2364285E+03_r8,     0.2351142E+03_r8,&  
         0.2337926E+03_r8,     0.2324647E+03_r8,     0.2311312E+03_r8,     0.2297930E+03_r8,     0.2284512E+03_r8,&  
         0.2271066E+03_r8,     0.2257600E+03_r8,     0.2244122E+03_r8,     0.2230641E+03_r8,     0.2217165E+03_r8,&  
         0.2203700E+03_r8,     0.2190254E+03_r8,     0.2176833E+03_r8,     0.2163441E+03_r8,     0.2150087E+03_r8,&  
         0.2136773E+03_r8,     0.2123504E+03_r8,     0.2110285E+03_r8,     0.2097118E+03_r8,     0.2084008E+03_r8,&  
         0.2070957E+03_r8,     0.2057968E+03_r8,     0.2045041E+03_r8,     0.2032181E+03_r8,     0.2019388E+03_r8,&  
         0.2006664E+03_r8,     0.1994010E+03_r8,     0.1981426E+03_r8,     0.1968914E+03_r8,     0.1956474E+03_r8,&  
         0.1944107E+03_r8,     0.1931813E+03_r8,     0.1919593E+03_r8,     0.1907447E+03_r8,     0.1895373E+03_r8,&  
         0.1883374E+03_r8,     0.1871449E+03_r8,     0.1859596E+03_r8,     0.1847817E+03_r8,     0.1836112E+03_r8,&  
         0.1824479E+03_r8,     0.1812919E+03_r8,     0.1801432E+03_r8,     0.1790016E+03_r8,     0.1778672E+03_r8,&  
         0.1767400E+03_r8,     0.1756199E+03_r8,     0.1745068E+03_r8,     0.1734007E+03_r8,     0.1723017E+03_r8,&  
         0.1712095E+03_r8,     0.1701243E+03_r8,     0.1690460E+03_r8,     0.1679744E+03_r8,     0.1669097E+03_r8,&  
         0.1658517E+03_r8,     0.1648004E+03_r8,     0.1637558E+03_r8,     0.1627177E+03_r8,     0.1616863E+03_r8,&  
         0.1606614E+03_r8,     0.1596430E+03_r8,     0.1586310E+03_r8,     0.1576255E+03_r8,     0.1566263E+03_r8,&  
         0.1556334E+03_r8,     0.1546469E+03_r8,     0.1536665E+03_r8,     0.1526925E+03_r8,     0.1517245E+03_r8,&  
         0.1507628E+03_r8,     0.1498071E+03_r8,     0.1488575E+03_r8,     0.1479139E+03_r8,     0.1469763E+03_r8,&  
         0.1460446E+03_r8,     0.1451188E+03_r8,     0.1441989E+03_r8,     0.1432848E+03_r8,     0.1423765E+03_r8/)  

    psaditmk(1:150,   72)= (/ &
         0.3019030E+03_r8,     0.3011588E+03_r8,     0.3004132E+03_r8,     0.2996656E+03_r8,     0.2989163E+03_r8,&  
         0.2981650E+03_r8,     0.2974112E+03_r8,     0.2966555E+03_r8,     0.2958972E+03_r8,     0.2951363E+03_r8,&  
         0.2943727E+03_r8,     0.2936064E+03_r8,     0.2928372E+03_r8,     0.2920648E+03_r8,     0.2912890E+03_r8,&  
         0.2905099E+03_r8,     0.2897273E+03_r8,     0.2889410E+03_r8,     0.2881508E+03_r8,     0.2873563E+03_r8,&  
         0.2865575E+03_r8,     0.2857542E+03_r8,     0.2849465E+03_r8,     0.2841336E+03_r8,     0.2833159E+03_r8,&  
         0.2824926E+03_r8,     0.2816637E+03_r8,     0.2808289E+03_r8,     0.2799883E+03_r8,     0.2791411E+03_r8,&  
         0.2782875E+03_r8,     0.2774269E+03_r8,     0.2765594E+03_r8,     0.2756841E+03_r8,     0.2748013E+03_r8,&  
         0.2739103E+03_r8,     0.2730110E+03_r8,     0.2721030E+03_r8,     0.2711859E+03_r8,     0.2702596E+03_r8,&  
         0.2693234E+03_r8,     0.2683773E+03_r8,     0.2674209E+03_r8,     0.2664534E+03_r8,     0.2654749E+03_r8,&  
         0.2644853E+03_r8,     0.2634836E+03_r8,     0.2624698E+03_r8,     0.2614436E+03_r8,     0.2604046E+03_r8,&  
         0.2593525E+03_r8,     0.2582872E+03_r8,     0.2572083E+03_r8,     0.2561156E+03_r8,     0.2550089E+03_r8,&  
         0.2538880E+03_r8,     0.2527531E+03_r8,     0.2516039E+03_r8,     0.2504405E+03_r8,     0.2492630E+03_r8,&  
         0.2480715E+03_r8,     0.2468662E+03_r8,     0.2456475E+03_r8,     0.2444156E+03_r8,     0.2431711E+03_r8,&  
         0.2419144E+03_r8,     0.2406462E+03_r8,     0.2393670E+03_r8,     0.2380776E+03_r8,     0.2367786E+03_r8,&  
         0.2354710E+03_r8,     0.2341557E+03_r8,     0.2328334E+03_r8,     0.2315050E+03_r8,     0.2301714E+03_r8,&  
         0.2288336E+03_r8,     0.2274925E+03_r8,     0.2261489E+03_r8,     0.2248037E+03_r8,     0.2234578E+03_r8,&  
         0.2221118E+03_r8,     0.2207666E+03_r8,     0.2194229E+03_r8,     0.2180813E+03_r8,     0.2167424E+03_r8,&  
         0.2154068E+03_r8,     0.2140751E+03_r8,     0.2127477E+03_r8,     0.2114249E+03_r8,     0.2101073E+03_r8,&  
         0.2087951E+03_r8,     0.2074887E+03_r8,     0.2061882E+03_r8,     0.2048940E+03_r8,     0.2036064E+03_r8,&  
         0.2023252E+03_r8,     0.2010510E+03_r8,     0.1997836E+03_r8,     0.1985232E+03_r8,     0.1972700E+03_r8,&  
         0.1960240E+03_r8,     0.1947852E+03_r8,     0.1935537E+03_r8,     0.1923295E+03_r8,     0.1911127E+03_r8,&  
         0.1899032E+03_r8,     0.1887010E+03_r8,     0.1875063E+03_r8,     0.1863189E+03_r8,     0.1851388E+03_r8,&  
         0.1839660E+03_r8,     0.1828006E+03_r8,     0.1816424E+03_r8,     0.1804914E+03_r8,     0.1793477E+03_r8,&  
         0.1782112E+03_r8,     0.1770818E+03_r8,     0.1759595E+03_r8,     0.1748443E+03_r8,     0.1737361E+03_r8,&  
         0.1726349E+03_r8,     0.1715407E+03_r8,     0.1704534E+03_r8,     0.1693729E+03_r8,     0.1682994E+03_r8,&  
         0.1672326E+03_r8,     0.1661725E+03_r8,     0.1651192E+03_r8,     0.1640725E+03_r8,     0.1630325E+03_r8,&  
         0.1619991E+03_r8,     0.1609722E+03_r8,     0.1599518E+03_r8,     0.1589379E+03_r8,     0.1579303E+03_r8,&  
         0.1569292E+03_r8,     0.1559345E+03_r8,     0.1549460E+03_r8,     0.1539638E+03_r8,     0.1529878E+03_r8,&  
         0.1520181E+03_r8,     0.1510544E+03_r8,     0.1500969E+03_r8,     0.1491454E+03_r8,     0.1482000E+03_r8,&  
         0.1472606E+03_r8,     0.1463271E+03_r8,     0.1453995E+03_r8,     0.1444778E+03_r8,     0.1435620E+03_r8/)  

    psaditmk(1:150,   73)= (/ &
         0.3024511E+03_r8,     0.3017130E+03_r8,     0.3009735E+03_r8,     0.3002324E+03_r8,     0.2994895E+03_r8,&  
         0.2987444E+03_r8,     0.2979978E+03_r8,     0.2972487E+03_r8,     0.2964976E+03_r8,     0.2957441E+03_r8,&  
         0.2949882E+03_r8,     0.2942296E+03_r8,     0.2934682E+03_r8,     0.2927041E+03_r8,     0.2919370E+03_r8,&  
         0.2911668E+03_r8,     0.2903929E+03_r8,     0.2896155E+03_r8,     0.2888346E+03_r8,     0.2880502E+03_r8,&  
         0.2872617E+03_r8,     0.2864688E+03_r8,     0.2856715E+03_r8,     0.2848695E+03_r8,     0.2840631E+03_r8,&  
         0.2832517E+03_r8,     0.2824348E+03_r8,     0.2816125E+03_r8,     0.2807845E+03_r8,     0.2799507E+03_r8,&  
         0.2791105E+03_r8,     0.2782642E+03_r8,     0.2774108E+03_r8,     0.2765508E+03_r8,     0.2756831E+03_r8,&  
         0.2748081E+03_r8,     0.2739251E+03_r8,     0.2730340E+03_r8,     0.2721345E+03_r8,     0.2712259E+03_r8,&  
         0.2703080E+03_r8,     0.2693809E+03_r8,     0.2684438E+03_r8,     0.2674964E+03_r8,     0.2665386E+03_r8,&  
         0.2655699E+03_r8,     0.2645898E+03_r8,     0.2635981E+03_r8,     0.2625944E+03_r8,     0.2615785E+03_r8,&  
         0.2605498E+03_r8,     0.2595084E+03_r8,     0.2584537E+03_r8,     0.2573855E+03_r8,     0.2563036E+03_r8,&  
         0.2552078E+03_r8,     0.2540979E+03_r8,     0.2529739E+03_r8,     0.2518356E+03_r8,     0.2506829E+03_r8,&  
         0.2495161E+03_r8,     0.2483353E+03_r8,     0.2471404E+03_r8,     0.2459318E+03_r8,     0.2447099E+03_r8,&  
         0.2434748E+03_r8,     0.2422275E+03_r8,     0.2409680E+03_r8,     0.2396973E+03_r8,     0.2384158E+03_r8,&  
         0.2371244E+03_r8,     0.2358238E+03_r8,     0.2345148E+03_r8,     0.2331984E+03_r8,     0.2318754E+03_r8,&  
         0.2305468E+03_r8,     0.2292134E+03_r8,     0.2278761E+03_r8,     0.2265358E+03_r8,     0.2251934E+03_r8,&  
         0.2238498E+03_r8,     0.2225057E+03_r8,     0.2211620E+03_r8,     0.2198194E+03_r8,     0.2184785E+03_r8,&  
         0.2171400E+03_r8,     0.2158046E+03_r8,     0.2144726E+03_r8,     0.2131447E+03_r8,     0.2118213E+03_r8,&  
         0.2105027E+03_r8,     0.2091895E+03_r8,     0.2078818E+03_r8,     0.2065800E+03_r8,     0.2052843E+03_r8,&  
         0.2039949E+03_r8,     0.2027121E+03_r8,     0.2014360E+03_r8,     0.2001667E+03_r8,     0.1989045E+03_r8,&  
         0.1976492E+03_r8,     0.1964011E+03_r8,     0.1951602E+03_r8,     0.1939266E+03_r8,     0.1927003E+03_r8,&  
         0.1914813E+03_r8,     0.1902696E+03_r8,     0.1890653E+03_r8,     0.1878683E+03_r8,     0.1866788E+03_r8,&  
         0.1854965E+03_r8,     0.1843215E+03_r8,     0.1831539E+03_r8,     0.1819935E+03_r8,     0.1808404E+03_r8,&  
         0.1796945E+03_r8,     0.1785558E+03_r8,     0.1774242E+03_r8,     0.1762998E+03_r8,     0.1751824E+03_r8,&  
         0.1740721E+03_r8,     0.1729688E+03_r8,     0.1718725E+03_r8,     0.1707831E+03_r8,     0.1697006E+03_r8,&  
         0.1686249E+03_r8,     0.1675560E+03_r8,     0.1664940E+03_r8,     0.1654386E+03_r8,     0.1643899E+03_r8,&  
         0.1633479E+03_r8,     0.1623124E+03_r8,     0.1612835E+03_r8,     0.1602612E+03_r8,     0.1592453E+03_r8,&  
         0.1582358E+03_r8,     0.1572328E+03_r8,     0.1562361E+03_r8,     0.1552457E+03_r8,     0.1542616E+03_r8,&  
         0.1532838E+03_r8,     0.1523121E+03_r8,     0.1513466E+03_r8,     0.1503873E+03_r8,     0.1494339E+03_r8,&  
         0.1484867E+03_r8,     0.1475454E+03_r8,     0.1466101E+03_r8,     0.1456808E+03_r8,     0.1447573E+03_r8/)  

    psaditmk(1:150,   74)= (/ &
         0.3029872E+03_r8,     0.3022549E+03_r8,     0.3015215E+03_r8,     0.3007861E+03_r8,     0.3000493E+03_r8,&  
         0.2993106E+03_r8,     0.2985704E+03_r8,     0.2978282E+03_r8,     0.2970839E+03_r8,     0.2963376E+03_r8,&  
         0.2955887E+03_r8,     0.2948374E+03_r8,     0.2940839E+03_r8,     0.2933277E+03_r8,     0.2925685E+03_r8,&  
         0.2918063E+03_r8,     0.2910412E+03_r8,     0.2902728E+03_r8,     0.2895012E+03_r8,     0.2887255E+03_r8,&  
         0.2879466E+03_r8,     0.2871637E+03_r8,     0.2863767E+03_r8,     0.2855853E+03_r8,     0.2847893E+03_r8,&  
         0.2839891E+03_r8,     0.2831837E+03_r8,     0.2823731E+03_r8,     0.2815575E+03_r8,     0.2807362E+03_r8,&  
         0.2799091E+03_r8,     0.2790758E+03_r8,     0.2782365E+03_r8,     0.2773903E+03_r8,     0.2765376E+03_r8,&  
         0.2756775E+03_r8,     0.2748102E+03_r8,     0.2739351E+03_r8,     0.2730520E+03_r8,     0.2721605E+03_r8,&  
         0.2712603E+03_r8,     0.2703512E+03_r8,     0.2694328E+03_r8,     0.2685046E+03_r8,     0.2675663E+03_r8,&  
         0.2666179E+03_r8,     0.2656584E+03_r8,     0.2646880E+03_r8,     0.2637062E+03_r8,     0.2627127E+03_r8,&  
         0.2617070E+03_r8,     0.2606888E+03_r8,     0.2596577E+03_r8,     0.2586136E+03_r8,     0.2575561E+03_r8,&  
         0.2564850E+03_r8,     0.2554000E+03_r8,     0.2543010E+03_r8,     0.2531880E+03_r8,     0.2520606E+03_r8,&  
         0.2509189E+03_r8,     0.2497629E+03_r8,     0.2485926E+03_r8,     0.2474082E+03_r8,     0.2462099E+03_r8,&  
         0.2449980E+03_r8,     0.2437728E+03_r8,     0.2425347E+03_r8,     0.2412845E+03_r8,     0.2400223E+03_r8,&  
         0.2387491E+03_r8,     0.2374654E+03_r8,     0.2361722E+03_r8,     0.2348700E+03_r8,     0.2335598E+03_r8,&  
         0.2322425E+03_r8,     0.2309190E+03_r8,     0.2295902E+03_r8,     0.2282570E+03_r8,     0.2269203E+03_r8,&  
         0.2255810E+03_r8,     0.2242400E+03_r8,     0.2228981E+03_r8,     0.2215561E+03_r8,     0.2202148E+03_r8,&  
         0.2188749E+03_r8,     0.2175370E+03_r8,     0.2162017E+03_r8,     0.2148697E+03_r8,     0.2135415E+03_r8,&  
         0.2122175E+03_r8,     0.2108982E+03_r8,     0.2095839E+03_r8,     0.2082751E+03_r8,     0.2069720E+03_r8,&  
         0.2056748E+03_r8,     0.2043839E+03_r8,     0.2030994E+03_r8,     0.2018215E+03_r8,     0.2005504E+03_r8,&  
         0.1992862E+03_r8,     0.1980289E+03_r8,     0.1967788E+03_r8,     0.1955359E+03_r8,     0.1943002E+03_r8,&  
         0.1930717E+03_r8,     0.1918506E+03_r8,     0.1906367E+03_r8,     0.1894302E+03_r8,     0.1882311E+03_r8,&  
         0.1870393E+03_r8,     0.1858548E+03_r8,     0.1846777E+03_r8,     0.1835079E+03_r8,     0.1823453E+03_r8,&  
         0.1811900E+03_r8,     0.1800419E+03_r8,     0.1789010E+03_r8,     0.1777673E+03_r8,     0.1766407E+03_r8,&  
         0.1755212E+03_r8,     0.1744088E+03_r8,     0.1733034E+03_r8,     0.1722049E+03_r8,     0.1711134E+03_r8,&  
         0.1700288E+03_r8,     0.1689511E+03_r8,     0.1678802E+03_r8,     0.1668160E+03_r8,     0.1657586E+03_r8,&  
         0.1647079E+03_r8,     0.1636639E+03_r8,     0.1626264E+03_r8,     0.1615956E+03_r8,     0.1605712E+03_r8,&  
         0.1595533E+03_r8,     0.1585420E+03_r8,     0.1575370E+03_r8,     0.1565384E+03_r8,     0.1555461E+03_r8,&  
         0.1545601E+03_r8,     0.1535803E+03_r8,     0.1526068E+03_r8,     0.1516394E+03_r8,     0.1506782E+03_r8,&  
         0.1497230E+03_r8,     0.1487739E+03_r8,     0.1478309E+03_r8,     0.1468938E+03_r8,     0.1459626E+03_r8/)  

    psaditmk(1:150,   75)= (/ &
         0.3035119E+03_r8,     0.3027850E+03_r8,     0.3020572E+03_r8,     0.3013277E+03_r8,     0.3005967E+03_r8,&  
         0.2998644E+03_r8,     0.2991304E+03_r8,     0.2983941E+03_r8,     0.2976567E+03_r8,     0.2969169E+03_r8,&  
         0.2961749E+03_r8,     0.2954309E+03_r8,     0.2946847E+03_r8,     0.2939357E+03_r8,     0.2931842E+03_r8,&  
         0.2924301E+03_r8,     0.2916733E+03_r8,     0.2909130E+03_r8,     0.2901498E+03_r8,     0.2893834E+03_r8,&  
         0.2886138E+03_r8,     0.2878398E+03_r8,     0.2870625E+03_r8,     0.2862812E+03_r8,     0.2854957E+03_r8,&  
         0.2847058E+03_r8,     0.2839114E+03_r8,     0.2831122E+03_r8,     0.2823080E+03_r8,     0.2814987E+03_r8,&  
         0.2806839E+03_r8,     0.2798634E+03_r8,     0.2790370E+03_r8,     0.2782047E+03_r8,     0.2773655E+03_r8,&  
         0.2765201E+03_r8,     0.2756674E+03_r8,     0.2748077E+03_r8,     0.2739402E+03_r8,     0.2730651E+03_r8,&  
         0.2721816E+03_r8,     0.2712898E+03_r8,     0.2703891E+03_r8,     0.2694793E+03_r8,     0.2685599E+03_r8,&  
         0.2676308E+03_r8,     0.2666913E+03_r8,     0.2657414E+03_r8,     0.2647807E+03_r8,     0.2638086E+03_r8,&  
         0.2628248E+03_r8,     0.2618290E+03_r8,     0.2608210E+03_r8,     0.2598004E+03_r8,     0.2587668E+03_r8,&  
         0.2577200E+03_r8,     0.2566597E+03_r8,     0.2555856E+03_r8,     0.2544976E+03_r8,     0.2533954E+03_r8,&  
         0.2522789E+03_r8,     0.2511480E+03_r8,     0.2500029E+03_r8,     0.2488433E+03_r8,     0.2476696E+03_r8,&  
         0.2464817E+03_r8,     0.2452800E+03_r8,     0.2440648E+03_r8,     0.2428364E+03_r8,     0.2415953E+03_r8,&  
         0.2403422E+03_r8,     0.2390774E+03_r8,     0.2378018E+03_r8,     0.2365159E+03_r8,     0.2352208E+03_r8,&  
         0.2339172E+03_r8,     0.2326060E+03_r8,     0.2312879E+03_r8,     0.2299641E+03_r8,     0.2286353E+03_r8,&  
         0.2273024E+03_r8,     0.2259665E+03_r8,     0.2246284E+03_r8,     0.2232889E+03_r8,     0.2219489E+03_r8,&  
         0.2206090E+03_r8,     0.2192702E+03_r8,     0.2179331E+03_r8,     0.2165983E+03_r8,     0.2152664E+03_r8,&  
         0.2139380E+03_r8,     0.2126136E+03_r8,     0.2112935E+03_r8,     0.2099784E+03_r8,     0.2086685E+03_r8,&  
         0.2073641E+03_r8,     0.2060656E+03_r8,     0.2047731E+03_r8,     0.2034870E+03_r8,     0.2022074E+03_r8,&  
         0.2009345E+03_r8,     0.1996683E+03_r8,     0.1984092E+03_r8,     0.1971571E+03_r8,     0.1959121E+03_r8,&  
         0.1946743E+03_r8,     0.1934437E+03_r8,     0.1922204E+03_r8,     0.1910045E+03_r8,     0.1897958E+03_r8,&  
         0.1885945E+03_r8,     0.1874005E+03_r8,     0.1862139E+03_r8,     0.1850345E+03_r8,     0.1838625E+03_r8,&  
         0.1826977E+03_r8,     0.1815402E+03_r8,     0.1803900E+03_r8,     0.1792469E+03_r8,     0.1781110E+03_r8,&  
         0.1769823E+03_r8,     0.1758607E+03_r8,     0.1747461E+03_r8,     0.1736385E+03_r8,     0.1725380E+03_r8,&  
         0.1714444E+03_r8,     0.1703577E+03_r8,     0.1692779E+03_r8,     0.1682049E+03_r8,     0.1671387E+03_r8,&  
         0.1660793E+03_r8,     0.1650265E+03_r8,     0.1639804E+03_r8,     0.1629410E+03_r8,     0.1619081E+03_r8,&  
         0.1608818E+03_r8,     0.1598620E+03_r8,     0.1588486E+03_r8,     0.1578417E+03_r8,     0.1568412E+03_r8,&  
         0.1558470E+03_r8,     0.1548591E+03_r8,     0.1538774E+03_r8,     0.1529020E+03_r8,     0.1519327E+03_r8,&  
         0.1509696E+03_r8,     0.1500126E+03_r8,     0.1490617E+03_r8,     0.1481168E+03_r8,     0.1471779E+03_r8/)  

    psaditmk(1:150,   76)= (/ &
         0.3040254E+03_r8,     0.3033041E+03_r8,     0.3025814E+03_r8,     0.3018573E+03_r8,     0.3011320E+03_r8,&  
         0.3004055E+03_r8,     0.2996770E+03_r8,     0.2989473E+03_r8,     0.2982160E+03_r8,     0.2974827E+03_r8,&  
         0.2967473E+03_r8,     0.2960102E+03_r8,     0.2952708E+03_r8,     0.2945291E+03_r8,     0.2937851E+03_r8,&  
         0.2930385E+03_r8,     0.2922891E+03_r8,     0.2915372E+03_r8,     0.2907822E+03_r8,     0.2900242E+03_r8,&  
         0.2892627E+03_r8,     0.2884984E+03_r8,     0.2877302E+03_r8,     0.2869583E+03_r8,     0.2861826E+03_r8,&  
         0.2854026E+03_r8,     0.2846188E+03_r8,     0.2838302E+03_r8,     0.2830371E+03_r8,     0.2822392E+03_r8,&  
         0.2814361E+03_r8,     0.2806277E+03_r8,     0.2798139E+03_r8,     0.2789942E+03_r8,     0.2781687E+03_r8,&  
         0.2773365E+03_r8,     0.2764982E+03_r8,     0.2756527E+03_r8,     0.2748006E+03_r8,     0.2739408E+03_r8,&  
         0.2730734E+03_r8,     0.2721979E+03_r8,     0.2713142E+03_r8,     0.2704218E+03_r8,     0.2695204E+03_r8,&  
         0.2686099E+03_r8,     0.2676895E+03_r8,     0.2667592E+03_r8,     0.2658184E+03_r8,     0.2648671E+03_r8,&  
         0.2639044E+03_r8,     0.2629306E+03_r8,     0.2619449E+03_r8,     0.2609471E+03_r8,     0.2599368E+03_r8,&  
         0.2589135E+03_r8,     0.2578774E+03_r8,     0.2568277E+03_r8,     0.2557644E+03_r8,     0.2546872E+03_r8,&  
         0.2535959E+03_r8,     0.2524904E+03_r8,     0.2513705E+03_r8,     0.2502362E+03_r8,     0.2490875E+03_r8,&  
         0.2479244E+03_r8,     0.2467471E+03_r8,     0.2455557E+03_r8,     0.2443506E+03_r8,     0.2431320E+03_r8,&  
         0.2419005E+03_r8,     0.2406564E+03_r8,     0.2394004E+03_r8,     0.2381331E+03_r8,     0.2368552E+03_r8,&  
         0.2355674E+03_r8,     0.2342706E+03_r8,     0.2329656E+03_r8,     0.2316534E+03_r8,     0.2303347E+03_r8,&  
         0.2290106E+03_r8,     0.2276820E+03_r8,     0.2263497E+03_r8,     0.2250147E+03_r8,     0.2236779E+03_r8,&  
         0.2223400E+03_r8,     0.2210019E+03_r8,     0.2196645E+03_r8,     0.2183283E+03_r8,     0.2169940E+03_r8,&  
         0.2156624E+03_r8,     0.2143340E+03_r8,     0.2130092E+03_r8,     0.2116887E+03_r8,     0.2103728E+03_r8,&  
         0.2090619E+03_r8,     0.2077564E+03_r8,     0.2064565E+03_r8,     0.2051626E+03_r8,     0.2038750E+03_r8,&  
         0.2025937E+03_r8,     0.2013190E+03_r8,     0.2000510E+03_r8,     0.1987900E+03_r8,     0.1975359E+03_r8,&  
         0.1962888E+03_r8,     0.1950490E+03_r8,     0.1938163E+03_r8,     0.1925909E+03_r8,     0.1913728E+03_r8,&  
         0.1901620E+03_r8,     0.1889585E+03_r8,     0.1877623E+03_r8,     0.1865735E+03_r8,     0.1853920E+03_r8,&  
         0.1842178E+03_r8,     0.1830508E+03_r8,     0.1818911E+03_r8,     0.1807387E+03_r8,     0.1795935E+03_r8,&  
         0.1784554E+03_r8,     0.1773245E+03_r8,     0.1762007E+03_r8,     0.1750840E+03_r8,     0.1739744E+03_r8,&  
         0.1728717E+03_r8,     0.1717760E+03_r8,     0.1706872E+03_r8,     0.1696053E+03_r8,     0.1685302E+03_r8,&  
         0.1674620E+03_r8,     0.1664005E+03_r8,     0.1653457E+03_r8,     0.1642976E+03_r8,     0.1632562E+03_r8,&  
         0.1622213E+03_r8,     0.1611930E+03_r8,     0.1601712E+03_r8,     0.1591559E+03_r8,     0.1581470E+03_r8,&  
         0.1571445E+03_r8,     0.1561484E+03_r8,     0.1551586E+03_r8,     0.1541751E+03_r8,     0.1531978E+03_r8,&  
         0.1522266E+03_r8,     0.1512617E+03_r8,     0.1503028E+03_r8,     0.1493501E+03_r8,     0.1484034E+03_r8/)  

    psaditmk(1:150,   77)= (/ &
         0.3045285E+03_r8,     0.3038120E+03_r8,     0.3030945E+03_r8,     0.3023757E+03_r8,     0.3016561E+03_r8,&  
         0.3009348E+03_r8,     0.3002123E+03_r8,     0.2994883E+03_r8,     0.2987628E+03_r8,     0.2980355E+03_r8,&  
         0.2973069E+03_r8,     0.2965761E+03_r8,     0.2958431E+03_r8,     0.2951085E+03_r8,     0.2943715E+03_r8,&  
         0.2936319E+03_r8,     0.2928901E+03_r8,     0.2921460E+03_r8,     0.2913986E+03_r8,     0.2906486E+03_r8,&  
         0.2898958E+03_r8,     0.2891397E+03_r8,     0.2883803E+03_r8,     0.2876176E+03_r8,     0.2868512E+03_r8,&  
         0.2860809E+03_r8,     0.2853068E+03_r8,     0.2845284E+03_r8,     0.2837458E+03_r8,     0.2829586E+03_r8,&  
         0.2821667E+03_r8,     0.2813699E+03_r8,     0.2805678E+03_r8,     0.2797603E+03_r8,     0.2789474E+03_r8,&  
         0.2781286E+03_r8,     0.2773033E+03_r8,     0.2764721E+03_r8,     0.2756337E+03_r8,     0.2747889E+03_r8,&  
         0.2739367E+03_r8,     0.2730769E+03_r8,     0.2722094E+03_r8,     0.2713337E+03_r8,     0.2704494E+03_r8,&  
         0.2695565E+03_r8,     0.2686544E+03_r8,     0.2677427E+03_r8,     0.2668214E+03_r8,     0.2658898E+03_r8,&  
         0.2649478E+03_r8,     0.2639948E+03_r8,     0.2630305E+03_r8,     0.2620546E+03_r8,     0.2610668E+03_r8,&  
         0.2600666E+03_r8,     0.2590539E+03_r8,     0.2580282E+03_r8,     0.2569891E+03_r8,     0.2559365E+03_r8,&  
         0.2548702E+03_r8,     0.2537898E+03_r8,     0.2526953E+03_r8,     0.2515863E+03_r8,     0.2504629E+03_r8,&  
         0.2493250E+03_r8,     0.2481727E+03_r8,     0.2470061E+03_r8,     0.2458251E+03_r8,     0.2446302E+03_r8,&  
         0.2434217E+03_r8,     0.2421998E+03_r8,     0.2409651E+03_r8,     0.2397180E+03_r8,     0.2384593E+03_r8,&  
         0.2371895E+03_r8,     0.2359094E+03_r8,     0.2346198E+03_r8,     0.2333214E+03_r8,     0.2320152E+03_r8,&  
         0.2307021E+03_r8,     0.2293830E+03_r8,     0.2280588E+03_r8,     0.2267304E+03_r8,     0.2253988E+03_r8,&  
         0.2240649E+03_r8,     0.2227295E+03_r8,     0.2213934E+03_r8,     0.2200575E+03_r8,     0.2187225E+03_r8,&  
         0.2173891E+03_r8,     0.2160579E+03_r8,     0.2147296E+03_r8,     0.2134047E+03_r8,     0.2120837E+03_r8,&  
         0.2107671E+03_r8,     0.2094554E+03_r8,     0.2081488E+03_r8,     0.2068477E+03_r8,     0.2055524E+03_r8,&  
         0.2042633E+03_r8,     0.2029803E+03_r8,     0.2017039E+03_r8,     0.2004342E+03_r8,     0.1991712E+03_r8,&  
         0.1979152E+03_r8,     0.1966662E+03_r8,     0.1954243E+03_r8,     0.1941895E+03_r8,     0.1929620E+03_r8,&  
         0.1917418E+03_r8,     0.1905288E+03_r8,     0.1893232E+03_r8,     0.1881248E+03_r8,     0.1869338E+03_r8,&  
         0.1857501E+03_r8,     0.1845737E+03_r8,     0.1834045E+03_r8,     0.1822427E+03_r8,     0.1810881E+03_r8,&  
         0.1799407E+03_r8,     0.1788005E+03_r8,     0.1776674E+03_r8,     0.1765415E+03_r8,     0.1754226E+03_r8,&  
         0.1743108E+03_r8,     0.1732060E+03_r8,     0.1721082E+03_r8,     0.1710173E+03_r8,     0.1699333E+03_r8,&  
         0.1688562E+03_r8,     0.1677859E+03_r8,     0.1667224E+03_r8,     0.1656656E+03_r8,     0.1646154E+03_r8,&  
         0.1635720E+03_r8,     0.1625351E+03_r8,     0.1615048E+03_r8,     0.1604811E+03_r8,     0.1594638E+03_r8,&  
         0.1584530E+03_r8,     0.1574485E+03_r8,     0.1564505E+03_r8,     0.1554588E+03_r8,     0.1544733E+03_r8,&  
         0.1534941E+03_r8,     0.1525211E+03_r8,     0.1515543E+03_r8,     0.1505936E+03_r8,     0.1496390E+03_r8/)  

    psaditmk(1:150,   78)= (/ &
         0.3050211E+03_r8,     0.3043095E+03_r8,     0.3035968E+03_r8,     0.3028830E+03_r8,     0.3021688E+03_r8,&  
         0.3014527E+03_r8,     0.3007356E+03_r8,     0.3000173E+03_r8,     0.2992975E+03_r8,     0.2985759E+03_r8,&  
         0.2978533E+03_r8,     0.2971288E+03_r8,     0.2964024E+03_r8,     0.2956741E+03_r8,     0.2949438E+03_r8,&  
         0.2942113E+03_r8,     0.2934765E+03_r8,     0.2927394E+03_r8,     0.2919998E+03_r8,     0.2912575E+03_r8,&  
         0.2905125E+03_r8,     0.2897643E+03_r8,     0.2890137E+03_r8,     0.2882596E+03_r8,     0.2875020E+03_r8,&  
         0.2867410E+03_r8,     0.2859762E+03_r8,     0.2852076E+03_r8,     0.2844347E+03_r8,     0.2836580E+03_r8,&  
         0.2828768E+03_r8,     0.2820909E+03_r8,     0.2813000E+03_r8,     0.2805042E+03_r8,     0.2797032E+03_r8,&  
         0.2788968E+03_r8,     0.2780845E+03_r8,     0.2772661E+03_r8,     0.2764417E+03_r8,     0.2756104E+03_r8,&  
         0.2747729E+03_r8,     0.2739280E+03_r8,     0.2730758E+03_r8,     0.2722160E+03_r8,     0.2713480E+03_r8,&  
         0.2704721E+03_r8,     0.2695874E+03_r8,     0.2686936E+03_r8,     0.2677909E+03_r8,     0.2668781E+03_r8,&  
         0.2659556E+03_r8,     0.2650227E+03_r8,     0.2640789E+03_r8,     0.2631243E+03_r8,     0.2621583E+03_r8,&  
         0.2611804E+03_r8,     0.2601904E+03_r8,     0.2591878E+03_r8,     0.2581725E+03_r8,     0.2571440E+03_r8,&  
         0.2561021E+03_r8,     0.2550465E+03_r8,     0.2539770E+03_r8,     0.2528933E+03_r8,     0.2517953E+03_r8,&  
         0.2506828E+03_r8,     0.2495559E+03_r8,     0.2484144E+03_r8,     0.2472584E+03_r8,     0.2460881E+03_r8,&  
         0.2449037E+03_r8,     0.2437052E+03_r8,     0.2424933E+03_r8,     0.2412681E+03_r8,     0.2400303E+03_r8,&  
         0.2387803E+03_r8,     0.2375189E+03_r8,     0.2362466E+03_r8,     0.2349644E+03_r8,     0.2336730E+03_r8,&  
         0.2323732E+03_r8,     0.2310659E+03_r8,     0.2297521E+03_r8,     0.2284326E+03_r8,     0.2271085E+03_r8,&  
         0.2257806E+03_r8,     0.2244499E+03_r8,     0.2231172E+03_r8,     0.2217833E+03_r8,     0.2204492E+03_r8,&  
         0.2191156E+03_r8,     0.2177831E+03_r8,     0.2164526E+03_r8,     0.2151245E+03_r8,     0.2137997E+03_r8,&  
         0.2124784E+03_r8,     0.2111613E+03_r8,     0.2098488E+03_r8,     0.2085412E+03_r8,     0.2072390E+03_r8,&  
         0.2059425E+03_r8,     0.2046518E+03_r8,     0.2033673E+03_r8,     0.2020893E+03_r8,     0.2008178E+03_r8,&  
         0.1995529E+03_r8,     0.1982950E+03_r8,     0.1970440E+03_r8,     0.1958001E+03_r8,     0.1945633E+03_r8,&  
         0.1933337E+03_r8,     0.1921113E+03_r8,     0.1908963E+03_r8,     0.1896884E+03_r8,     0.1884879E+03_r8,&  
         0.1872947E+03_r8,     0.1861088E+03_r8,     0.1849302E+03_r8,     0.1837589E+03_r8,     0.1825949E+03_r8,&  
         0.1814381E+03_r8,     0.1802885E+03_r8,     0.1791461E+03_r8,     0.1780109E+03_r8,     0.1768828E+03_r8,&  
         0.1757618E+03_r8,     0.1746479E+03_r8,     0.1735410E+03_r8,     0.1724411E+03_r8,     0.1713481E+03_r8,&  
         0.1702621E+03_r8,     0.1691828E+03_r8,     0.1681104E+03_r8,     0.1670449E+03_r8,     0.1659860E+03_r8,&  
         0.1649338E+03_r8,     0.1638884E+03_r8,     0.1628495E+03_r8,     0.1618172E+03_r8,     0.1607915E+03_r8,&  
         0.1597722E+03_r8,     0.1587595E+03_r8,     0.1577531E+03_r8,     0.1567531E+03_r8,     0.1557595E+03_r8,&  
         0.1547721E+03_r8,     0.1537910E+03_r8,     0.1528161E+03_r8,     0.1518474E+03_r8,     0.1508849E+03_r8/)  

    psaditmk(1:150,   79)= (/ &
         0.3055036E+03_r8,     0.3047968E+03_r8,     0.3040886E+03_r8,     0.3033801E+03_r8,     0.3026703E+03_r8,&  
         0.3019595E+03_r8,     0.3012477E+03_r8,     0.3005348E+03_r8,     0.2998205E+03_r8,     0.2991047E+03_r8,&  
         0.2983877E+03_r8,     0.2976689E+03_r8,     0.2969490E+03_r8,     0.2962270E+03_r8,     0.2955028E+03_r8,&  
         0.2947769E+03_r8,     0.2940491E+03_r8,     0.2933188E+03_r8,     0.2925864E+03_r8,     0.2918516E+03_r8,&  
         0.2911139E+03_r8,     0.2903739E+03_r8,     0.2896310E+03_r8,     0.2888849E+03_r8,     0.2881359E+03_r8,&  
         0.2873835E+03_r8,     0.2866277E+03_r8,     0.2858683E+03_r8,     0.2851053E+03_r8,     0.2843383E+03_r8,&  
         0.2835670E+03_r8,     0.2827915E+03_r8,     0.2820115E+03_r8,     0.2812268E+03_r8,     0.2804370E+03_r8,&  
         0.2796423E+03_r8,     0.2788422E+03_r8,     0.2780365E+03_r8,     0.2772248E+03_r8,     0.2764073E+03_r8,&  
         0.2755829E+03_r8,     0.2747525E+03_r8,     0.2739149E+03_r8,     0.2730701E+03_r8,     0.2722179E+03_r8,&  
         0.2713579E+03_r8,     0.2704897E+03_r8,     0.2696132E+03_r8,     0.2687277E+03_r8,     0.2678334E+03_r8,&  
         0.2669294E+03_r8,     0.2660157E+03_r8,     0.2650919E+03_r8,     0.2641576E+03_r8,     0.2632124E+03_r8,&  
         0.2622559E+03_r8,     0.2612878E+03_r8,     0.2603077E+03_r8,     0.2593154E+03_r8,     0.2583104E+03_r8,&  
         0.2572924E+03_r8,     0.2562611E+03_r8,     0.2552162E+03_r8,     0.2541575E+03_r8,     0.2530847E+03_r8,&  
         0.2519975E+03_r8,     0.2508960E+03_r8,     0.2497800E+03_r8,     0.2486494E+03_r8,     0.2475043E+03_r8,&  
         0.2463447E+03_r8,     0.2451707E+03_r8,     0.2439825E+03_r8,     0.2427807E+03_r8,     0.2415652E+03_r8,&  
         0.2403369E+03_r8,     0.2390959E+03_r8,     0.2378431E+03_r8,     0.2365791E+03_r8,     0.2353045E+03_r8,&  
         0.2340203E+03_r8,     0.2327272E+03_r8,     0.2314261E+03_r8,     0.2301179E+03_r8,     0.2288035E+03_r8,&  
         0.2274839E+03_r8,     0.2261600E+03_r8,     0.2248327E+03_r8,     0.2235029E+03_r8,     0.2221716E+03_r8,&  
         0.2208395E+03_r8,     0.2195075E+03_r8,     0.2181762E+03_r8,     0.2168464E+03_r8,     0.2155189E+03_r8,&  
         0.2141941E+03_r8,     0.2128728E+03_r8,     0.2115553E+03_r8,     0.2102421E+03_r8,     0.2089337E+03_r8,&  
         0.2076304E+03_r8,     0.2063326E+03_r8,     0.2050406E+03_r8,     0.2037547E+03_r8,     0.2024750E+03_r8,&  
         0.2012017E+03_r8,     0.1999352E+03_r8,     0.1986753E+03_r8,     0.1974224E+03_r8,     0.1961765E+03_r8,&  
         0.1949377E+03_r8,     0.1937060E+03_r8,     0.1924815E+03_r8,     0.1912643E+03_r8,     0.1900544E+03_r8,&  
         0.1888517E+03_r8,     0.1876563E+03_r8,     0.1864682E+03_r8,     0.1852875E+03_r8,     0.1841140E+03_r8,&  
         0.1829477E+03_r8,     0.1817888E+03_r8,     0.1806370E+03_r8,     0.1794925E+03_r8,     0.1783551E+03_r8,&  
         0.1772249E+03_r8,     0.1761017E+03_r8,     0.1749857E+03_r8,     0.1738766E+03_r8,     0.1727746E+03_r8,&  
         0.1716795E+03_r8,     0.1705914E+03_r8,     0.1695101E+03_r8,     0.1684356E+03_r8,     0.1673680E+03_r8,&  
         0.1663071E+03_r8,     0.1652529E+03_r8,     0.1642054E+03_r8,     0.1631645E+03_r8,     0.1621302E+03_r8,&  
         0.1611025E+03_r8,     0.1600813E+03_r8,     0.1590666E+03_r8,     0.1580583E+03_r8,     0.1570563E+03_r8,&  
         0.1560608E+03_r8,     0.1550715E+03_r8,     0.1540885E+03_r8,     0.1531118E+03_r8,     0.1521412E+03_r8/)  

    psaditmk(1:150,   80)= (/ &
         0.3059766E+03_r8,     0.3052743E+03_r8,     0.3045706E+03_r8,     0.3038666E+03_r8,     0.3031618E+03_r8,&  
         0.3024558E+03_r8,     0.3017490E+03_r8,     0.3010413E+03_r8,     0.3003321E+03_r8,     0.2996221E+03_r8,&  
         0.2989103E+03_r8,     0.2981974E+03_r8,     0.2974829E+03_r8,     0.2967669E+03_r8,     0.2960490E+03_r8,&  
         0.2953295E+03_r8,     0.2946081E+03_r8,     0.2938846E+03_r8,     0.2931592E+03_r8,     0.2924310E+03_r8,&  
         0.2917009E+03_r8,     0.2909683E+03_r8,     0.2902325E+03_r8,     0.2894947E+03_r8,     0.2887538E+03_r8,&  
         0.2880096E+03_r8,     0.2872625E+03_r8,     0.2865118E+03_r8,     0.2857578E+03_r8,     0.2850000E+03_r8,&  
         0.2842384E+03_r8,     0.2834728E+03_r8,     0.2827031E+03_r8,     0.2819288E+03_r8,     0.2811501E+03_r8,&  
         0.2803664E+03_r8,     0.2795779E+03_r8,     0.2787841E+03_r8,     0.2779845E+03_r8,     0.2771796E+03_r8,&  
         0.2763688E+03_r8,     0.2755513E+03_r8,     0.2747278E+03_r8,     0.2738975E+03_r8,     0.2730600E+03_r8,&  
         0.2722154E+03_r8,     0.2713631E+03_r8,     0.2705026E+03_r8,     0.2696341E+03_r8,     0.2687570E+03_r8,&  
         0.2678708E+03_r8,     0.2669757E+03_r8,     0.2660705E+03_r8,     0.2651558E+03_r8,     0.2642304E+03_r8,&  
         0.2632945E+03_r8,     0.2623476E+03_r8,     0.2613892E+03_r8,     0.2604190E+03_r8,     0.2594367E+03_r8,&  
         0.2584419E+03_r8,     0.2574343E+03_r8,     0.2564135E+03_r8,     0.2553792E+03_r8,     0.2543313E+03_r8,&  
         0.2532693E+03_r8,     0.2521932E+03_r8,     0.2511026E+03_r8,     0.2499975E+03_r8,     0.2488779E+03_r8,&  
         0.2477436E+03_r8,     0.2465947E+03_r8,     0.2454314E+03_r8,     0.2442537E+03_r8,     0.2430619E+03_r8,&  
         0.2418565E+03_r8,     0.2406377E+03_r8,     0.2394060E+03_r8,     0.2381621E+03_r8,     0.2369064E+03_r8,&  
         0.2356400E+03_r8,     0.2343633E+03_r8,     0.2330771E+03_r8,     0.2317825E+03_r8,     0.2304801E+03_r8,&  
         0.2291712E+03_r8,     0.2278563E+03_r8,     0.2265367E+03_r8,     0.2252132E+03_r8,     0.2238866E+03_r8,&  
         0.2225581E+03_r8,     0.2212282E+03_r8,     0.2198980E+03_r8,     0.2185681E+03_r8,     0.2172395E+03_r8,&  
         0.2159125E+03_r8,     0.2145881E+03_r8,     0.2132667E+03_r8,     0.2119489E+03_r8,     0.2106352E+03_r8,&  
         0.2093261E+03_r8,     0.2080218E+03_r8,     0.2067229E+03_r8,     0.2054297E+03_r8,     0.2041423E+03_r8,&  
         0.2028610E+03_r8,     0.2015862E+03_r8,     0.2003178E+03_r8,     0.1990561E+03_r8,     0.1978013E+03_r8,&  
         0.1965534E+03_r8,     0.1953126E+03_r8,     0.1940789E+03_r8,     0.1928523E+03_r8,     0.1916329E+03_r8,&  
         0.1904209E+03_r8,     0.1892161E+03_r8,     0.1880185E+03_r8,     0.1868283E+03_r8,     0.1856453E+03_r8,&  
         0.1844697E+03_r8,     0.1833013E+03_r8,     0.1821401E+03_r8,     0.1809862E+03_r8,     0.1798395E+03_r8,&  
         0.1786999E+03_r8,     0.1775676E+03_r8,     0.1764423E+03_r8,     0.1753241E+03_r8,     0.1742129E+03_r8,&  
         0.1731088E+03_r8,     0.1720116E+03_r8,     0.1709213E+03_r8,     0.1698380E+03_r8,     0.1687614E+03_r8,&  
         0.1676917E+03_r8,     0.1666288E+03_r8,     0.1655726E+03_r8,     0.1645230E+03_r8,     0.1634801E+03_r8,&  
         0.1624439E+03_r8,     0.1614142E+03_r8,     0.1603910E+03_r8,     0.1593743E+03_r8,     0.1583640E+03_r8,&  
         0.1573602E+03_r8,     0.1563627E+03_r8,     0.1553715E+03_r8,     0.1543866E+03_r8,     0.1534079E+03_r8/)  

    psaditmk(1:150,   81)= (/ &
         0.3064403E+03_r8,     0.3057421E+03_r8,     0.3050430E+03_r8,     0.3043433E+03_r8,     0.3036431E+03_r8,&  
         0.3029418E+03_r8,     0.3022397E+03_r8,     0.3015370E+03_r8,     0.3008330E+03_r8,     0.3001282E+03_r8,&  
         0.2994218E+03_r8,     0.2987141E+03_r8,     0.2980052E+03_r8,     0.2972948E+03_r8,     0.2965831E+03_r8,&  
         0.2958697E+03_r8,     0.2951542E+03_r8,     0.2944371E+03_r8,     0.2937180E+03_r8,     0.2929968E+03_r8,&  
         0.2922735E+03_r8,     0.2915477E+03_r8,     0.2908196E+03_r8,     0.2900892E+03_r8,     0.2893559E+03_r8,&  
         0.2886198E+03_r8,     0.2878807E+03_r8,     0.2871384E+03_r8,     0.2863931E+03_r8,     0.2856443E+03_r8,&  
         0.2848919E+03_r8,     0.2841356E+03_r8,     0.2833755E+03_r8,     0.2826114E+03_r8,     0.2818430E+03_r8,&  
         0.2810699E+03_r8,     0.2802924E+03_r8,     0.2795099E+03_r8,     0.2787223E+03_r8,     0.2779291E+03_r8,&  
         0.2771306E+03_r8,     0.2763264E+03_r8,     0.2755157E+03_r8,     0.2746991E+03_r8,     0.2738758E+03_r8,&  
         0.2730454E+03_r8,     0.2722082E+03_r8,     0.2713635E+03_r8,     0.2705107E+03_r8,     0.2696500E+03_r8,&  
         0.2687810E+03_r8,     0.2679030E+03_r8,     0.2670163E+03_r8,     0.2661200E+03_r8,     0.2652137E+03_r8,&  
         0.2642977E+03_r8,     0.2633710E+03_r8,     0.2624333E+03_r8,     0.2614846E+03_r8,     0.2605242E+03_r8,&  
         0.2595519E+03_r8,     0.2585671E+03_r8,     0.2575698E+03_r8,     0.2565595E+03_r8,     0.2555358E+03_r8,&  
         0.2544986E+03_r8,     0.2534473E+03_r8,     0.2523821E+03_r8,     0.2513025E+03_r8,     0.2502084E+03_r8,&  
         0.2490997E+03_r8,     0.2479763E+03_r8,     0.2468382E+03_r8,     0.2456855E+03_r8,     0.2445184E+03_r8,&  
         0.2433371E+03_r8,     0.2421417E+03_r8,     0.2409328E+03_r8,     0.2397105E+03_r8,     0.2384757E+03_r8,&  
         0.2372288E+03_r8,     0.2359705E+03_r8,     0.2347015E+03_r8,     0.2334226E+03_r8,     0.2321348E+03_r8,&  
         0.2308387E+03_r8,     0.2295353E+03_r8,     0.2282257E+03_r8,     0.2269106E+03_r8,     0.2255912E+03_r8,&  
         0.2242681E+03_r8,     0.2229426E+03_r8,     0.2216152E+03_r8,     0.2202871E+03_r8,     0.2189589E+03_r8,&  
         0.2176314E+03_r8,     0.2163053E+03_r8,     0.2149814E+03_r8,     0.2136601E+03_r8,     0.2123422E+03_r8,&  
         0.2110282E+03_r8,     0.2097184E+03_r8,     0.2084133E+03_r8,     0.2071134E+03_r8,     0.2058189E+03_r8,&  
         0.2045301E+03_r8,     0.2032474E+03_r8,     0.2019709E+03_r8,     0.2007008E+03_r8,     0.1994374E+03_r8,&  
         0.1981807E+03_r8,     0.1969309E+03_r8,     0.1956881E+03_r8,     0.1944523E+03_r8,     0.1932237E+03_r8,&  
         0.1920023E+03_r8,     0.1907880E+03_r8,     0.1895811E+03_r8,     0.1883814E+03_r8,     0.1871890E+03_r8,&  
         0.1860039E+03_r8,     0.1848260E+03_r8,     0.1836554E+03_r8,     0.1824921E+03_r8,     0.1813360E+03_r8,&  
         0.1801871E+03_r8,     0.1790454E+03_r8,     0.1779109E+03_r8,     0.1767834E+03_r8,     0.1756631E+03_r8,&  
         0.1745498E+03_r8,     0.1734436E+03_r8,     0.1723442E+03_r8,     0.1712519E+03_r8,     0.1701665E+03_r8,&  
         0.1690879E+03_r8,     0.1680161E+03_r8,     0.1669511E+03_r8,     0.1658928E+03_r8,     0.1648413E+03_r8,&  
         0.1637964E+03_r8,     0.1627581E+03_r8,     0.1617264E+03_r8,     0.1607012E+03_r8,     0.1596826E+03_r8,&  
         0.1586704E+03_r8,     0.1576646E+03_r8,     0.1566651E+03_r8,     0.1556721E+03_r8,     0.1546852E+03_r8/)  

    psaditmk(1:150,   82)= (/ &
         0.3068949E+03_r8,     0.3062008E+03_r8,     0.3055060E+03_r8,     0.3048105E+03_r8,     0.3041147E+03_r8,&  
         0.3034179E+03_r8,     0.3027204E+03_r8,     0.3020226E+03_r8,     0.3013232E+03_r8,     0.3006235E+03_r8,&  
         0.2999220E+03_r8,     0.2992198E+03_r8,     0.2985163E+03_r8,     0.2978111E+03_r8,     0.2971050E+03_r8,&  
         0.2963975E+03_r8,     0.2956879E+03_r8,     0.2949770E+03_r8,     0.2942640E+03_r8,     0.2935495E+03_r8,&  
         0.2928329E+03_r8,     0.2921136E+03_r8,     0.2913926E+03_r8,     0.2906690E+03_r8,     0.2899431E+03_r8,&  
         0.2892147E+03_r8,     0.2884833E+03_r8,     0.2877492E+03_r8,     0.2870120E+03_r8,     0.2862716E+03_r8,&  
         0.2855278E+03_r8,     0.2847807E+03_r8,     0.2840297E+03_r8,     0.2832754E+03_r8,     0.2825168E+03_r8,&  
         0.2817538E+03_r8,     0.2809866E+03_r8,     0.2802149E+03_r8,     0.2794383E+03_r8,     0.2786566E+03_r8,&  
         0.2778699E+03_r8,     0.2770778E+03_r8,     0.2762801E+03_r8,     0.2754761E+03_r8,     0.2746662E+03_r8,&  
         0.2738499E+03_r8,     0.2730266E+03_r8,     0.2721966E+03_r8,     0.2713593E+03_r8,     0.2705142E+03_r8,&  
         0.2696613E+03_r8,     0.2688002E+03_r8,     0.2679305E+03_r8,     0.2670518E+03_r8,     0.2661643E+03_r8,&  
         0.2652667E+03_r8,     0.2643594E+03_r8,     0.2634418E+03_r8,     0.2625135E+03_r8,     0.2615742E+03_r8,&  
         0.2606234E+03_r8,     0.2596609E+03_r8,     0.2586862E+03_r8,     0.2576991E+03_r8,     0.2566991E+03_r8,&  
         0.2556859E+03_r8,     0.2546593E+03_r8,     0.2536188E+03_r8,     0.2525644E+03_r8,     0.2514957E+03_r8,&  
         0.2504125E+03_r8,     0.2493148E+03_r8,     0.2482023E+03_r8,     0.2470751E+03_r8,     0.2459332E+03_r8,&  
         0.2447768E+03_r8,     0.2436059E+03_r8,     0.2424207E+03_r8,     0.2412217E+03_r8,     0.2400093E+03_r8,&  
         0.2387837E+03_r8,     0.2375458E+03_r8,     0.2362960E+03_r8,     0.2350351E+03_r8,     0.2337638E+03_r8,&  
         0.2324830E+03_r8,     0.2311933E+03_r8,     0.2298960E+03_r8,     0.2285918E+03_r8,     0.2272816E+03_r8,&  
         0.2259664E+03_r8,     0.2246473E+03_r8,     0.2233249E+03_r8,     0.2220005E+03_r8,     0.2206746E+03_r8,&  
         0.2193482E+03_r8,     0.2180222E+03_r8,     0.2166972E+03_r8,     0.2153739E+03_r8,     0.2140530E+03_r8,&  
         0.2127352E+03_r8,     0.2114208E+03_r8,     0.2101105E+03_r8,     0.2088047E+03_r8,     0.2075038E+03_r8,&  
         0.2062082E+03_r8,     0.2049182E+03_r8,     0.2036340E+03_r8,     0.2023560E+03_r8,     0.2010843E+03_r8,&  
         0.1998191E+03_r8,     0.1985605E+03_r8,     0.1973088E+03_r8,     0.1960641E+03_r8,     0.1948263E+03_r8,&  
         0.1935956E+03_r8,     0.1923721E+03_r8,     0.1911558E+03_r8,     0.1899467E+03_r8,     0.1887449E+03_r8,&  
         0.1875503E+03_r8,     0.1863630E+03_r8,     0.1851830E+03_r8,     0.1840102E+03_r8,     0.1828448E+03_r8,&  
         0.1816865E+03_r8,     0.1805354E+03_r8,     0.1793916E+03_r8,     0.1782549E+03_r8,     0.1771253E+03_r8,&  
         0.1760028E+03_r8,     0.1748874E+03_r8,     0.1737790E+03_r8,     0.1726776E+03_r8,     0.1715831E+03_r8,&  
         0.1704956E+03_r8,     0.1694149E+03_r8,     0.1683411E+03_r8,     0.1672740E+03_r8,     0.1662137E+03_r8,&  
         0.1651601E+03_r8,     0.1641132E+03_r8,     0.1630730E+03_r8,     0.1620392E+03_r8,     0.1610121E+03_r8,&  
         0.1599915E+03_r8,     0.1589773E+03_r8,     0.1579696E+03_r8,     0.1569682E+03_r8,     0.1559732E+03_r8/)  

    psaditmk(1:150,   83)= (/ &
         0.3073409E+03_r8,     0.3066503E+03_r8,     0.3059601E+03_r8,     0.3052687E+03_r8,     0.3045771E+03_r8,&  
         0.3038846E+03_r8,     0.3031917E+03_r8,     0.3024984E+03_r8,     0.3018034E+03_r8,     0.3011083E+03_r8,&  
         0.3004117E+03_r8,     0.2997146E+03_r8,     0.2990161E+03_r8,     0.2983166E+03_r8,     0.2976158E+03_r8,&  
         0.2969134E+03_r8,     0.2962097E+03_r8,     0.2955046E+03_r8,     0.2947976E+03_r8,     0.2940891E+03_r8,&  
         0.2933786E+03_r8,     0.2926663E+03_r8,     0.2919518E+03_r8,     0.2912349E+03_r8,     0.2905162E+03_r8,&  
         0.2897948E+03_r8,     0.2890710E+03_r8,     0.2883445E+03_r8,     0.2876151E+03_r8,     0.2868829E+03_r8,&  
         0.2861474E+03_r8,     0.2854089E+03_r8,     0.2846669E+03_r8,     0.2839212E+03_r8,     0.2831719E+03_r8,&  
         0.2824188E+03_r8,     0.2816615E+03_r8,     0.2809000E+03_r8,     0.2801340E+03_r8,     0.2793632E+03_r8,&  
         0.2785878E+03_r8,     0.2778074E+03_r8,     0.2770215E+03_r8,     0.2762300E+03_r8,     0.2754325E+03_r8,&  
         0.2746294E+03_r8,     0.2738199E+03_r8,     0.2730036E+03_r8,     0.2721807E+03_r8,     0.2713508E+03_r8,&  
         0.2705131E+03_r8,     0.2696678E+03_r8,     0.2688147E+03_r8,     0.2679531E+03_r8,     0.2670822E+03_r8,&  
         0.2662030E+03_r8,     0.2653142E+03_r8,     0.2644157E+03_r8,     0.2635070E+03_r8,     0.2625878E+03_r8,&  
         0.2616579E+03_r8,     0.2607167E+03_r8,     0.2597639E+03_r8,     0.2587991E+03_r8,     0.2578220E+03_r8,&  
         0.2568322E+03_r8,     0.2558294E+03_r8,     0.2548133E+03_r8,     0.2537836E+03_r8,     0.2527399E+03_r8,&  
         0.2516821E+03_r8,     0.2506098E+03_r8,     0.2495231E+03_r8,     0.2484216E+03_r8,     0.2473053E+03_r8,&  
         0.2461744E+03_r8,     0.2450286E+03_r8,     0.2438683E+03_r8,     0.2426936E+03_r8,     0.2415047E+03_r8,&  
         0.2403021E+03_r8,     0.2390862E+03_r8,     0.2378574E+03_r8,     0.2366163E+03_r8,     0.2353637E+03_r8,&  
         0.2341003E+03_r8,     0.2328268E+03_r8,     0.2315441E+03_r8,     0.2302530E+03_r8,     0.2289544E+03_r8,&  
         0.2276495E+03_r8,     0.2263391E+03_r8,     0.2250238E+03_r8,     0.2237051E+03_r8,     0.2223837E+03_r8,&  
         0.2210604E+03_r8,     0.2197361E+03_r8,     0.2184117E+03_r8,     0.2170880E+03_r8,     0.2157656E+03_r8,&  
         0.2144452E+03_r8,     0.2131275E+03_r8,     0.2118131E+03_r8,     0.2105025E+03_r8,     0.2091960E+03_r8,&  
         0.2078943E+03_r8,     0.2065977E+03_r8,     0.2053064E+03_r8,     0.2040209E+03_r8,     0.2027414E+03_r8,&  
         0.2014680E+03_r8,     0.2002012E+03_r8,     0.1989409E+03_r8,     0.1976873E+03_r8,     0.1964406E+03_r8,&  
         0.1952008E+03_r8,     0.1939681E+03_r8,     0.1927426E+03_r8,     0.1915241E+03_r8,     0.1903129E+03_r8,&  
         0.1891090E+03_r8,     0.1879123E+03_r8,     0.1867228E+03_r8,     0.1855406E+03_r8,     0.1843657E+03_r8,&  
         0.1831980E+03_r8,     0.1820376E+03_r8,     0.1808844E+03_r8,     0.1797384E+03_r8,     0.1785995E+03_r8,&  
         0.1774678E+03_r8,     0.1763431E+03_r8,     0.1752256E+03_r8,     0.1741151E+03_r8,     0.1730115E+03_r8,&  
         0.1719150E+03_r8,     0.1708254E+03_r8,     0.1697426E+03_r8,     0.1686667E+03_r8,     0.1675976E+03_r8,&  
         0.1665352E+03_r8,     0.1654796E+03_r8,     0.1644307E+03_r8,     0.1633884E+03_r8,     0.1623527E+03_r8,&  
         0.1613236E+03_r8,     0.1603009E+03_r8,     0.1592848E+03_r8,     0.1582751E+03_r8,     0.1572718E+03_r8/)  

    psaditmk(1:150,   84)= (/ &
         0.3077788E+03_r8,     0.3070919E+03_r8,     0.3064050E+03_r8,     0.3057179E+03_r8,     0.3050303E+03_r8,&  
         0.3043423E+03_r8,     0.3036534E+03_r8,     0.3029641E+03_r8,     0.3022738E+03_r8,     0.3015832E+03_r8,&  
         0.3008919E+03_r8,     0.3001991E+03_r8,     0.2995057E+03_r8,     0.2988110E+03_r8,     0.2981155E+03_r8,&  
         0.2974183E+03_r8,     0.2967200E+03_r8,     0.2960205E+03_r8,     0.2953192E+03_r8,     0.2946166E+03_r8,&  
         0.2939121E+03_r8,     0.2932062E+03_r8,     0.2924980E+03_r8,     0.2917876E+03_r8,     0.2910756E+03_r8,&  
         0.2903608E+03_r8,     0.2896443E+03_r8,     0.2889249E+03_r8,     0.2882032E+03_r8,     0.2874787E+03_r8,&  
         0.2867512E+03_r8,     0.2860208E+03_r8,     0.2852872E+03_r8,     0.2845502E+03_r8,     0.2838098E+03_r8,&  
         0.2830657E+03_r8,     0.2823180E+03_r8,     0.2815663E+03_r8,     0.2808103E+03_r8,     0.2800498E+03_r8,&  
         0.2792851E+03_r8,     0.2785156E+03_r8,     0.2777411E+03_r8,     0.2769613E+03_r8,     0.2761759E+03_r8,&  
         0.2753854E+03_r8,     0.2745887E+03_r8,     0.2737859E+03_r8,     0.2729766E+03_r8,     0.2721607E+03_r8,&  
         0.2713380E+03_r8,     0.2705076E+03_r8,     0.2696698E+03_r8,     0.2688242E+03_r8,     0.2679706E+03_r8,&  
         0.2671082E+03_r8,     0.2662370E+03_r8,     0.2653568E+03_r8,     0.2644666E+03_r8,     0.2635669E+03_r8,&  
         0.2626570E+03_r8,     0.2617360E+03_r8,     0.2608042E+03_r8,     0.2598610E+03_r8,     0.2589060E+03_r8,&  
         0.2579389E+03_r8,     0.2569593E+03_r8,     0.2559668E+03_r8,     0.2549612E+03_r8,     0.2539420E+03_r8,&  
         0.2529091E+03_r8,     0.2518621E+03_r8,     0.2508008E+03_r8,     0.2497249E+03_r8,     0.2486343E+03_r8,&  
         0.2475291E+03_r8,     0.2464090E+03_r8,     0.2452740E+03_r8,     0.2441244E+03_r8,     0.2429602E+03_r8,&  
         0.2417817E+03_r8,     0.2405891E+03_r8,     0.2393829E+03_r8,     0.2381635E+03_r8,     0.2369315E+03_r8,&  
         0.2356875E+03_r8,     0.2344321E+03_r8,     0.2331663E+03_r8,     0.2318906E+03_r8,     0.2306061E+03_r8,&  
         0.2293137E+03_r8,     0.2280142E+03_r8,     0.2267086E+03_r8,     0.2253979E+03_r8,     0.2240830E+03_r8,&  
         0.2227649E+03_r8,     0.2214445E+03_r8,     0.2201226E+03_r8,     0.2188000E+03_r8,     0.2174778E+03_r8,&  
         0.2161564E+03_r8,     0.2148368E+03_r8,     0.2135194E+03_r8,     0.2122050E+03_r8,     0.2108942E+03_r8,&  
         0.2095873E+03_r8,     0.2082848E+03_r8,     0.2069872E+03_r8,     0.2056949E+03_r8,     0.2044081E+03_r8,&  
         0.2031272E+03_r8,     0.2018523E+03_r8,     0.2005838E+03_r8,     0.1993217E+03_r8,     0.1980663E+03_r8,&  
         0.1968177E+03_r8,     0.1955760E+03_r8,     0.1943413E+03_r8,     0.1931137E+03_r8,     0.1918932E+03_r8,&  
         0.1906799E+03_r8,     0.1894738E+03_r8,     0.1882749E+03_r8,     0.1870833E+03_r8,     0.1858990E+03_r8,&  
         0.1847219E+03_r8,     0.1835520E+03_r8,     0.1823895E+03_r8,     0.1812341E+03_r8,     0.1800859E+03_r8,&  
         0.1789449E+03_r8,     0.1778110E+03_r8,     0.1766842E+03_r8,     0.1755645E+03_r8,     0.1744519E+03_r8,&  
         0.1733462E+03_r8,     0.1722476E+03_r8,     0.1711558E+03_r8,     0.1700710E+03_r8,     0.1689930E+03_r8,&  
         0.1679218E+03_r8,     0.1668574E+03_r8,     0.1657998E+03_r8,     0.1647488E+03_r8,     0.1637045E+03_r8,&  
         0.1626668E+03_r8,     0.1616357E+03_r8,     0.1606111E+03_r8,     0.1595930E+03_r8,     0.1585814E+03_r8/)  

    psaditmk(1:150,   85)= (/ &
         0.3082078E+03_r8,     0.3075247E+03_r8,     0.3068419E+03_r8,     0.3061585E+03_r8,     0.3054746E+03_r8,&  
         0.3047906E+03_r8,     0.3041060E+03_r8,     0.3034208E+03_r8,     0.3027349E+03_r8,     0.3020488E+03_r8,&  
         0.3013615E+03_r8,     0.3006735E+03_r8,     0.2999850E+03_r8,     0.2992951E+03_r8,     0.2986044E+03_r8,&  
         0.2979124E+03_r8,     0.2972192E+03_r8,     0.2965251E+03_r8,     0.2958292E+03_r8,     0.2951323E+03_r8,&  
         0.2944336E+03_r8,     0.2937333E+03_r8,     0.2930312E+03_r8,     0.2923275E+03_r8,     0.2916216E+03_r8,&  
         0.2909135E+03_r8,     0.2902039E+03_r8,     0.2894915E+03_r8,     0.2887769E+03_r8,     0.2880595E+03_r8,&  
         0.2873399E+03_r8,     0.2866169E+03_r8,     0.2858913E+03_r8,     0.2851627E+03_r8,     0.2844306E+03_r8,&  
         0.2836957E+03_r8,     0.2829569E+03_r8,     0.2822144E+03_r8,     0.2814679E+03_r8,     0.2807174E+03_r8,&  
         0.2799627E+03_r8,     0.2792036E+03_r8,     0.2784400E+03_r8,     0.2776714E+03_r8,     0.2768977E+03_r8,&  
         0.2761186E+03_r8,     0.2753342E+03_r8,     0.2745441E+03_r8,     0.2737479E+03_r8,     0.2729453E+03_r8,&  
         0.2721364E+03_r8,     0.2713207E+03_r8,     0.2704977E+03_r8,     0.2696673E+03_r8,     0.2688293E+03_r8,&  
         0.2679833E+03_r8,     0.2671289E+03_r8,     0.2662659E+03_r8,     0.2653938E+03_r8,     0.2645125E+03_r8,&  
         0.2636215E+03_r8,     0.2627200E+03_r8,     0.2618083E+03_r8,     0.2608859E+03_r8,     0.2599521E+03_r8,&  
         0.2590068E+03_r8,     0.2580496E+03_r8,     0.2570801E+03_r8,     0.2560978E+03_r8,     0.2551026E+03_r8,&  
         0.2540940E+03_r8,     0.2530717E+03_r8,     0.2520354E+03_r8,     0.2509849E+03_r8,     0.2499199E+03_r8,&  
         0.2488404E+03_r8,     0.2477461E+03_r8,     0.2466369E+03_r8,     0.2455129E+03_r8,     0.2443740E+03_r8,&  
         0.2432204E+03_r8,     0.2420523E+03_r8,     0.2408700E+03_r8,     0.2396735E+03_r8,     0.2384639E+03_r8,&  
         0.2372411E+03_r8,     0.2360058E+03_r8,     0.2347590E+03_r8,     0.2335010E+03_r8,     0.2322328E+03_r8,&  
         0.2309552E+03_r8,     0.2296691E+03_r8,     0.2283754E+03_r8,     0.2270750E+03_r8,     0.2257691E+03_r8,&  
         0.2244583E+03_r8,     0.2231438E+03_r8,     0.2218265E+03_r8,     0.2205071E+03_r8,     0.2191868E+03_r8,&  
         0.2178662E+03_r8,     0.2165461E+03_r8,     0.2152274E+03_r8,     0.2139106E+03_r8,     0.2125964E+03_r8,&  
         0.2112854E+03_r8,     0.2099781E+03_r8,     0.2086751E+03_r8,     0.2073767E+03_r8,     0.2060833E+03_r8,&  
         0.2047954E+03_r8,     0.2035131E+03_r8,     0.2022367E+03_r8,     0.2009666E+03_r8,     0.1997029E+03_r8,&  
         0.1984457E+03_r8,     0.1971952E+03_r8,     0.1959516E+03_r8,     0.1947149E+03_r8,     0.1934853E+03_r8,&  
         0.1922627E+03_r8,     0.1910473E+03_r8,     0.1898391E+03_r8,     0.1886381E+03_r8,     0.1874444E+03_r8,&  
         0.1862579E+03_r8,     0.1850786E+03_r8,     0.1839066E+03_r8,     0.1827419E+03_r8,     0.1815843E+03_r8,&  
         0.1804339E+03_r8,     0.1792908E+03_r8,     0.1781547E+03_r8,     0.1770258E+03_r8,     0.1759040E+03_r8,&  
         0.1747892E+03_r8,     0.1736815E+03_r8,     0.1725807E+03_r8,     0.1714868E+03_r8,     0.1703999E+03_r8,&  
         0.1693199E+03_r8,     0.1682466E+03_r8,     0.1671802E+03_r8,     0.1661205E+03_r8,     0.1650675E+03_r8,&  
         0.1640212E+03_r8,     0.1629815E+03_r8,     0.1619484E+03_r8,     0.1609218E+03_r8,     0.1599017E+03_r8/)  

    psaditmk(1:150,   86)= (/ &
         0.3086291E+03_r8,     0.3079499E+03_r8,     0.3072703E+03_r8,     0.3065905E+03_r8,     0.3059104E+03_r8,&  
         0.3052306E+03_r8,     0.3045495E+03_r8,     0.3038688E+03_r8,     0.3031870E+03_r8,     0.3025047E+03_r8,&  
         0.3018221E+03_r8,     0.3011383E+03_r8,     0.3004542E+03_r8,     0.2997690E+03_r8,     0.2990833E+03_r8,&  
         0.2983961E+03_r8,     0.2977081E+03_r8,     0.2970188E+03_r8,     0.2963284E+03_r8,     0.2956366E+03_r8,&  
         0.2949433E+03_r8,     0.2942489E+03_r8,     0.2935526E+03_r8,     0.2928550E+03_r8,     0.2921552E+03_r8,&  
         0.2914533E+03_r8,     0.2907500E+03_r8,     0.2900444E+03_r8,     0.2893365E+03_r8,     0.2886263E+03_r8,&  
         0.2879137E+03_r8,     0.2871983E+03_r8,     0.2864807E+03_r8,     0.2857596E+03_r8,     0.2850358E+03_r8,&  
         0.2843089E+03_r8,     0.2835786E+03_r8,     0.2828449E+03_r8,     0.2821077E+03_r8,     0.2813666E+03_r8,&  
         0.2806215E+03_r8,     0.2798724E+03_r8,     0.2791189E+03_r8,     0.2783609E+03_r8,     0.2775982E+03_r8,&  
         0.2768303E+03_r8,     0.2760577E+03_r8,     0.2752796E+03_r8,     0.2744958E+03_r8,     0.2737060E+03_r8,&  
         0.2729101E+03_r8,     0.2721080E+03_r8,     0.2712993E+03_r8,     0.2704834E+03_r8,     0.2696605E+03_r8,&  
         0.2688297E+03_r8,     0.2679914E+03_r8,     0.2671451E+03_r8,     0.2662899E+03_r8,     0.2654260E+03_r8,&  
         0.2645530E+03_r8,     0.2636702E+03_r8,     0.2627778E+03_r8,     0.2618753E+03_r8,     0.2609619E+03_r8,&  
         0.2600375E+03_r8,     0.2591018E+03_r8,     0.2581543E+03_r8,     0.2571946E+03_r8,     0.2562225E+03_r8,&  
         0.2552374E+03_r8,     0.2542392E+03_r8,     0.2532275E+03_r8,     0.2522019E+03_r8,     0.2511622E+03_r8,&  
         0.2501084E+03_r8,     0.2490398E+03_r8,     0.2479564E+03_r8,     0.2468582E+03_r8,     0.2457451E+03_r8,&  
         0.2446171E+03_r8,     0.2434742E+03_r8,     0.2423166E+03_r8,     0.2411446E+03_r8,     0.2399585E+03_r8,&  
         0.2387584E+03_r8,     0.2375451E+03_r8,     0.2363190E+03_r8,     0.2350807E+03_r8,     0.2338309E+03_r8,&  
         0.2325704E+03_r8,     0.2313000E+03_r8,     0.2300205E+03_r8,     0.2287330E+03_r8,     0.2274382E+03_r8,&  
         0.2261371E+03_r8,     0.2248309E+03_r8,     0.2235203E+03_r8,     0.2222063E+03_r8,     0.2208899E+03_r8,&  
         0.2195719E+03_r8,     0.2182532E+03_r8,     0.2169347E+03_r8,     0.2156170E+03_r8,     0.2143010E+03_r8,&  
         0.2129872E+03_r8,     0.2116762E+03_r8,     0.2103688E+03_r8,     0.2090652E+03_r8,     0.2077661E+03_r8,&  
         0.2064718E+03_r8,     0.2051828E+03_r8,     0.2038992E+03_r8,     0.2026215E+03_r8,     0.2013498E+03_r8,&  
         0.2000844E+03_r8,     0.1988255E+03_r8,     0.1975732E+03_r8,     0.1963277E+03_r8,     0.1950891E+03_r8,&  
         0.1938574E+03_r8,     0.1926328E+03_r8,     0.1914154E+03_r8,     0.1902050E+03_r8,     0.1890020E+03_r8,&  
         0.1878061E+03_r8,     0.1866174E+03_r8,     0.1854360E+03_r8,     0.1842618E+03_r8,     0.1830949E+03_r8,&  
         0.1819352E+03_r8,     0.1807827E+03_r8,     0.1796373E+03_r8,     0.1784991E+03_r8,     0.1773681E+03_r8,&  
         0.1762441E+03_r8,     0.1751272E+03_r8,     0.1740173E+03_r8,     0.1729144E+03_r8,     0.1718185E+03_r8,&  
         0.1707295E+03_r8,     0.1696474E+03_r8,     0.1685721E+03_r8,     0.1675036E+03_r8,     0.1664418E+03_r8,&  
         0.1653868E+03_r8,     0.1643384E+03_r8,     0.1632967E+03_r8,     0.1622616E+03_r8,     0.1612331E+03_r8/)  

    psaditmk(1:150,   87)= (/ &
         0.3090429E+03_r8,     0.3083672E+03_r8,     0.3076909E+03_r8,     0.3070149E+03_r8,     0.3063383E+03_r8,&  
         0.3056617E+03_r8,     0.3049852E+03_r8,     0.3043076E+03_r8,     0.3036301E+03_r8,     0.3029518E+03_r8,&  
         0.3022733E+03_r8,     0.3015938E+03_r8,     0.3009141E+03_r8,     0.3002333E+03_r8,     0.2995521E+03_r8,&  
         0.2988696E+03_r8,     0.2981863E+03_r8,     0.2975021E+03_r8,     0.2968166E+03_r8,     0.2961301E+03_r8,&  
         0.2954421E+03_r8,     0.2947531E+03_r8,     0.2940622E+03_r8,     0.2933703E+03_r8,     0.2926764E+03_r8,&  
         0.2919808E+03_r8,     0.2912838E+03_r8,     0.2905842E+03_r8,     0.2898829E+03_r8,     0.2891794E+03_r8,&  
         0.2884737E+03_r8,     0.2877654E+03_r8,     0.2870549E+03_r8,     0.2863414E+03_r8,     0.2856252E+03_r8,&  
         0.2849063E+03_r8,     0.2841842E+03_r8,     0.2834587E+03_r8,     0.2827303E+03_r8,     0.2819985E+03_r8,&  
         0.2812624E+03_r8,     0.2805227E+03_r8,     0.2797789E+03_r8,     0.2790310E+03_r8,     0.2782785E+03_r8,&  
         0.2775217E+03_r8,     0.2767600E+03_r8,     0.2759932E+03_r8,     0.2752213E+03_r8,     0.2744437E+03_r8,&  
         0.2736602E+03_r8,     0.2728711E+03_r8,     0.2720756E+03_r8,     0.2712738E+03_r8,     0.2704649E+03_r8,&  
         0.2696492E+03_r8,     0.2688257E+03_r8,     0.2679949E+03_r8,     0.2671562E+03_r8,     0.2663091E+03_r8,&  
         0.2654533E+03_r8,     0.2645884E+03_r8,     0.2637143E+03_r8,     0.2628306E+03_r8,     0.2619366E+03_r8,&  
         0.2610323E+03_r8,     0.2601172E+03_r8,     0.2591909E+03_r8,     0.2582531E+03_r8,     0.2573033E+03_r8,&  
         0.2563411E+03_r8,     0.2553663E+03_r8,     0.2543784E+03_r8,     0.2533771E+03_r8,     0.2523622E+03_r8,&  
         0.2513333E+03_r8,     0.2502901E+03_r8,     0.2492324E+03_r8,     0.2481600E+03_r8,     0.2470728E+03_r8,&  
         0.2459707E+03_r8,     0.2448536E+03_r8,     0.2437214E+03_r8,     0.2425746E+03_r8,     0.2414131E+03_r8,&  
         0.2402371E+03_r8,     0.2390471E+03_r8,     0.2378435E+03_r8,     0.2366266E+03_r8,     0.2353972E+03_r8,&  
         0.2341559E+03_r8,     0.2329034E+03_r8,     0.2316404E+03_r8,     0.2303679E+03_r8,     0.2290868E+03_r8,&  
         0.2277979E+03_r8,     0.2265021E+03_r8,     0.2252006E+03_r8,     0.2238942E+03_r8,     0.2225839E+03_r8,&  
         0.2212706E+03_r8,     0.2199553E+03_r8,     0.2186388E+03_r8,     0.2173220E+03_r8,     0.2160056E+03_r8,&  
         0.2146905E+03_r8,     0.2133773E+03_r8,     0.2120666E+03_r8,     0.2107590E+03_r8,     0.2094552E+03_r8,&  
         0.2081554E+03_r8,     0.2068603E+03_r8,     0.2055702E+03_r8,     0.2042855E+03_r8,     0.2030064E+03_r8,&  
         0.2017333E+03_r8,     0.2004663E+03_r8,     0.1992057E+03_r8,     0.1979516E+03_r8,     0.1967043E+03_r8,&  
         0.1954637E+03_r8,     0.1942301E+03_r8,     0.1930035E+03_r8,     0.1917840E+03_r8,     0.1905716E+03_r8,&  
         0.1893664E+03_r8,     0.1881684E+03_r8,     0.1869776E+03_r8,     0.1857941E+03_r8,     0.1846177E+03_r8,&  
         0.1834486E+03_r8,     0.1822867E+03_r8,     0.1811320E+03_r8,     0.1799846E+03_r8,     0.1788442E+03_r8,&  
         0.1777110E+03_r8,     0.1765849E+03_r8,     0.1754658E+03_r8,     0.1743539E+03_r8,     0.1732489E+03_r8,&  
         0.1721508E+03_r8,     0.1710597E+03_r8,     0.1699755E+03_r8,     0.1688981E+03_r8,     0.1678276E+03_r8,&  
         0.1667638E+03_r8,     0.1657067E+03_r8,     0.1646563E+03_r8,     0.1636126E+03_r8,     0.1625755E+03_r8/)  

    psaditmk(1:150,   88)= (/ &
         0.3094497E+03_r8,     0.3087765E+03_r8,     0.3081042E+03_r8,     0.3074310E+03_r8,     0.3067583E+03_r8,&  
         0.3060850E+03_r8,     0.3054122E+03_r8,     0.3047383E+03_r8,     0.3040648E+03_r8,     0.3033904E+03_r8,&  
         0.3027159E+03_r8,     0.3020405E+03_r8,     0.3013649E+03_r8,     0.3006886E+03_r8,     0.3000112E+03_r8,&  
         0.2993336E+03_r8,     0.2986547E+03_r8,     0.2979753E+03_r8,     0.2972946E+03_r8,     0.2966130E+03_r8,&  
         0.2959302E+03_r8,     0.2952462E+03_r8,     0.2945608E+03_r8,     0.2938742E+03_r8,     0.2931860E+03_r8,&  
         0.2924962E+03_r8,     0.2918050E+03_r8,     0.2911116E+03_r8,     0.2904164E+03_r8,     0.2897193E+03_r8,&  
         0.2890202E+03_r8,     0.2883186E+03_r8,     0.2876151E+03_r8,     0.2869088E+03_r8,     0.2862002E+03_r8,&  
         0.2854886E+03_r8,     0.2847744E+03_r8,     0.2840572E+03_r8,     0.2833367E+03_r8,     0.2826131E+03_r8,&  
         0.2818863E+03_r8,     0.2811555E+03_r8,     0.2804210E+03_r8,     0.2796826E+03_r8,     0.2789400E+03_r8,&  
         0.2781934E+03_r8,     0.2774421E+03_r8,     0.2766862E+03_r8,     0.2759254E+03_r8,     0.2751593E+03_r8,&  
         0.2743880E+03_r8,     0.2736109E+03_r8,     0.2728282E+03_r8,     0.2720394E+03_r8,     0.2712441E+03_r8,&  
         0.2704423E+03_r8,     0.2696336E+03_r8,     0.2688174E+03_r8,     0.2679940E+03_r8,     0.2671628E+03_r8,&  
         0.2663236E+03_r8,     0.2654757E+03_r8,     0.2646190E+03_r8,     0.2637534E+03_r8,     0.2628777E+03_r8,&  
         0.2619926E+03_r8,     0.2610973E+03_r8,     0.2601913E+03_r8,     0.2592744E+03_r8,     0.2583459E+03_r8,&  
         0.2574057E+03_r8,     0.2564537E+03_r8,     0.2554890E+03_r8,     0.2545113E+03_r8,     0.2535205E+03_r8,&  
         0.2525161E+03_r8,     0.2514978E+03_r8,     0.2504654E+03_r8,     0.2494186E+03_r8,     0.2483571E+03_r8,&  
         0.2472808E+03_r8,     0.2461897E+03_r8,     0.2450835E+03_r8,     0.2439623E+03_r8,     0.2428262E+03_r8,&  
         0.2416752E+03_r8,     0.2405096E+03_r8,     0.2393298E+03_r8,     0.2381359E+03_r8,     0.2369287E+03_r8,&  
         0.2357084E+03_r8,     0.2344758E+03_r8,     0.2332315E+03_r8,     0.2319763E+03_r8,     0.2307111E+03_r8,&  
         0.2294366E+03_r8,     0.2281539E+03_r8,     0.2268638E+03_r8,     0.2255673E+03_r8,     0.2242654E+03_r8,&  
         0.2229591E+03_r8,     0.2216492E+03_r8,     0.2203368E+03_r8,     0.2190228E+03_r8,     0.2177079E+03_r8,&  
         0.2163931E+03_r8,     0.2150791E+03_r8,     0.2137666E+03_r8,     0.2124563E+03_r8,     0.2111488E+03_r8,&  
         0.2098448E+03_r8,     0.2085446E+03_r8,     0.2072488E+03_r8,     0.2059578E+03_r8,     0.2046719E+03_r8,&  
         0.2033916E+03_r8,     0.2021171E+03_r8,     0.2008485E+03_r8,     0.1995863E+03_r8,     0.1983305E+03_r8,&  
         0.1970813E+03_r8,     0.1958389E+03_r8,     0.1946033E+03_r8,     0.1933748E+03_r8,     0.1921532E+03_r8,&  
         0.1909388E+03_r8,     0.1897314E+03_r8,     0.1885313E+03_r8,     0.1873384E+03_r8,     0.1861527E+03_r8,&  
         0.1849742E+03_r8,     0.1838030E+03_r8,     0.1826389E+03_r8,     0.1814821E+03_r8,     0.1803324E+03_r8,&  
         0.1791899E+03_r8,     0.1780546E+03_r8,     0.1769263E+03_r8,     0.1758051E+03_r8,     0.1746910E+03_r8,&  
         0.1735839E+03_r8,     0.1724837E+03_r8,     0.1713906E+03_r8,     0.1703042E+03_r8,     0.1692248E+03_r8,&  
         0.1681522E+03_r8,     0.1670863E+03_r8,     0.1660272E+03_r8,     0.1649748E+03_r8,     0.1639291E+03_r8/)  

    psaditmk(1:150,   89)= (/ &
         0.3098485E+03_r8,     0.3091789E+03_r8,     0.3085098E+03_r8,     0.3078398E+03_r8,     0.3071705E+03_r8,&  
         0.3065005E+03_r8,     0.3058311E+03_r8,     0.3051610E+03_r8,     0.3044910E+03_r8,     0.3038206E+03_r8,&  
         0.3031495E+03_r8,     0.3024784E+03_r8,     0.3018066E+03_r8,     0.3011346E+03_r8,     0.3004614E+03_r8,&  
         0.2997882E+03_r8,     0.2991136E+03_r8,     0.2984388E+03_r8,     0.2977626E+03_r8,     0.2970858E+03_r8,&  
         0.2964078E+03_r8,     0.2957289E+03_r8,     0.2950484E+03_r8,     0.2943672E+03_r8,     0.2936842E+03_r8,&  
         0.2930000E+03_r8,     0.2923144E+03_r8,     0.2916265E+03_r8,     0.2909376E+03_r8,     0.2902467E+03_r8,&  
         0.2895538E+03_r8,     0.2888588E+03_r8,     0.2881618E+03_r8,     0.2874625E+03_r8,     0.2867608E+03_r8,&  
         0.2860563E+03_r8,     0.2853496E+03_r8,     0.2846400E+03_r8,     0.2839276E+03_r8,     0.2832120E+03_r8,&  
         0.2824933E+03_r8,     0.2817711E+03_r8,     0.2810457E+03_r8,     0.2803165E+03_r8,     0.2795833E+03_r8,&  
         0.2788463E+03_r8,     0.2781050E+03_r8,     0.2773594E+03_r8,     0.2766091E+03_r8,     0.2758541E+03_r8,&  
         0.2750941E+03_r8,     0.2743286E+03_r8,     0.2735580E+03_r8,     0.2727816E+03_r8,     0.2719993E+03_r8,&  
         0.2712106E+03_r8,     0.2704156E+03_r8,     0.2696140E+03_r8,     0.2688048E+03_r8,     0.2679885E+03_r8,&  
         0.2671648E+03_r8,     0.2663332E+03_r8,     0.2654932E+03_r8,     0.2646446E+03_r8,     0.2637870E+03_r8,&  
         0.2629201E+03_r8,     0.2620436E+03_r8,     0.2611569E+03_r8,     0.2602598E+03_r8,     0.2593521E+03_r8,&  
         0.2584331E+03_r8,     0.2575025E+03_r8,     0.2565601E+03_r8,     0.2556052E+03_r8,     0.2546376E+03_r8,&  
         0.2536571E+03_r8,     0.2526633E+03_r8,     0.2516558E+03_r8,     0.2506341E+03_r8,     0.2495981E+03_r8,&  
         0.2485475E+03_r8,     0.2474822E+03_r8,     0.2464020E+03_r8,     0.2453068E+03_r8,     0.2441965E+03_r8,&  
         0.2430712E+03_r8,     0.2419310E+03_r8,     0.2407759E+03_r8,     0.2396062E+03_r8,     0.2384226E+03_r8,&  
         0.2372250E+03_r8,     0.2360140E+03_r8,     0.2347904E+03_r8,     0.2335546E+03_r8,     0.2323075E+03_r8,&  
         0.2310498E+03_r8,     0.2297824E+03_r8,     0.2285061E+03_r8,     0.2272220E+03_r8,     0.2259308E+03_r8,&  
         0.2246338E+03_r8,     0.2233316E+03_r8,     0.2220255E+03_r8,     0.2207163E+03_r8,     0.2194050E+03_r8,&  
         0.2180923E+03_r8,     0.2167793E+03_r8,     0.2154666E+03_r8,     0.2141550E+03_r8,     0.2128453E+03_r8,&  
         0.2115381E+03_r8,     0.2102339E+03_r8,     0.2089335E+03_r8,     0.2076371E+03_r8,     0.2063453E+03_r8,&  
         0.2050585E+03_r8,     0.2037769E+03_r8,     0.2025010E+03_r8,     0.2012311E+03_r8,     0.1999673E+03_r8,&  
         0.1987098E+03_r8,     0.1974588E+03_r8,     0.1962146E+03_r8,     0.1949771E+03_r8,     0.1937466E+03_r8,&  
         0.1925230E+03_r8,     0.1913065E+03_r8,     0.1900971E+03_r8,     0.1888949E+03_r8,     0.1876998E+03_r8,&  
         0.1865120E+03_r8,     0.1853314E+03_r8,     0.1841579E+03_r8,     0.1829918E+03_r8,     0.1818327E+03_r8,&  
         0.1806809E+03_r8,     0.1795363E+03_r8,     0.1783988E+03_r8,     0.1772684E+03_r8,     0.1761451E+03_r8,&  
         0.1750288E+03_r8,     0.1739196E+03_r8,     0.1728173E+03_r8,     0.1717220E+03_r8,     0.1706336E+03_r8,&  
         0.1695521E+03_r8,     0.1684774E+03_r8,     0.1674095E+03_r8,     0.1663484E+03_r8,     0.1652939E+03_r8/)  

    psaditmk(1:150,   90)= (/ &
         0.3102406E+03_r8,     0.3095744E+03_r8,     0.3089076E+03_r8,     0.3082416E+03_r8,     0.3075750E+03_r8,&  
         0.3069088E+03_r8,     0.3062424E+03_r8,     0.3055762E+03_r8,     0.3049091E+03_r8,     0.3042425E+03_r8,&  
         0.3035751E+03_r8,     0.3029078E+03_r8,     0.3022399E+03_r8,     0.3015717E+03_r8,     0.3009025E+03_r8,&  
         0.3002334E+03_r8,     0.2995633E+03_r8,     0.2988926E+03_r8,     0.2982210E+03_r8,     0.2975490E+03_r8,&  
         0.2968754E+03_r8,     0.2962010E+03_r8,     0.2955260E+03_r8,     0.2948492E+03_r8,     0.2941719E+03_r8,&  
         0.2934926E+03_r8,     0.2928126E+03_r8,     0.2921302E+03_r8,     0.2914474E+03_r8,     0.2907619E+03_r8,&  
         0.2900750E+03_r8,     0.2893863E+03_r8,     0.2886956E+03_r8,     0.2880028E+03_r8,     0.2873078E+03_r8,&  
         0.2866104E+03_r8,     0.2859106E+03_r8,     0.2852084E+03_r8,     0.2845034E+03_r8,     0.2837957E+03_r8,&  
         0.2830848E+03_r8,     0.2823710E+03_r8,     0.2816538E+03_r8,     0.2809333E+03_r8,     0.2802090E+03_r8,&  
         0.2794812E+03_r8,     0.2787494E+03_r8,     0.2780136E+03_r8,     0.2772734E+03_r8,     0.2765288E+03_r8,&  
         0.2757792E+03_r8,     0.2750252E+03_r8,     0.2742660E+03_r8,     0.2735015E+03_r8,     0.2727313E+03_r8,&  
         0.2719552E+03_r8,     0.2711733E+03_r8,     0.2703850E+03_r8,     0.2695901E+03_r8,     0.2687881E+03_r8,&  
         0.2679789E+03_r8,     0.2671627E+03_r8,     0.2663382E+03_r8,     0.2655062E+03_r8,     0.2646654E+03_r8,&  
         0.2638159E+03_r8,     0.2629574E+03_r8,     0.2620890E+03_r8,     0.2612112E+03_r8,     0.2603230E+03_r8,&  
         0.2594242E+03_r8,     0.2585144E+03_r8,     0.2575933E+03_r8,     0.2566607E+03_r8,     0.2557157E+03_r8,&  
         0.2547582E+03_r8,     0.2537878E+03_r8,     0.2528042E+03_r8,     0.2518071E+03_r8,     0.2507961E+03_r8,&  
         0.2497708E+03_r8,     0.2487314E+03_r8,     0.2476770E+03_r8,     0.2466078E+03_r8,     0.2455235E+03_r8,&  
         0.2444242E+03_r8,     0.2433098E+03_r8,     0.2421803E+03_r8,     0.2410358E+03_r8,     0.2398769E+03_r8,&  
         0.2387032E+03_r8,     0.2375155E+03_r8,     0.2363141E+03_r8,     0.2350996E+03_r8,     0.2338727E+03_r8,&  
         0.2326339E+03_r8,     0.2313840E+03_r8,     0.2301239E+03_r8,     0.2288544E+03_r8,     0.2275765E+03_r8,&  
         0.2262911E+03_r8,     0.2249990E+03_r8,     0.2237015E+03_r8,     0.2223994E+03_r8,     0.2210937E+03_r8,&  
         0.2197852E+03_r8,     0.2184751E+03_r8,     0.2171640E+03_r8,     0.2158529E+03_r8,     0.2145425E+03_r8,&  
         0.2132336E+03_r8,     0.2119268E+03_r8,     0.2106227E+03_r8,     0.2093221E+03_r8,     0.2080252E+03_r8,&  
         0.2067327E+03_r8,     0.2054450E+03_r8,     0.2041624E+03_r8,     0.2028853E+03_r8,     0.2016139E+03_r8,&  
         0.2003485E+03_r8,     0.1990895E+03_r8,     0.1978368E+03_r8,     0.1965907E+03_r8,     0.1953514E+03_r8,&  
         0.1941189E+03_r8,     0.1928933E+03_r8,     0.1916748E+03_r8,     0.1904633E+03_r8,     0.1892590E+03_r8,&  
         0.1880619E+03_r8,     0.1868719E+03_r8,     0.1856891E+03_r8,     0.1845136E+03_r8,     0.1833452E+03_r8,&  
         0.1821841E+03_r8,     0.1810301E+03_r8,     0.1798833E+03_r8,     0.1787436E+03_r8,     0.1776111E+03_r8,&  
         0.1764856E+03_r8,     0.1753673E+03_r8,     0.1742559E+03_r8,     0.1731515E+03_r8,     0.1720541E+03_r8,&  
         0.1709637E+03_r8,     0.1698800E+03_r8,     0.1688033E+03_r8,     0.1677333E+03_r8,     0.1666701E+03_r8/)  

    psaditmk(1:150,   91)= (/ &
         0.3106259E+03_r8,     0.3099627E+03_r8,     0.3092991E+03_r8,     0.3086358E+03_r8,     0.3079726E+03_r8,&  
         0.3073096E+03_r8,     0.3066460E+03_r8,     0.3059833E+03_r8,     0.3053198E+03_r8,     0.3046565E+03_r8,&  
         0.3039929E+03_r8,     0.3033293E+03_r8,     0.3026648E+03_r8,     0.3020003E+03_r8,     0.3013354E+03_r8,&  
         0.3006698E+03_r8,     0.3000041E+03_r8,     0.2993374E+03_r8,     0.2986702E+03_r8,     0.2980021E+03_r8,&  
         0.2973334E+03_r8,     0.2966637E+03_r8,     0.2959932E+03_r8,     0.2953214E+03_r8,     0.2946487E+03_r8,&  
         0.2939748E+03_r8,     0.2932996E+03_r8,     0.2926230E+03_r8,     0.2919451E+03_r8,     0.2912656E+03_r8,&  
         0.2905845E+03_r8,     0.2899016E+03_r8,     0.2892169E+03_r8,     0.2885303E+03_r8,     0.2878421E+03_r8,&  
         0.2871511E+03_r8,     0.2864580E+03_r8,     0.2857626E+03_r8,     0.2850649E+03_r8,     0.2843646E+03_r8,&  
         0.2836613E+03_r8,     0.2829551E+03_r8,     0.2822460E+03_r8,     0.2815337E+03_r8,     0.2808182E+03_r8,&  
         0.2800990E+03_r8,     0.2793761E+03_r8,     0.2786495E+03_r8,     0.2779190E+03_r8,     0.2771843E+03_r8,&  
         0.2764452E+03_r8,     0.2757017E+03_r8,     0.2749533E+03_r8,     0.2742000E+03_r8,     0.2734413E+03_r8,&  
         0.2726774E+03_r8,     0.2719077E+03_r8,     0.2711321E+03_r8,     0.2703505E+03_r8,     0.2695624E+03_r8,&  
         0.2687673E+03_r8,     0.2679653E+03_r8,     0.2671561E+03_r8,     0.2663392E+03_r8,     0.2655145E+03_r8,&  
         0.2646816E+03_r8,     0.2638399E+03_r8,     0.2629893E+03_r8,     0.2621296E+03_r8,     0.2612602E+03_r8,&  
         0.2603808E+03_r8,     0.2594909E+03_r8,     0.2585904E+03_r8,     0.2576786E+03_r8,     0.2567552E+03_r8,&  
         0.2558199E+03_r8,     0.2548725E+03_r8,     0.2539124E+03_r8,     0.2529391E+03_r8,     0.2519524E+03_r8,&  
         0.2509519E+03_r8,     0.2499374E+03_r8,     0.2489085E+03_r8,     0.2478650E+03_r8,     0.2468068E+03_r8,&  
         0.2457337E+03_r8,     0.2446453E+03_r8,     0.2435419E+03_r8,     0.2424232E+03_r8,     0.2412895E+03_r8,&  
         0.2401409E+03_r8,     0.2389777E+03_r8,     0.2378001E+03_r8,     0.2366085E+03_r8,     0.2354035E+03_r8,&  
         0.2341855E+03_r8,     0.2329553E+03_r8,     0.2317135E+03_r8,     0.2304610E+03_r8,     0.2291987E+03_r8,&  
         0.2279272E+03_r8,     0.2266477E+03_r8,     0.2253613E+03_r8,     0.2240685E+03_r8,     0.2227707E+03_r8,&  
         0.2214688E+03_r8,     0.2201636E+03_r8,     0.2188561E+03_r8,     0.2175473E+03_r8,     0.2162380E+03_r8,&  
         0.2149289E+03_r8,     0.2136209E+03_r8,     0.2123148E+03_r8,     0.2110110E+03_r8,     0.2097103E+03_r8,&  
         0.2084132E+03_r8,     0.2071201E+03_r8,     0.2058316E+03_r8,     0.2045480E+03_r8,     0.2032697E+03_r8,&  
         0.2019970E+03_r8,     0.2007302E+03_r8,     0.1994695E+03_r8,     0.1982151E+03_r8,     0.1969673E+03_r8,&  
         0.1957261E+03_r8,     0.1944917E+03_r8,     0.1932642E+03_r8,     0.1920436E+03_r8,     0.1908301E+03_r8,&  
         0.1896238E+03_r8,     0.1884245E+03_r8,     0.1872324E+03_r8,     0.1860475E+03_r8,     0.1848698E+03_r8,&  
         0.1836993E+03_r8,     0.1825360E+03_r8,     0.1813799E+03_r8,     0.1802309E+03_r8,     0.1790891E+03_r8,&  
         0.1779544E+03_r8,     0.1768269E+03_r8,     0.1757063E+03_r8,     0.1745928E+03_r8,     0.1734864E+03_r8,&  
         0.1723869E+03_r8,     0.1712943E+03_r8,     0.1702086E+03_r8,     0.1691298E+03_r8,     0.1680578E+03_r8/)  

    psaditmk(1:150,   92)= (/ &
         0.3110053E+03_r8,     0.3103439E+03_r8,     0.3096838E+03_r8,     0.3090231E+03_r8,     0.3083634E+03_r8,&  
         0.3077029E+03_r8,     0.3070431E+03_r8,     0.3063829E+03_r8,     0.3057230E+03_r8,     0.3050628E+03_r8,&  
         0.3044029E+03_r8,     0.3037424E+03_r8,     0.3030819E+03_r8,     0.3024209E+03_r8,     0.3017597E+03_r8,&  
         0.3010978E+03_r8,     0.3004361E+03_r8,     0.2997735E+03_r8,     0.2991106E+03_r8,     0.2984465E+03_r8,&  
         0.2977822E+03_r8,     0.2971168E+03_r8,     0.2964508E+03_r8,     0.2957836E+03_r8,     0.2951156E+03_r8,&  
         0.2944466E+03_r8,     0.2937764E+03_r8,     0.2931049E+03_r8,     0.2924321E+03_r8,     0.2917581E+03_r8,&  
         0.2910825E+03_r8,     0.2904052E+03_r8,     0.2897263E+03_r8,     0.2890457E+03_r8,     0.2883635E+03_r8,&  
         0.2876791E+03_r8,     0.2869923E+03_r8,     0.2863036E+03_r8,     0.2856126E+03_r8,     0.2849189E+03_r8,&  
         0.2842232E+03_r8,     0.2835243E+03_r8,     0.2828229E+03_r8,     0.2821186E+03_r8,     0.2814111E+03_r8,&  
         0.2807004E+03_r8,     0.2799863E+03_r8,     0.2792686E+03_r8,     0.2785471E+03_r8,     0.2778219E+03_r8,&  
         0.2770926E+03_r8,     0.2763587E+03_r8,     0.2756208E+03_r8,     0.2748782E+03_r8,     0.2741306E+03_r8,&  
         0.2733780E+03_r8,     0.2726200E+03_r8,     0.2718565E+03_r8,     0.2710873E+03_r8,     0.2703121E+03_r8,&  
         0.2695305E+03_r8,     0.2687421E+03_r8,     0.2679473E+03_r8,     0.2671452E+03_r8,     0.2663356E+03_r8,&  
         0.2655183E+03_r8,     0.2646929E+03_r8,     0.2638592E+03_r8,     0.2630167E+03_r8,     0.2621653E+03_r8,&  
         0.2613041E+03_r8,     0.2604333E+03_r8,     0.2595522E+03_r8,     0.2586606E+03_r8,     0.2577580E+03_r8,&  
         0.2568443E+03_r8,     0.2559186E+03_r8,     0.2549808E+03_r8,     0.2540305E+03_r8,     0.2530673E+03_r8,&  
         0.2520911E+03_r8,     0.2511014E+03_r8,     0.2500974E+03_r8,     0.2490793E+03_r8,     0.2480467E+03_r8,&  
         0.2469993E+03_r8,     0.2459371E+03_r8,     0.2448598E+03_r8,     0.2437673E+03_r8,     0.2426597E+03_r8,&  
         0.2415367E+03_r8,     0.2403988E+03_r8,     0.2392460E+03_r8,     0.2380785E+03_r8,     0.2368970E+03_r8,&  
         0.2357016E+03_r8,     0.2344929E+03_r8,     0.2332715E+03_r8,     0.2320381E+03_r8,     0.2307935E+03_r8,&  
         0.2295385E+03_r8,     0.2282739E+03_r8,     0.2270008E+03_r8,     0.2257200E+03_r8,     0.2244325E+03_r8,&  
         0.2231392E+03_r8,     0.2218414E+03_r8,     0.2205397E+03_r8,     0.2192352E+03_r8,     0.2179289E+03_r8,&  
         0.2166216E+03_r8,     0.2153141E+03_r8,     0.2140074E+03_r8,     0.2127020E+03_r8,     0.2113986E+03_r8,&  
         0.2100981E+03_r8,     0.2088007E+03_r8,     0.2075072E+03_r8,     0.2062180E+03_r8,     0.2049336E+03_r8,&  
         0.2036542E+03_r8,     0.2023802E+03_r8,     0.2011120E+03_r8,     0.1998499E+03_r8,     0.1985939E+03_r8,&  
         0.1973443E+03_r8,     0.1961013E+03_r8,     0.1948650E+03_r8,     0.1936356E+03_r8,     0.1924131E+03_r8,&  
         0.1911975E+03_r8,     0.1899891E+03_r8,     0.1887877E+03_r8,     0.1875936E+03_r8,     0.1864065E+03_r8,&  
         0.1852267E+03_r8,     0.1840541E+03_r8,     0.1828886E+03_r8,     0.1817303E+03_r8,     0.1805792E+03_r8,&  
         0.1794353E+03_r8,     0.1782985E+03_r8,     0.1771687E+03_r8,     0.1760461E+03_r8,     0.1749305E+03_r8,&  
         0.1738219E+03_r8,     0.1727202E+03_r8,     0.1716256E+03_r8,     0.1705378E+03_r8,     0.1694569E+03_r8/)  

    psaditmk(1:150,   93)= (/ &
         0.3113774E+03_r8,     0.3107196E+03_r8,     0.3100614E+03_r8,     0.3094043E+03_r8,     0.3087465E+03_r8,&  
         0.3080897E+03_r8,     0.3074327E+03_r8,     0.3067760E+03_r8,     0.3061193E+03_r8,     0.3054621E+03_r8,&  
         0.3048053E+03_r8,     0.3041483E+03_r8,     0.3034911E+03_r8,     0.3028333E+03_r8,     0.3021763E+03_r8,&  
         0.3015179E+03_r8,     0.3008599E+03_r8,     0.3002011E+03_r8,     0.2995420E+03_r8,     0.2988820E+03_r8,&  
         0.2982218E+03_r8,     0.2975608E+03_r8,     0.2968992E+03_r8,     0.2962363E+03_r8,     0.2955730E+03_r8,&  
         0.2949084E+03_r8,     0.2942430E+03_r8,     0.2935764E+03_r8,     0.2929087E+03_r8,     0.2922396E+03_r8,&  
         0.2915694E+03_r8,     0.2908976E+03_r8,     0.2902243E+03_r8,     0.2895494E+03_r8,     0.2888727E+03_r8,&  
         0.2881942E+03_r8,     0.2875138E+03_r8,     0.2868316E+03_r8,     0.2861470E+03_r8,     0.2854602E+03_r8,&  
         0.2847713E+03_r8,     0.2840799E+03_r8,     0.2833856E+03_r8,     0.2826887E+03_r8,     0.2819888E+03_r8,&  
         0.2812863E+03_r8,     0.2805802E+03_r8,     0.2798709E+03_r8,     0.2791583E+03_r8,     0.2784420E+03_r8,&  
         0.2777219E+03_r8,     0.2769977E+03_r8,     0.2762695E+03_r8,     0.2755369E+03_r8,     0.2747999E+03_r8,&  
         0.2740579E+03_r8,     0.2733112E+03_r8,     0.2725591E+03_r8,     0.2718018E+03_r8,     0.2710388E+03_r8,&  
         0.2702700E+03_r8,     0.2694948E+03_r8,     0.2687132E+03_r8,     0.2679253E+03_r8,     0.2671301E+03_r8,&  
         0.2663277E+03_r8,     0.2655177E+03_r8,     0.2646999E+03_r8,     0.2638737E+03_r8,     0.2630393E+03_r8,&  
         0.2621957E+03_r8,     0.2613431E+03_r8,     0.2604806E+03_r8,     0.2596082E+03_r8,     0.2587255E+03_r8,&  
         0.2578321E+03_r8,     0.2569274E+03_r8,     0.2560111E+03_r8,     0.2550834E+03_r8,     0.2541430E+03_r8,&  
         0.2531898E+03_r8,     0.2522237E+03_r8,     0.2512440E+03_r8,     0.2502507E+03_r8,     0.2492436E+03_r8,&  
         0.2482218E+03_r8,     0.2471853E+03_r8,     0.2461340E+03_r8,     0.2450676E+03_r8,     0.2439862E+03_r8,&  
         0.2428894E+03_r8,     0.2417775E+03_r8,     0.2406502E+03_r8,     0.2395080E+03_r8,     0.2383511E+03_r8,&  
         0.2371796E+03_r8,     0.2359940E+03_r8,     0.2347947E+03_r8,     0.2335824E+03_r8,     0.2323578E+03_r8,&  
         0.2311213E+03_r8,     0.2298739E+03_r8,     0.2286165E+03_r8,     0.2273500E+03_r8,     0.2260752E+03_r8,&  
         0.2247932E+03_r8,     0.2235049E+03_r8,     0.2222113E+03_r8,     0.2209135E+03_r8,     0.2196123E+03_r8,&  
         0.2183088E+03_r8,     0.2170037E+03_r8,     0.2156981E+03_r8,     0.2143927E+03_r8,     0.2130883E+03_r8,&  
         0.2117856E+03_r8,     0.2104853E+03_r8,     0.2091879E+03_r8,     0.2078941E+03_r8,     0.2066044E+03_r8,&  
         0.2053191E+03_r8,     0.2040388E+03_r8,     0.2027637E+03_r8,     0.2014942E+03_r8,     0.2002305E+03_r8,&  
         0.1989730E+03_r8,     0.1977217E+03_r8,     0.1964770E+03_r8,     0.1952389E+03_r8,     0.1940075E+03_r8,&  
         0.1927830E+03_r8,     0.1915655E+03_r8,     0.1903550E+03_r8,     0.1891516E+03_r8,     0.1879553E+03_r8,&  
         0.1867662E+03_r8,     0.1855842E+03_r8,     0.1844094E+03_r8,     0.1832418E+03_r8,     0.1820814E+03_r8,&  
         0.1809282E+03_r8,     0.1797821E+03_r8,     0.1786431E+03_r8,     0.1775112E+03_r8,     0.1763864E+03_r8,&  
         0.1752687E+03_r8,     0.1741580E+03_r8,     0.1730543E+03_r8,     0.1719575E+03_r8,     0.1708676E+03_r8/)  

    psaditmk(1:150,   94)= (/ &
         0.3117440E+03_r8,     0.3110881E+03_r8,     0.3104333E+03_r8,     0.3097787E+03_r8,     0.3091238E+03_r8,&  
         0.3084697E+03_r8,     0.3078155E+03_r8,     0.3071615E+03_r8,     0.3065078E+03_r8,     0.3058542E+03_r8,&  
         0.3052001E+03_r8,     0.3045467E+03_r8,     0.3038926E+03_r8,     0.3032388E+03_r8,     0.3025845E+03_r8,&  
         0.3019302E+03_r8,     0.3012754E+03_r8,     0.3006205E+03_r8,     0.2999649E+03_r8,     0.2993093E+03_r8,&  
         0.2986527E+03_r8,     0.2979958E+03_r8,     0.2973382E+03_r8,     0.2966800E+03_r8,     0.2960210E+03_r8,&  
         0.2953604E+03_r8,     0.2947000E+03_r8,     0.2940380E+03_r8,     0.2933752E+03_r8,     0.2927109E+03_r8,&  
         0.2920457E+03_r8,     0.2913790E+03_r8,     0.2907111E+03_r8,     0.2900416E+03_r8,     0.2893707E+03_r8,&  
         0.2886978E+03_r8,     0.2880236E+03_r8,     0.2873472E+03_r8,     0.2866689E+03_r8,     0.2859886E+03_r8,&  
         0.2853062E+03_r8,     0.2846214E+03_r8,     0.2839342E+03_r8,     0.2832444E+03_r8,     0.2825517E+03_r8,&  
         0.2818567E+03_r8,     0.2811586E+03_r8,     0.2804576E+03_r8,     0.2797531E+03_r8,     0.2790452E+03_r8,&  
         0.2783337E+03_r8,     0.2776189E+03_r8,     0.2769001E+03_r8,     0.2761771E+03_r8,     0.2754500E+03_r8,&  
         0.2747184E+03_r8,     0.2739819E+03_r8,     0.2732410E+03_r8,     0.2724950E+03_r8,     0.2717436E+03_r8,&  
         0.2709868E+03_r8,     0.2702241E+03_r8,     0.2694553E+03_r8,     0.2686803E+03_r8,     0.2678991E+03_r8,&  
         0.2671109E+03_r8,     0.2663156E+03_r8,     0.2655128E+03_r8,     0.2647024E+03_r8,     0.2638839E+03_r8,&  
         0.2630570E+03_r8,     0.2622213E+03_r8,     0.2613770E+03_r8,     0.2605229E+03_r8,     0.2596591E+03_r8,&  
         0.2587851E+03_r8,     0.2579005E+03_r8,     0.2570052E+03_r8,     0.2560984E+03_r8,     0.2551798E+03_r8,&  
         0.2542490E+03_r8,     0.2533061E+03_r8,     0.2523503E+03_r8,     0.2513809E+03_r8,     0.2503980E+03_r8,&  
         0.2494012E+03_r8,     0.2483901E+03_r8,     0.2473647E+03_r8,     0.2463243E+03_r8,     0.2452689E+03_r8,&  
         0.2441984E+03_r8,     0.2431127E+03_r8,     0.2420116E+03_r8,     0.2408954E+03_r8,     0.2397638E+03_r8,&  
         0.2386173E+03_r8,     0.2374561E+03_r8,     0.2362805E+03_r8,     0.2350909E+03_r8,     0.2338879E+03_r8,&  
         0.2326721E+03_r8,     0.2314441E+03_r8,     0.2302048E+03_r8,     0.2289548E+03_r8,     0.2276952E+03_r8,&  
         0.2264268E+03_r8,     0.2251505E+03_r8,     0.2238675E+03_r8,     0.2225785E+03_r8,     0.2212848E+03_r8,&  
         0.2199872E+03_r8,     0.2186867E+03_r8,     0.2173842E+03_r8,     0.2160806E+03_r8,     0.2147769E+03_r8,&  
         0.2134737E+03_r8,     0.2121718E+03_r8,     0.2108719E+03_r8,     0.2095747E+03_r8,     0.2082807E+03_r8,&  
         0.2069906E+03_r8,     0.2057047E+03_r8,     0.2044234E+03_r8,     0.2031472E+03_r8,     0.2018765E+03_r8,&  
         0.2006115E+03_r8,     0.1993524E+03_r8,     0.1980995E+03_r8,     0.1968530E+03_r8,     0.1956131E+03_r8,&  
         0.1943799E+03_r8,     0.1931535E+03_r8,     0.1919340E+03_r8,     0.1907215E+03_r8,     0.1895160E+03_r8,&  
         0.1883177E+03_r8,     0.1871264E+03_r8,     0.1859424E+03_r8,     0.1847654E+03_r8,     0.1835957E+03_r8,&  
         0.1824331E+03_r8,     0.1812777E+03_r8,     0.1801295E+03_r8,     0.1789884E+03_r8,     0.1778544E+03_r8,&  
         0.1767275E+03_r8,     0.1756076E+03_r8,     0.1744948E+03_r8,     0.1733889E+03_r8,     0.1722900E+03_r8/)  

    psaditmk(1:150,   95)= (/ &
         0.3121037E+03_r8,     0.3114512E+03_r8,     0.3107985E+03_r8,     0.3101461E+03_r8,     0.3094948E+03_r8,&  
         0.3088430E+03_r8,     0.3081920E+03_r8,     0.3075406E+03_r8,     0.3068899E+03_r8,     0.3062391E+03_r8,&  
         0.3055885E+03_r8,     0.3049375E+03_r8,     0.3042870E+03_r8,     0.3036365E+03_r8,     0.3029855E+03_r8,&  
         0.3023349E+03_r8,     0.3016833E+03_r8,     0.3010320E+03_r8,     0.3003799E+03_r8,     0.2997280E+03_r8,&  
         0.2990753E+03_r8,     0.2984222E+03_r8,     0.2977685E+03_r8,     0.2971144E+03_r8,     0.2964595E+03_r8,&  
         0.2958038E+03_r8,     0.2951473E+03_r8,     0.2944900E+03_r8,     0.2938313E+03_r8,     0.2931720E+03_r8,&  
         0.2925120E+03_r8,     0.2918502E+03_r8,     0.2911874E+03_r8,     0.2905230E+03_r8,     0.2898575E+03_r8,&  
         0.2891901E+03_r8,     0.2885214E+03_r8,     0.2878508E+03_r8,     0.2871786E+03_r8,     0.2865045E+03_r8,&  
         0.2858279E+03_r8,     0.2851498E+03_r8,     0.2844693E+03_r8,     0.2837866E+03_r8,     0.2831012E+03_r8,&  
         0.2824131E+03_r8,     0.2817224E+03_r8,     0.2810289E+03_r8,     0.2803324E+03_r8,     0.2796328E+03_r8,&  
         0.2789297E+03_r8,     0.2782234E+03_r8,     0.2775134E+03_r8,     0.2767995E+03_r8,     0.2760818E+03_r8,&  
         0.2753599E+03_r8,     0.2746337E+03_r8,     0.2739031E+03_r8,     0.2731678E+03_r8,     0.2724274E+03_r8,&  
         0.2716820E+03_r8,     0.2709312E+03_r8,     0.2701747E+03_r8,     0.2694123E+03_r8,     0.2686437E+03_r8,&  
         0.2678692E+03_r8,     0.2670877E+03_r8,     0.2662993E+03_r8,     0.2655037E+03_r8,     0.2647005E+03_r8,&  
         0.2638895E+03_r8,     0.2630703E+03_r8,     0.2622426E+03_r8,     0.2614059E+03_r8,     0.2605602E+03_r8,&  
         0.2597046E+03_r8,     0.2588393E+03_r8,     0.2579637E+03_r8,     0.2570772E+03_r8,     0.2561795E+03_r8,&  
         0.2552708E+03_r8,     0.2543497E+03_r8,     0.2534163E+03_r8,     0.2524702E+03_r8,     0.2515114E+03_r8,&  
         0.2505391E+03_r8,     0.2495527E+03_r8,     0.2485523E+03_r8,     0.2475375E+03_r8,     0.2465079E+03_r8,&  
         0.2454637E+03_r8,     0.2444042E+03_r8,     0.2433294E+03_r8,     0.2422393E+03_r8,     0.2411339E+03_r8,&  
         0.2400133E+03_r8,     0.2388773E+03_r8,     0.2377265E+03_r8,     0.2365611E+03_r8,     0.2353814E+03_r8,&  
         0.2341880E+03_r8,     0.2329812E+03_r8,     0.2317620E+03_r8,     0.2305309E+03_r8,     0.2292887E+03_r8,&  
         0.2280362E+03_r8,     0.2267745E+03_r8,     0.2255043E+03_r8,     0.2242268E+03_r8,     0.2229429E+03_r8,&  
         0.2216535E+03_r8,     0.2203598E+03_r8,     0.2190626E+03_r8,     0.2177629E+03_r8,     0.2164617E+03_r8,&  
         0.2151598E+03_r8,     0.2138579E+03_r8,     0.2125571E+03_r8,     0.2112578E+03_r8,     0.2099609E+03_r8,&  
         0.2086670E+03_r8,     0.2073765E+03_r8,     0.2060900E+03_r8,     0.2048080E+03_r8,     0.2035309E+03_r8,&  
         0.2022590E+03_r8,     0.2009926E+03_r8,     0.1997321E+03_r8,     0.1984777E+03_r8,     0.1972295E+03_r8,&  
         0.1959879E+03_r8,     0.1947528E+03_r8,     0.1935245E+03_r8,     0.1923031E+03_r8,     0.1910885E+03_r8,&  
         0.1898810E+03_r8,     0.1886806E+03_r8,     0.1874873E+03_r8,     0.1863011E+03_r8,     0.1851221E+03_r8,&  
         0.1839502E+03_r8,     0.1827855E+03_r8,     0.1816280E+03_r8,     0.1804776E+03_r8,     0.1793343E+03_r8,&  
         0.1781982E+03_r8,     0.1770691E+03_r8,     0.1759471E+03_r8,     0.1748322E+03_r8,     0.1737242E+03_r8/)  

    psaditmk(1:150,   96)= (/ &
         0.3124586E+03_r8,     0.3118080E+03_r8,     0.3111579E+03_r8,     0.3105085E+03_r8,     0.3098590E+03_r8,&  
         0.3092102E+03_r8,     0.3085615E+03_r8,     0.3079134E+03_r8,     0.3072655E+03_r8,     0.3066171E+03_r8,&  
         0.3059699E+03_r8,     0.3053219E+03_r8,     0.3046746E+03_r8,     0.3040268E+03_r8,     0.3033795E+03_r8,&  
         0.3027314E+03_r8,     0.3020839E+03_r8,     0.3014356E+03_r8,     0.3007878E+03_r8,     0.3001389E+03_r8,&  
         0.2994899E+03_r8,     0.2988406E+03_r8,     0.2981907E+03_r8,     0.2975405E+03_r8,     0.2968895E+03_r8,&  
         0.2962381E+03_r8,     0.2955857E+03_r8,     0.2949327E+03_r8,     0.2942787E+03_r8,     0.2936238E+03_r8,&  
         0.2929681E+03_r8,     0.2923113E+03_r8,     0.2916534E+03_r8,     0.2909940E+03_r8,     0.2903336E+03_r8,&  
         0.2896715E+03_r8,     0.2890083E+03_r8,     0.2883431E+03_r8,     0.2876766E+03_r8,     0.2870083E+03_r8,&  
         0.2863380E+03_r8,     0.2856660E+03_r8,     0.2849919E+03_r8,     0.2843154E+03_r8,     0.2836364E+03_r8,&  
         0.2829557E+03_r8,     0.2822717E+03_r8,     0.2815858E+03_r8,     0.2808969E+03_r8,     0.2802049E+03_r8,&  
         0.2795098E+03_r8,     0.2788117E+03_r8,     0.2781102E+03_r8,     0.2774052E+03_r8,     0.2766964E+03_r8,&  
         0.2759838E+03_r8,     0.2752672E+03_r8,     0.2745464E+03_r8,     0.2738212E+03_r8,     0.2730913E+03_r8,&  
         0.2723570E+03_r8,     0.2716171E+03_r8,     0.2708722E+03_r8,     0.2701218E+03_r8,     0.2693656E+03_r8,&  
         0.2686034E+03_r8,     0.2678354E+03_r8,     0.2670606E+03_r8,     0.2662790E+03_r8,     0.2654903E+03_r8,&  
         0.2646943E+03_r8,     0.2638908E+03_r8,     0.2630790E+03_r8,     0.2622589E+03_r8,     0.2614301E+03_r8,&  
         0.2605925E+03_r8,     0.2597454E+03_r8,     0.2588885E+03_r8,     0.2580215E+03_r8,     0.2571438E+03_r8,&  
         0.2562556E+03_r8,     0.2553557E+03_r8,     0.2544443E+03_r8,     0.2535207E+03_r8,     0.2525848E+03_r8,&  
         0.2516357E+03_r8,     0.2506733E+03_r8,     0.2496977E+03_r8,     0.2487080E+03_r8,     0.2477039E+03_r8,&  
         0.2466852E+03_r8,     0.2456517E+03_r8,     0.2446033E+03_r8,     0.2435395E+03_r8,     0.2424604E+03_r8,&  
         0.2413660E+03_r8,     0.2402560E+03_r8,     0.2391310E+03_r8,     0.2379908E+03_r8,     0.2368356E+03_r8,&  
         0.2356660E+03_r8,     0.2344823E+03_r8,     0.2332850E+03_r8,     0.2320747E+03_r8,     0.2308521E+03_r8,&  
         0.2296178E+03_r8,     0.2283729E+03_r8,     0.2271181E+03_r8,     0.2258544E+03_r8,     0.2245827E+03_r8,&  
         0.2233041E+03_r8,     0.2220195E+03_r8,     0.2207299E+03_r8,     0.2194362E+03_r8,     0.2181396E+03_r8,&  
         0.2168410E+03_r8,     0.2155412E+03_r8,     0.2142410E+03_r8,     0.2129414E+03_r8,     0.2116430E+03_r8,&  
         0.2103465E+03_r8,     0.2090528E+03_r8,     0.2077621E+03_r8,     0.2064752E+03_r8,     0.2051926E+03_r8,&  
         0.2039146E+03_r8,     0.2026416E+03_r8,     0.2013741E+03_r8,     0.2001122E+03_r8,     0.1988562E+03_r8,&  
         0.1976065E+03_r8,     0.1963631E+03_r8,     0.1951262E+03_r8,     0.1938960E+03_r8,     0.1926727E+03_r8,&  
         0.1914562E+03_r8,     0.1902467E+03_r8,     0.1890442E+03_r8,     0.1878488E+03_r8,     0.1866605E+03_r8,&  
         0.1854793E+03_r8,     0.1843053E+03_r8,     0.1831385E+03_r8,     0.1819788E+03_r8,     0.1808263E+03_r8,&  
         0.1796809E+03_r8,     0.1785426E+03_r8,     0.1774114E+03_r8,     0.1762873E+03_r8,     0.1751702E+03_r8/)  

    psaditmk(1:150,   97)= (/ &
         0.3128071E+03_r8,     0.3121589E+03_r8,     0.3115114E+03_r8,     0.3108642E+03_r8,     0.3102178E+03_r8,&  
         0.3095710E+03_r8,     0.3089251E+03_r8,     0.3082799E+03_r8,     0.3076342E+03_r8,     0.3069896E+03_r8,&  
         0.3063443E+03_r8,     0.3056997E+03_r8,     0.3050549E+03_r8,     0.3044106E+03_r8,     0.3037659E+03_r8,&  
         0.3031212E+03_r8,     0.3024767E+03_r8,     0.3018320E+03_r8,     0.3011873E+03_r8,     0.3005420E+03_r8,&  
         0.2998965E+03_r8,     0.2992505E+03_r8,     0.2986049E+03_r8,     0.2979580E+03_r8,     0.2973111E+03_r8,&  
         0.2966635E+03_r8,     0.2960150E+03_r8,     0.2953665E+03_r8,     0.2947167E+03_r8,     0.2940664E+03_r8,&  
         0.2934148E+03_r8,     0.2927626E+03_r8,     0.2921093E+03_r8,     0.2914551E+03_r8,     0.2907995E+03_r8,&  
         0.2901423E+03_r8,     0.2894843E+03_r8,     0.2888244E+03_r8,     0.2881634E+03_r8,     0.2875005E+03_r8,&  
         0.2868360E+03_r8,     0.2861699E+03_r8,     0.2855016E+03_r8,     0.2848315E+03_r8,     0.2841594E+03_r8,&  
         0.2834851E+03_r8,     0.2828081E+03_r8,     0.2821289E+03_r8,     0.2814473E+03_r8,     0.2807625E+03_r8,&  
         0.2800750E+03_r8,     0.2793848E+03_r8,     0.2786912E+03_r8,     0.2779944E+03_r8,     0.2772940E+03_r8,&  
         0.2765905E+03_r8,     0.2758830E+03_r8,     0.2751715E+03_r8,     0.2744560E+03_r8,     0.2737361E+03_r8,&  
         0.2730118E+03_r8,     0.2722828E+03_r8,     0.2715489E+03_r8,     0.2708099E+03_r8,     0.2700655E+03_r8,&  
         0.2693154E+03_r8,     0.2685595E+03_r8,     0.2677978E+03_r8,     0.2670295E+03_r8,     0.2662547E+03_r8,&  
         0.2654729E+03_r8,     0.2646840E+03_r8,     0.2638877E+03_r8,     0.2630832E+03_r8,     0.2622709E+03_r8,&  
         0.2614498E+03_r8,     0.2606200E+03_r8,     0.2597811E+03_r8,     0.2589326E+03_r8,     0.2580739E+03_r8,&  
         0.2572054E+03_r8,     0.2563258E+03_r8,     0.2554353E+03_r8,     0.2545333E+03_r8,     0.2536191E+03_r8,&  
         0.2526929E+03_r8,     0.2517540E+03_r8,     0.2508020E+03_r8,     0.2498364E+03_r8,     0.2488571E+03_r8,&  
         0.2478639E+03_r8,     0.2468560E+03_r8,     0.2458333E+03_r8,     0.2447957E+03_r8,     0.2437431E+03_r8,&  
         0.2426750E+03_r8,     0.2415915E+03_r8,     0.2404926E+03_r8,     0.2393783E+03_r8,     0.2382486E+03_r8,&  
         0.2371041E+03_r8,     0.2359446E+03_r8,     0.2347709E+03_r8,     0.2335831E+03_r8,     0.2323820E+03_r8,&  
         0.2311681E+03_r8,     0.2299423E+03_r8,     0.2287052E+03_r8,     0.2274576E+03_r8,     0.2262006E+03_r8,&  
         0.2249351E+03_r8,     0.2236620E+03_r8,     0.2223825E+03_r8,     0.2210973E+03_r8,     0.2198076E+03_r8,&  
         0.2185144E+03_r8,     0.2172186E+03_r8,     0.2159210E+03_r8,     0.2146228E+03_r8,     0.2133246E+03_r8,&  
         0.2120272E+03_r8,     0.2107315E+03_r8,     0.2094380E+03_r8,     0.2081474E+03_r8,     0.2068602E+03_r8,&  
         0.2055770E+03_r8,     0.2042982E+03_r8,     0.2030243E+03_r8,     0.2017556E+03_r8,     0.2004924E+03_r8,&  
         0.1992350E+03_r8,     0.1979837E+03_r8,     0.1967386E+03_r8,     0.1955000E+03_r8,     0.1942680E+03_r8,&  
         0.1930427E+03_r8,     0.1918243E+03_r8,     0.1906128E+03_r8,     0.1894083E+03_r8,     0.1882108E+03_r8,&  
         0.1870204E+03_r8,     0.1858372E+03_r8,     0.1846611E+03_r8,     0.1834921E+03_r8,     0.1823303E+03_r8,&  
         0.1811756E+03_r8,     0.1800281E+03_r8,     0.1788877E+03_r8,     0.1777543E+03_r8,     0.1766281E+03_r8/)  

    psaditmk(1:150,   98)= (/ &
         0.3131504E+03_r8,     0.3125046E+03_r8,     0.3118590E+03_r8,     0.3112147E+03_r8,     0.3105700E+03_r8,&  
         0.3099263E+03_r8,     0.3092831E+03_r8,     0.3086396E+03_r8,     0.3079973E+03_r8,     0.3073546E+03_r8,&  
         0.3067126E+03_r8,     0.3060706E+03_r8,     0.3054288E+03_r8,     0.3047875E+03_r8,     0.3041456E+03_r8,&  
         0.3035043E+03_r8,     0.3028625E+03_r8,     0.3022212E+03_r8,     0.3015793E+03_r8,     0.3009377E+03_r8,&  
         0.3002956E+03_r8,     0.2996530E+03_r8,     0.2990109E+03_r8,     0.2983677E+03_r8,     0.2977245E+03_r8,&  
         0.2970807E+03_r8,     0.2964365E+03_r8,     0.2957914E+03_r8,     0.2951460E+03_r8,     0.2944995E+03_r8,&  
         0.2938525E+03_r8,     0.2932048E+03_r8,     0.2925558E+03_r8,     0.2919060E+03_r8,     0.2912550E+03_r8,&  
         0.2906030E+03_r8,     0.2899500E+03_r8,     0.2892951E+03_r8,     0.2886395E+03_r8,     0.2879818E+03_r8,&  
         0.2873231E+03_r8,     0.2866624E+03_r8,     0.2859998E+03_r8,     0.2853359E+03_r8,     0.2846699E+03_r8,&  
         0.2840015E+03_r8,     0.2833312E+03_r8,     0.2826585E+03_r8,     0.2819833E+03_r8,     0.2813060E+03_r8,&  
         0.2806261E+03_r8,     0.2799431E+03_r8,     0.2792573E+03_r8,     0.2785683E+03_r8,     0.2778763E+03_r8,&  
         0.2771809E+03_r8,     0.2764820E+03_r8,     0.2757794E+03_r8,     0.2750730E+03_r8,     0.2743626E+03_r8,&  
         0.2736482E+03_r8,     0.2729292E+03_r8,     0.2722058E+03_r8,     0.2714776E+03_r8,     0.2707443E+03_r8,&  
         0.2700058E+03_r8,     0.2692617E+03_r8,     0.2685121E+03_r8,     0.2677566E+03_r8,     0.2669948E+03_r8,&  
         0.2662265E+03_r8,     0.2654515E+03_r8,     0.2646694E+03_r8,     0.2638802E+03_r8,     0.2630831E+03_r8,&  
         0.2622781E+03_r8,     0.2614651E+03_r8,     0.2606429E+03_r8,     0.2598118E+03_r8,     0.2589716E+03_r8,&  
         0.2581216E+03_r8,     0.2572614E+03_r8,     0.2563911E+03_r8,     0.2555093E+03_r8,     0.2546165E+03_r8,&  
         0.2537119E+03_r8,     0.2527954E+03_r8,     0.2518661E+03_r8,     0.2509243E+03_r8,     0.2499692E+03_r8,&  
         0.2490002E+03_r8,     0.2480173E+03_r8,     0.2470203E+03_r8,     0.2460085E+03_r8,     0.2449818E+03_r8,&  
         0.2439400E+03_r8,     0.2428830E+03_r8,     0.2418104E+03_r8,     0.2407225E+03_r8,     0.2396191E+03_r8,&  
         0.2385003E+03_r8,     0.2373663E+03_r8,     0.2362171E+03_r8,     0.2350536E+03_r8,     0.2338756E+03_r8,&  
         0.2326840E+03_r8,     0.2314791E+03_r8,     0.2302617E+03_r8,     0.2290327E+03_r8,     0.2277927E+03_r8,&  
         0.2265428E+03_r8,     0.2252838E+03_r8,     0.2240166E+03_r8,     0.2227423E+03_r8,     0.2214619E+03_r8,&  
         0.2201764E+03_r8,     0.2188868E+03_r8,     0.2175941E+03_r8,     0.2162993E+03_r8,     0.2150031E+03_r8,&  
         0.2137066E+03_r8,     0.2124105E+03_r8,     0.2111156E+03_r8,     0.2098226E+03_r8,     0.2085321E+03_r8,&  
         0.2072448E+03_r8,     0.2059612E+03_r8,     0.2046818E+03_r8,     0.2034070E+03_r8,     0.2021372E+03_r8,&  
         0.2008728E+03_r8,     0.1996141E+03_r8,     0.1983612E+03_r8,     0.1971146E+03_r8,     0.1958743E+03_r8,&  
         0.1946405E+03_r8,     0.1934133E+03_r8,     0.1921930E+03_r8,     0.1909795E+03_r8,     0.1897730E+03_r8,&  
         0.1885734E+03_r8,     0.1873810E+03_r8,     0.1861957E+03_r8,     0.1850174E+03_r8,     0.1838464E+03_r8,&  
         0.1826824E+03_r8,     0.1815256E+03_r8,     0.1803759E+03_r8,     0.1792334E+03_r8,     0.1780979E+03_r8/)  

    psaditmk(1:150,   99)= (/ &
         0.3134879E+03_r8,     0.3128444E+03_r8,     0.3122015E+03_r8,     0.3115587E+03_r8,     0.3109171E+03_r8,&  
         0.3102758E+03_r8,     0.3096345E+03_r8,     0.3089943E+03_r8,     0.3083537E+03_r8,     0.3077141E+03_r8,&  
         0.3070745E+03_r8,     0.3064355E+03_r8,     0.3057964E+03_r8,     0.3051576E+03_r8,     0.3045187E+03_r8,&  
         0.3038800E+03_r8,     0.3032418E+03_r8,     0.3026031E+03_r8,     0.3019646E+03_r8,     0.3013264E+03_r8,&  
         0.3006873E+03_r8,     0.3000485E+03_r8,     0.2994093E+03_r8,     0.2987702E+03_r8,     0.2981299E+03_r8,&  
         0.2974900E+03_r8,     0.2968494E+03_r8,     0.2962082E+03_r8,     0.2955669E+03_r8,     0.2949243E+03_r8,&  
         0.2942814E+03_r8,     0.2936377E+03_r8,     0.2929928E+03_r8,     0.2923480E+03_r8,     0.2917013E+03_r8,&  
         0.2910540E+03_r8,     0.2904054E+03_r8,     0.2897557E+03_r8,     0.2891048E+03_r8,     0.2884525E+03_r8,&  
         0.2877991E+03_r8,     0.2871436E+03_r8,     0.2864869E+03_r8,     0.2858284E+03_r8,     0.2851682E+03_r8,&  
         0.2845058E+03_r8,     0.2838417E+03_r8,     0.2831756E+03_r8,     0.2825071E+03_r8,     0.2818364E+03_r8,&  
         0.2811631E+03_r8,     0.2804872E+03_r8,     0.2798087E+03_r8,     0.2791274E+03_r8,     0.2784431E+03_r8,&  
         0.2777556E+03_r8,     0.2770651E+03_r8,     0.2763710E+03_r8,     0.2756734E+03_r8,     0.2749721E+03_r8,&  
         0.2742668E+03_r8,     0.2735575E+03_r8,     0.2728439E+03_r8,     0.2721259E+03_r8,     0.2714032E+03_r8,&  
         0.2706757E+03_r8,     0.2699431E+03_r8,     0.2692049E+03_r8,     0.2684613E+03_r8,     0.2677117E+03_r8,&  
         0.2669563E+03_r8,     0.2661946E+03_r8,     0.2654263E+03_r8,     0.2646510E+03_r8,     0.2638688E+03_r8,&  
         0.2630788E+03_r8,     0.2622809E+03_r8,     0.2614755E+03_r8,     0.2606611E+03_r8,     0.2598381E+03_r8,&  
         0.2590059E+03_r8,     0.2581643E+03_r8,     0.2573126E+03_r8,     0.2564506E+03_r8,     0.2555780E+03_r8,&  
         0.2546943E+03_r8,     0.2537990E+03_r8,     0.2528921E+03_r8,     0.2519729E+03_r8,     0.2510407E+03_r8,&  
         0.2500956E+03_r8,     0.2491371E+03_r8,     0.2481647E+03_r8,     0.2471781E+03_r8,     0.2461770E+03_r8,&  
         0.2451614E+03_r8,     0.2441305E+03_r8,     0.2430844E+03_r8,     0.2420230E+03_r8,     0.2409462E+03_r8,&  
         0.2398536E+03_r8,     0.2387456E+03_r8,     0.2376223E+03_r8,     0.2364838E+03_r8,     0.2353303E+03_r8,&  
         0.2341624E+03_r8,     0.2329803E+03_r8,     0.2317847E+03_r8,     0.2305762E+03_r8,     0.2293555E+03_r8,&  
         0.2281235E+03_r8,     0.2268808E+03_r8,     0.2256286E+03_r8,     0.2243676E+03_r8,     0.2230990E+03_r8,&  
         0.2218237E+03_r8,     0.2205427E+03_r8,     0.2192570E+03_r8,     0.2179678E+03_r8,     0.2166758E+03_r8,&  
         0.2153820E+03_r8,     0.2140874E+03_r8,     0.2127927E+03_r8,     0.2114989E+03_r8,     0.2102066E+03_r8,&  
         0.2089164E+03_r8,     0.2076292E+03_r8,     0.2063453E+03_r8,     0.2050654E+03_r8,     0.2037898E+03_r8,&  
         0.2025191E+03_r8,     0.2012536E+03_r8,     0.1999935E+03_r8,     0.1987393E+03_r8,     0.1974910E+03_r8,&  
         0.1962490E+03_r8,     0.1950135E+03_r8,     0.1937845E+03_r8,     0.1925623E+03_r8,     0.1913468E+03_r8,&  
         0.1901383E+03_r8,     0.1889368E+03_r8,     0.1877423E+03_r8,     0.1865548E+03_r8,     0.1853745E+03_r8,&  
         0.1842013E+03_r8,     0.1830352E+03_r8,     0.1818763E+03_r8,     0.1807245E+03_r8,     0.1795798E+03_r8/)  

    psaditmk(1:150,  100)= (/ &
         0.3138205E+03_r8,     0.3131791E+03_r8,     0.3125380E+03_r8,     0.3118984E+03_r8,     0.3112584E+03_r8,&  
         0.3106193E+03_r8,     0.3099809E+03_r8,     0.3093425E+03_r8,     0.3087050E+03_r8,     0.3080672E+03_r8,&  
         0.3074305E+03_r8,     0.3067943E+03_r8,     0.3061576E+03_r8,     0.3055215E+03_r8,     0.3048855E+03_r8,&  
         0.3042497E+03_r8,     0.3036142E+03_r8,     0.3029786E+03_r8,     0.3023435E+03_r8,     0.3017075E+03_r8,&  
         0.3010721E+03_r8,     0.3004362E+03_r8,     0.2998004E+03_r8,     0.2991643E+03_r8,     0.2985278E+03_r8,&  
         0.2978913E+03_r8,     0.2972543E+03_r8,     0.2966171E+03_r8,     0.2959790E+03_r8,     0.2953407E+03_r8,&  
         0.2947018E+03_r8,     0.2940620E+03_r8,     0.2934217E+03_r8,     0.2927804E+03_r8,     0.2921385E+03_r8,&  
         0.2914957E+03_r8,     0.2908515E+03_r8,     0.2902067E+03_r8,     0.2895603E+03_r8,     0.2889127E+03_r8,&  
         0.2882644E+03_r8,     0.2876142E+03_r8,     0.2869630E+03_r8,     0.2863098E+03_r8,     0.2856549E+03_r8,&  
         0.2849984E+03_r8,     0.2843405E+03_r8,     0.2836804E+03_r8,     0.2830181E+03_r8,     0.2823534E+03_r8,&  
         0.2816867E+03_r8,     0.2810179E+03_r8,     0.2803465E+03_r8,     0.2796724E+03_r8,     0.2789953E+03_r8,&  
         0.2783156E+03_r8,     0.2776328E+03_r8,     0.2769467E+03_r8,     0.2762573E+03_r8,     0.2755647E+03_r8,&  
         0.2748683E+03_r8,     0.2741682E+03_r8,     0.2734639E+03_r8,     0.2727558E+03_r8,     0.2720428E+03_r8,&  
         0.2713256E+03_r8,     0.2706035E+03_r8,     0.2698766E+03_r8,     0.2691444E+03_r8,     0.2684067E+03_r8,&  
         0.2676637E+03_r8,     0.2669141E+03_r8,     0.2661589E+03_r8,     0.2653971E+03_r8,     0.2646284E+03_r8,&  
         0.2638532E+03_r8,     0.2630704E+03_r8,     0.2622798E+03_r8,     0.2614815E+03_r8,     0.2606749E+03_r8,&  
         0.2598593E+03_r8,     0.2590353E+03_r8,     0.2582016E+03_r8,     0.2573585E+03_r8,     0.2565049E+03_r8,&  
         0.2556413E+03_r8,     0.2547664E+03_r8,     0.2538806E+03_r8,     0.2529831E+03_r8,     0.2520732E+03_r8,&  
         0.2511510E+03_r8,     0.2502160E+03_r8,     0.2492675E+03_r8,     0.2483056E+03_r8,     0.2473297E+03_r8,&  
         0.2463392E+03_r8,     0.2453342E+03_r8,     0.2443144E+03_r8,     0.2432793E+03_r8,     0.2422289E+03_r8,&  
         0.2411630E+03_r8,     0.2400815E+03_r8,     0.2389844E+03_r8,     0.2378719E+03_r8,     0.2367440E+03_r8,&  
         0.2356009E+03_r8,     0.2344432E+03_r8,     0.2332709E+03_r8,     0.2320849E+03_r8,     0.2308854E+03_r8,&  
         0.2296734E+03_r8,     0.2284494E+03_r8,     0.2272144E+03_r8,     0.2259693E+03_r8,     0.2247149E+03_r8,&  
         0.2234522E+03_r8,     0.2221822E+03_r8,     0.2209061E+03_r8,     0.2196247E+03_r8,     0.2183391E+03_r8,&  
         0.2170502E+03_r8,     0.2157591E+03_r8,     0.2144666E+03_r8,     0.2131737E+03_r8,     0.2118811E+03_r8,&  
         0.2105897E+03_r8,     0.2093000E+03_r8,     0.2080130E+03_r8,     0.2067290E+03_r8,     0.2054486E+03_r8,&  
         0.2041725E+03_r8,     0.2029009E+03_r8,     0.2016343E+03_r8,     0.2003730E+03_r8,     0.1991174E+03_r8,&  
         0.1978677E+03_r8,     0.1966241E+03_r8,     0.1953868E+03_r8,     0.1941561E+03_r8,     0.1929320E+03_r8,&  
         0.1917146E+03_r8,     0.1905041E+03_r8,     0.1893006E+03_r8,     0.1881040E+03_r8,     0.1869145E+03_r8,&  
         0.1857321E+03_r8,     0.1845568E+03_r8,     0.1833886E+03_r8,     0.1822276E+03_r8,     0.1810736E+03_r8/)  

    psaditmk(1:150,  101)= (/ &
         0.3141481E+03_r8,     0.3135089E+03_r8,     0.3128702E+03_r8,     0.3122318E+03_r8,     0.3115943E+03_r8,&  
         0.3109576E+03_r8,     0.3103210E+03_r8,     0.3096855E+03_r8,     0.3090501E+03_r8,     0.3084149E+03_r8,&  
         0.3077809E+03_r8,     0.3071463E+03_r8,     0.3065129E+03_r8,     0.3058790E+03_r8,     0.3052459E+03_r8,&  
         0.3046131E+03_r8,     0.3039799E+03_r8,     0.3033476E+03_r8,     0.3027147E+03_r8,     0.3020822E+03_r8,&  
         0.3014498E+03_r8,     0.3008171E+03_r8,     0.3001846E+03_r8,     0.2995515E+03_r8,     0.2989186E+03_r8,&  
         0.2982850E+03_r8,     0.2976519E+03_r8,     0.2970179E+03_r8,     0.2963835E+03_r8,     0.2957490E+03_r8,&  
         0.2951137E+03_r8,     0.2944779E+03_r8,     0.2938414E+03_r8,     0.2932043E+03_r8,     0.2925666E+03_r8,&  
         0.2919277E+03_r8,     0.2912882E+03_r8,     0.2906477E+03_r8,     0.2900061E+03_r8,     0.2893637E+03_r8,&  
         0.2887196E+03_r8,     0.2880746E+03_r8,     0.2874283E+03_r8,     0.2867802E+03_r8,     0.2861312E+03_r8,&  
         0.2854802E+03_r8,     0.2848279E+03_r8,     0.2841734E+03_r8,     0.2835168E+03_r8,     0.2828586E+03_r8,&  
         0.2821984E+03_r8,     0.2815357E+03_r8,     0.2808711E+03_r8,     0.2802035E+03_r8,     0.2795338E+03_r8,&  
         0.2788613E+03_r8,     0.2781858E+03_r8,     0.2775075E+03_r8,     0.2768262E+03_r8,     0.2761416E+03_r8,&  
         0.2754535E+03_r8,     0.2747619E+03_r8,     0.2740666E+03_r8,     0.2733676E+03_r8,     0.2726645E+03_r8,&  
         0.2719570E+03_r8,     0.2712451E+03_r8,     0.2705285E+03_r8,     0.2698069E+03_r8,     0.2690807E+03_r8,&  
         0.2683490E+03_r8,     0.2676117E+03_r8,     0.2668689E+03_r8,     0.2661196E+03_r8,     0.2653641E+03_r8,&  
         0.2646022E+03_r8,     0.2638336E+03_r8,     0.2630577E+03_r8,     0.2622743E+03_r8,     0.2614831E+03_r8,&  
         0.2606840E+03_r8,     0.2598763E+03_r8,     0.2590599E+03_r8,     0.2582345E+03_r8,     0.2573994E+03_r8,&  
         0.2565545E+03_r8,     0.2556993E+03_r8,     0.2548335E+03_r8,     0.2539567E+03_r8,     0.2530683E+03_r8,&  
         0.2521680E+03_r8,     0.2512555E+03_r8,     0.2503304E+03_r8,     0.2493923E+03_r8,     0.2484403E+03_r8,&  
         0.2474749E+03_r8,     0.2464952E+03_r8,     0.2455009E+03_r8,     0.2444917E+03_r8,     0.2434677E+03_r8,&  
         0.2424282E+03_r8,     0.2413733E+03_r8,     0.2403028E+03_r8,     0.2392169E+03_r8,     0.2381152E+03_r8,&  
         0.2369979E+03_r8,     0.2358655E+03_r8,     0.2347180E+03_r8,     0.2335558E+03_r8,     0.2323794E+03_r8,&  
         0.2311892E+03_r8,     0.2299861E+03_r8,     0.2287706E+03_r8,     0.2275435E+03_r8,     0.2263057E+03_r8,&  
         0.2250582E+03_r8,     0.2238018E+03_r8,     0.2225375E+03_r8,     0.2212665E+03_r8,     0.2199896E+03_r8,&  
         0.2187080E+03_r8,     0.2174226E+03_r8,     0.2161344E+03_r8,     0.2148443E+03_r8,     0.2135533E+03_r8,&  
         0.2122622E+03_r8,     0.2109718E+03_r8,     0.2096830E+03_r8,     0.2083962E+03_r8,     0.2071122E+03_r8,&  
         0.2058317E+03_r8,     0.2045550E+03_r8,     0.2032826E+03_r8,     0.2020151E+03_r8,     0.2007527E+03_r8,&  
         0.1994958E+03_r8,     0.1982446E+03_r8,     0.1969995E+03_r8,     0.1957606E+03_r8,     0.1945281E+03_r8,&  
         0.1933021E+03_r8,     0.1920829E+03_r8,     0.1908705E+03_r8,     0.1896649E+03_r8,     0.1884664E+03_r8,&  
         0.1872748E+03_r8,     0.1860903E+03_r8,     0.1849129E+03_r8,     0.1837426E+03_r8,     0.1825794E+03_r8/)  

    psaditmk(1:150,  102)= (/ &
         0.3144706E+03_r8,     0.3138331E+03_r8,     0.3131963E+03_r8,     0.3125606E+03_r8,     0.3119251E+03_r8,&  
         0.3112902E+03_r8,     0.3106564E+03_r8,     0.3100224E+03_r8,     0.3093898E+03_r8,     0.3087573E+03_r8,&  
         0.3081248E+03_r8,     0.3074932E+03_r8,     0.3068618E+03_r8,     0.3062308E+03_r8,     0.3056004E+03_r8,&  
         0.3049701E+03_r8,     0.3043398E+03_r8,     0.3037099E+03_r8,     0.3030800E+03_r8,     0.3024505E+03_r8,&  
         0.3018207E+03_r8,     0.3011914E+03_r8,     0.3005616E+03_r8,     0.2999321E+03_r8,     0.2993021E+03_r8,&  
         0.2986717E+03_r8,     0.2980417E+03_r8,     0.2974110E+03_r8,     0.2967805E+03_r8,     0.2961492E+03_r8,&  
         0.2955176E+03_r8,     0.2948860E+03_r8,     0.2942530E+03_r8,     0.2936200E+03_r8,     0.2929861E+03_r8,&  
         0.2923513E+03_r8,     0.2917163E+03_r8,     0.2910798E+03_r8,     0.2904427E+03_r8,     0.2898046E+03_r8,&  
         0.2891653E+03_r8,     0.2885251E+03_r8,     0.2878833E+03_r8,     0.2872405E+03_r8,     0.2865965E+03_r8,&  
         0.2859507E+03_r8,     0.2853036E+03_r8,     0.2846546E+03_r8,     0.2840040E+03_r8,     0.2833520E+03_r8,&  
         0.2826975E+03_r8,     0.2820411E+03_r8,     0.2813824E+03_r8,     0.2807218E+03_r8,     0.2800587E+03_r8,&  
         0.2793931E+03_r8,     0.2787245E+03_r8,     0.2780539E+03_r8,     0.2773799E+03_r8,     0.2767032E+03_r8,&  
         0.2760232E+03_r8,     0.2753399E+03_r8,     0.2746530E+03_r8,     0.2739628E+03_r8,     0.2732686E+03_r8,&  
         0.2725705E+03_r8,     0.2718682E+03_r8,     0.2711616E+03_r8,     0.2704505E+03_r8,     0.2697345E+03_r8,&  
         0.2690138E+03_r8,     0.2682878E+03_r8,     0.2675565E+03_r8,     0.2668196E+03_r8,     0.2660768E+03_r8,&  
         0.2653275E+03_r8,     0.2645721E+03_r8,     0.2638102E+03_r8,     0.2630411E+03_r8,     0.2622646E+03_r8,&  
         0.2614806E+03_r8,     0.2606888E+03_r8,     0.2598887E+03_r8,     0.2590800E+03_r8,     0.2582625E+03_r8,&  
         0.2574352E+03_r8,     0.2565990E+03_r8,     0.2557523E+03_r8,     0.2548953E+03_r8,     0.2540272E+03_r8,&  
         0.2531480E+03_r8,     0.2522572E+03_r8,     0.2513542E+03_r8,     0.2504388E+03_r8,     0.2495106E+03_r8,&  
         0.2485692E+03_r8,     0.2476138E+03_r8,     0.2466446E+03_r8,     0.2456611E+03_r8,     0.2446627E+03_r8,&  
         0.2436495E+03_r8,     0.2426210E+03_r8,     0.2415772E+03_r8,     0.2405177E+03_r8,     0.2394427E+03_r8,&  
         0.2383519E+03_r8,     0.2372456E+03_r8,     0.2361237E+03_r8,     0.2349866E+03_r8,     0.2338347E+03_r8,&  
         0.2326682E+03_r8,     0.2314875E+03_r8,     0.2302935E+03_r8,     0.2290867E+03_r8,     0.2278678E+03_r8,&  
         0.2266378E+03_r8,     0.2253974E+03_r8,     0.2241476E+03_r8,     0.2228893E+03_r8,     0.2216237E+03_r8,&  
         0.2203517E+03_r8,     0.2190744E+03_r8,     0.2177927E+03_r8,     0.2165077E+03_r8,     0.2152203E+03_r8,&  
         0.2139315E+03_r8,     0.2126421E+03_r8,     0.2113531E+03_r8,     0.2100650E+03_r8,     0.2087788E+03_r8,&  
         0.2074951E+03_r8,     0.2062144E+03_r8,     0.2049373E+03_r8,     0.2036643E+03_r8,     0.2023959E+03_r8,&  
         0.2011325E+03_r8,     0.1998744E+03_r8,     0.1986218E+03_r8,     0.1973752E+03_r8,     0.1961347E+03_r8,&  
         0.1949005E+03_r8,     0.1936728E+03_r8,     0.1924517E+03_r8,     0.1912374E+03_r8,     0.1900299E+03_r8,&  
         0.1888293E+03_r8,     0.1876357E+03_r8,     0.1864491E+03_r8,     0.1852696E+03_r8,     0.1840972E+03_r8/)  

    psaditmk(1:150,  103)= (/ &
         0.3147881E+03_r8,     0.3141526E+03_r8,     0.3135178E+03_r8,     0.3128837E+03_r8,     0.3122507E+03_r8,&  
         0.3116180E+03_r8,     0.3109859E+03_r8,     0.3103548E+03_r8,     0.3097237E+03_r8,     0.3090932E+03_r8,&  
         0.3084637E+03_r8,     0.3078343E+03_r8,     0.3072055E+03_r8,     0.3065770E+03_r8,     0.3059487E+03_r8,&  
         0.3053212E+03_r8,     0.3046935E+03_r8,     0.3040662E+03_r8,     0.3034393E+03_r8,     0.3028120E+03_r8,&  
         0.3021855E+03_r8,     0.3015584E+03_r8,     0.3009320E+03_r8,     0.3003055E+03_r8,     0.2996783E+03_r8,&  
         0.2990516E+03_r8,     0.2984244E+03_r8,     0.2977974E+03_r8,     0.2971699E+03_r8,     0.2965420E+03_r8,&  
         0.2959142E+03_r8,     0.2952855E+03_r8,     0.2946565E+03_r8,     0.2940271E+03_r8,     0.2933971E+03_r8,&  
         0.2927666E+03_r8,     0.2921353E+03_r8,     0.2915032E+03_r8,     0.2908705E+03_r8,     0.2902363E+03_r8,&  
         0.2896019E+03_r8,     0.2889659E+03_r8,     0.2883290E+03_r8,     0.2876911E+03_r8,     0.2870517E+03_r8,&  
         0.2864109E+03_r8,     0.2857691E+03_r8,     0.2851255E+03_r8,     0.2844804E+03_r8,     0.2838337E+03_r8,&  
         0.2831848E+03_r8,     0.2825345E+03_r8,     0.2818820E+03_r8,     0.2812275E+03_r8,     0.2805711E+03_r8,&  
         0.2799117E+03_r8,     0.2792504E+03_r8,     0.2785864E+03_r8,     0.2779198E+03_r8,     0.2772503E+03_r8,&  
         0.2765779E+03_r8,     0.2759024E+03_r8,     0.2752236E+03_r8,     0.2745418E+03_r8,     0.2738563E+03_r8,&  
         0.2731671E+03_r8,     0.2724738E+03_r8,     0.2717768E+03_r8,     0.2710755E+03_r8,     0.2703697E+03_r8,&  
         0.2696591E+03_r8,     0.2689436E+03_r8,     0.2682235E+03_r8,     0.2674980E+03_r8,     0.2667671E+03_r8,&  
         0.2660302E+03_r8,     0.2652872E+03_r8,     0.2645382E+03_r8,     0.2637829E+03_r8,     0.2630205E+03_r8,&  
         0.2622509E+03_r8,     0.2614739E+03_r8,     0.2606894E+03_r8,     0.2598966E+03_r8,     0.2590955E+03_r8,&  
         0.2582857E+03_r8,     0.2574669E+03_r8,     0.2566382E+03_r8,     0.2558000E+03_r8,     0.2549516E+03_r8,&  
         0.2540925E+03_r8,     0.2532223E+03_r8,     0.2523406E+03_r8,     0.2514471E+03_r8,     0.2505415E+03_r8,&  
         0.2496231E+03_r8,     0.2486916E+03_r8,     0.2477468E+03_r8,     0.2467878E+03_r8,     0.2458147E+03_r8,&  
         0.2448273E+03_r8,     0.2438248E+03_r8,     0.2428074E+03_r8,     0.2417745E+03_r8,     0.2407261E+03_r8,&  
         0.2396620E+03_r8,     0.2385822E+03_r8,     0.2374868E+03_r8,     0.2363757E+03_r8,     0.2352493E+03_r8,&  
         0.2341076E+03_r8,     0.2329511E+03_r8,     0.2317803E+03_r8,     0.2305956E+03_r8,     0.2293977E+03_r8,&  
         0.2281874E+03_r8,     0.2269653E+03_r8,     0.2257323E+03_r8,     0.2244894E+03_r8,     0.2232375E+03_r8,&  
         0.2219776E+03_r8,     0.2207108E+03_r8,     0.2194380E+03_r8,     0.2181604E+03_r8,     0.2168789E+03_r8,&  
         0.2155944E+03_r8,     0.2143081E+03_r8,     0.2130207E+03_r8,     0.2117331E+03_r8,     0.2104462E+03_r8,&  
         0.2091607E+03_r8,     0.2078773E+03_r8,     0.2065966E+03_r8,     0.2053194E+03_r8,     0.2040459E+03_r8,&  
         0.2027767E+03_r8,     0.2015124E+03_r8,     0.2002531E+03_r8,     0.1989993E+03_r8,     0.1977512E+03_r8,&  
         0.1965092E+03_r8,     0.1952733E+03_r8,     0.1940439E+03_r8,     0.1928209E+03_r8,     0.1916047E+03_r8,&  
         0.1903953E+03_r8,     0.1891927E+03_r8,     0.1879971E+03_r8,     0.1868085E+03_r8,     0.1856269E+03_r8/)  

    psaditmk(1:150,  104)= (/ &
         0.3151007E+03_r8,     0.3144673E+03_r8,     0.3138347E+03_r8,     0.3132028E+03_r8,     0.3125714E+03_r8,&  
         0.3119409E+03_r8,     0.3113111E+03_r8,     0.3106812E+03_r8,     0.3100530E+03_r8,     0.3094247E+03_r8,&  
         0.3087971E+03_r8,     0.3081703E+03_r8,     0.3075434E+03_r8,     0.3069172E+03_r8,     0.3062919E+03_r8,&  
         0.3056663E+03_r8,     0.3050415E+03_r8,     0.3044165E+03_r8,     0.3037917E+03_r8,     0.3031677E+03_r8,&  
         0.3025434E+03_r8,     0.3019198E+03_r8,     0.3012957E+03_r8,     0.3006718E+03_r8,     0.3000481E+03_r8,&  
         0.2994242E+03_r8,     0.2988004E+03_r8,     0.2981764E+03_r8,     0.2975519E+03_r8,     0.2969279E+03_r8,&  
         0.2963029E+03_r8,     0.2956778E+03_r8,     0.2950528E+03_r8,     0.2944266E+03_r8,     0.2938003E+03_r8,&  
         0.2931736E+03_r8,     0.2925460E+03_r8,     0.2919182E+03_r8,     0.2912891E+03_r8,     0.2906595E+03_r8,&  
         0.2900292E+03_r8,     0.2893975E+03_r8,     0.2887652E+03_r8,     0.2881318E+03_r8,     0.2874969E+03_r8,&  
         0.2868615E+03_r8,     0.2862242E+03_r8,     0.2855857E+03_r8,     0.2849460E+03_r8,     0.2843044E+03_r8,&  
         0.2836612E+03_r8,     0.2830165E+03_r8,     0.2823699E+03_r8,     0.2817213E+03_r8,     0.2810705E+03_r8,&  
         0.2804181E+03_r8,     0.2797630E+03_r8,     0.2791058E+03_r8,     0.2784460E+03_r8,     0.2777836E+03_r8,&  
         0.2771183E+03_r8,     0.2764504E+03_r8,     0.2757795E+03_r8,     0.2751054E+03_r8,     0.2744279E+03_r8,&  
         0.2737471E+03_r8,     0.2730625E+03_r8,     0.2723746E+03_r8,     0.2716824E+03_r8,     0.2709861E+03_r8,&  
         0.2702856E+03_r8,     0.2695806E+03_r8,     0.2688709E+03_r8,     0.2681561E+03_r8,     0.2674362E+03_r8,&  
         0.2667112E+03_r8,     0.2659804E+03_r8,     0.2652436E+03_r8,     0.2645008E+03_r8,     0.2637519E+03_r8,&  
         0.2629960E+03_r8,     0.2622332E+03_r8,     0.2614631E+03_r8,     0.2606855E+03_r8,     0.2599002E+03_r8,&  
         0.2591066E+03_r8,     0.2583044E+03_r8,     0.2574934E+03_r8,     0.2566730E+03_r8,     0.2558430E+03_r8,&  
         0.2550029E+03_r8,     0.2541524E+03_r8,     0.2532912E+03_r8,     0.2524186E+03_r8,     0.2515345E+03_r8,&  
         0.2506384E+03_r8,     0.2497297E+03_r8,     0.2488083E+03_r8,     0.2478732E+03_r8,     0.2469248E+03_r8,&  
         0.2459623E+03_r8,     0.2449852E+03_r8,     0.2439938E+03_r8,     0.2429871E+03_r8,     0.2419653E+03_r8,&  
         0.2409280E+03_r8,     0.2398749E+03_r8,     0.2388061E+03_r8,     0.2377216E+03_r8,     0.2366215E+03_r8,&  
         0.2355056E+03_r8,     0.2343744E+03_r8,     0.2332283E+03_r8,     0.2320673E+03_r8,     0.2308922E+03_r8,&  
         0.2297035E+03_r8,     0.2285019E+03_r8,     0.2272880E+03_r8,     0.2260628E+03_r8,     0.2248271E+03_r8,&  
         0.2235818E+03_r8,     0.2223280E+03_r8,     0.2210668E+03_r8,     0.2197989E+03_r8,     0.2185256E+03_r8,&  
         0.2172478E+03_r8,     0.2159666E+03_r8,     0.2146830E+03_r8,     0.2133978E+03_r8,     0.2121120E+03_r8,&  
         0.2108264E+03_r8,     0.2095418E+03_r8,     0.2082589E+03_r8,     0.2069785E+03_r8,     0.2057011E+03_r8,&  
         0.2044272E+03_r8,     0.2031575E+03_r8,     0.2018922E+03_r8,     0.2006320E+03_r8,     0.1993769E+03_r8,&  
         0.1981275E+03_r8,     0.1968840E+03_r8,     0.1956465E+03_r8,     0.1944154E+03_r8,     0.1931907E+03_r8,&  
         0.1919726E+03_r8,     0.1907613E+03_r8,     0.1895567E+03_r8,     0.1883591E+03_r8,     0.1871685E+03_r8/)  

    psaditmk(1:150,  105)= (/ &
         0.3154095E+03_r8,     0.3147774E+03_r8,     0.3141463E+03_r8,     0.3135166E+03_r8,     0.3128871E+03_r8,&  
         0.3122584E+03_r8,     0.3116304E+03_r8,     0.3110033E+03_r8,     0.3103764E+03_r8,     0.3097506E+03_r8,&  
         0.3091255E+03_r8,     0.3085001E+03_r8,     0.3078756E+03_r8,     0.3072524E+03_r8,     0.3066284E+03_r8,&  
         0.3060059E+03_r8,     0.3053832E+03_r8,     0.3047609E+03_r8,     0.3041389E+03_r8,     0.3035172E+03_r8,&  
         0.3028957E+03_r8,     0.3022744E+03_r8,     0.3016532E+03_r8,     0.3010323E+03_r8,     0.3004113E+03_r8,&  
         0.2997904E+03_r8,     0.2991697E+03_r8,     0.2985485E+03_r8,     0.2979274E+03_r8,     0.2973063E+03_r8,&  
         0.2966844E+03_r8,     0.2960628E+03_r8,     0.2954408E+03_r8,     0.2948187E+03_r8,     0.2941957E+03_r8,&  
         0.2935725E+03_r8,     0.2929489E+03_r8,     0.2923246E+03_r8,     0.2916996E+03_r8,     0.2910741E+03_r8,&  
         0.2904475E+03_r8,     0.2898203E+03_r8,     0.2891926E+03_r8,     0.2885631E+03_r8,     0.2879331E+03_r8,&  
         0.2873020E+03_r8,     0.2866693E+03_r8,     0.2860358E+03_r8,     0.2854011E+03_r8,     0.2847645E+03_r8,&  
         0.2841266E+03_r8,     0.2834874E+03_r8,     0.2828465E+03_r8,     0.2822034E+03_r8,     0.2815588E+03_r8,&  
         0.2809121E+03_r8,     0.2802632E+03_r8,     0.2796124E+03_r8,     0.2789592E+03_r8,     0.2783034E+03_r8,&  
         0.2776452E+03_r8,     0.2769845E+03_r8,     0.2763208E+03_r8,     0.2756541E+03_r8,     0.2749846E+03_r8,&  
         0.2743118E+03_r8,     0.2736357E+03_r8,     0.2729561E+03_r8,     0.2722728E+03_r8,     0.2715856E+03_r8,&  
         0.2708944E+03_r8,     0.2701988E+03_r8,     0.2694991E+03_r8,     0.2687947E+03_r8,     0.2680858E+03_r8,&  
         0.2673714E+03_r8,     0.2666519E+03_r8,     0.2659272E+03_r8,     0.2651965E+03_r8,     0.2644598E+03_r8,&  
         0.2637170E+03_r8,     0.2629676E+03_r8,     0.2622116E+03_r8,     0.2614482E+03_r8,     0.2606779E+03_r8,&  
         0.2598996E+03_r8,     0.2591133E+03_r8,     0.2583186E+03_r8,     0.2575153E+03_r8,     0.2567029E+03_r8,&  
         0.2558810E+03_r8,     0.2550495E+03_r8,     0.2542074E+03_r8,     0.2533550E+03_r8,     0.2524914E+03_r8,&  
         0.2516164E+03_r8,     0.2507296E+03_r8,     0.2498306E+03_r8,     0.2489188E+03_r8,     0.2479942E+03_r8,&  
         0.2470556E+03_r8,     0.2461037E+03_r8,     0.2451372E+03_r8,     0.2441562E+03_r8,     0.2431605E+03_r8,&  
         0.2421495E+03_r8,     0.2411232E+03_r8,     0.2400812E+03_r8,     0.2390234E+03_r8,     0.2379499E+03_r8,&  
         0.2368607E+03_r8,     0.2357557E+03_r8,     0.2346352E+03_r8,     0.2334993E+03_r8,     0.2323485E+03_r8,&  
         0.2311832E+03_r8,     0.2300039E+03_r8,     0.2288113E+03_r8,     0.2276060E+03_r8,     0.2263887E+03_r8,&  
         0.2251606E+03_r8,     0.2239223E+03_r8,     0.2226750E+03_r8,     0.2214194E+03_r8,     0.2201567E+03_r8,&  
         0.2188881E+03_r8,     0.2176144E+03_r8,     0.2163367E+03_r8,     0.2150560E+03_r8,     0.2137733E+03_r8,&  
         0.2124895E+03_r8,     0.2112054E+03_r8,     0.2099220E+03_r8,     0.2086398E+03_r8,     0.2073598E+03_r8,&  
         0.2060824E+03_r8,     0.2048083E+03_r8,     0.2035381E+03_r8,     0.2022721E+03_r8,     0.2010109E+03_r8,&  
         0.1997548E+03_r8,     0.1985041E+03_r8,     0.1972591E+03_r8,     0.1960201E+03_r8,     0.1947873E+03_r8,&  
         0.1935609E+03_r8,     0.1923410E+03_r8,     0.1911278E+03_r8,     0.1899213E+03_r8,     0.1887217E+03_r8/)  

    psaditmk(1:150,  106)= (/ &
         0.3157130E+03_r8,     0.3150834E+03_r8,     0.3144542E+03_r8,     0.3138255E+03_r8,     0.3131980E+03_r8,&  
         0.3125716E+03_r8,     0.3119454E+03_r8,     0.3113199E+03_r8,     0.3106958E+03_r8,     0.3100717E+03_r8,&  
         0.3094480E+03_r8,     0.3088256E+03_r8,     0.3082034E+03_r8,     0.3075814E+03_r8,     0.3069606E+03_r8,&  
         0.3063394E+03_r8,     0.3057195E+03_r8,     0.3050998E+03_r8,     0.3044800E+03_r8,     0.3038611E+03_r8,&  
         0.3032421E+03_r8,     0.3026230E+03_r8,     0.3020051E+03_r8,     0.3013864E+03_r8,     0.3007679E+03_r8,&  
         0.3001503E+03_r8,     0.2995319E+03_r8,     0.2989139E+03_r8,     0.2982960E+03_r8,     0.2976775E+03_r8,&  
         0.2970594E+03_r8,     0.2964408E+03_r8,     0.2958221E+03_r8,     0.2952031E+03_r8,     0.2945836E+03_r8,&  
         0.2939640E+03_r8,     0.2933439E+03_r8,     0.2927231E+03_r8,     0.2921021E+03_r8,     0.2914801E+03_r8,&  
         0.2908577E+03_r8,     0.2902348E+03_r8,     0.2896105E+03_r8,     0.2889857E+03_r8,     0.2883603E+03_r8,&  
         0.2877330E+03_r8,     0.2871052E+03_r8,     0.2864765E+03_r8,     0.2858466E+03_r8,     0.2852146E+03_r8,&  
         0.2845821E+03_r8,     0.2839478E+03_r8,     0.2833122E+03_r8,     0.2826746E+03_r8,     0.2820356E+03_r8,&  
         0.2813945E+03_r8,     0.2807517E+03_r8,     0.2801068E+03_r8,     0.2794598E+03_r8,     0.2788106E+03_r8,&  
         0.2781588E+03_r8,     0.2775050E+03_r8,     0.2768482E+03_r8,     0.2761890E+03_r8,     0.2755269E+03_r8,&  
         0.2748616E+03_r8,     0.2741934E+03_r8,     0.2735219E+03_r8,     0.2728470E+03_r8,     0.2721683E+03_r8,&  
         0.2714861E+03_r8,     0.2708000E+03_r8,     0.2701095E+03_r8,     0.2694148E+03_r8,     0.2687157E+03_r8,&  
         0.2680120E+03_r8,     0.2673034E+03_r8,     0.2665895E+03_r8,     0.2658704E+03_r8,     0.2651461E+03_r8,&  
         0.2644153E+03_r8,     0.2636785E+03_r8,     0.2629357E+03_r8,     0.2621861E+03_r8,     0.2614295E+03_r8,&  
         0.2606660E+03_r8,     0.2598948E+03_r8,     0.2591157E+03_r8,     0.2583283E+03_r8,     0.2575327E+03_r8,&  
         0.2567282E+03_r8,     0.2559142E+03_r8,     0.2550906E+03_r8,     0.2542573E+03_r8,     0.2534134E+03_r8,&  
         0.2525587E+03_r8,     0.2516927E+03_r8,     0.2508154E+03_r8,     0.2499256E+03_r8,     0.2490238E+03_r8,&  
         0.2481087E+03_r8,     0.2471807E+03_r8,     0.2462386E+03_r8,     0.2452830E+03_r8,     0.2443125E+03_r8,&  
         0.2433275E+03_r8,     0.2423275E+03_r8,     0.2413119E+03_r8,     0.2402809E+03_r8,     0.2392343E+03_r8,&  
         0.2381719E+03_r8,     0.2370936E+03_r8,     0.2359994E+03_r8,     0.2348896E+03_r8,     0.2337643E+03_r8,&  
         0.2326238E+03_r8,     0.2314684E+03_r8,     0.2302988E+03_r8,     0.2291154E+03_r8,     0.2279189E+03_r8,&  
         0.2267100E+03_r8,     0.2254897E+03_r8,     0.2242586E+03_r8,     0.2230179E+03_r8,     0.2217685E+03_r8,&  
         0.2205115E+03_r8,     0.2192477E+03_r8,     0.2179784E+03_r8,     0.2167045E+03_r8,     0.2154271E+03_r8,&  
         0.2141471E+03_r8,     0.2128655E+03_r8,     0.2115833E+03_r8,     0.2103011E+03_r8,     0.2090199E+03_r8,&  
         0.2077404E+03_r8,     0.2064633E+03_r8,     0.2051891E+03_r8,     0.2039185E+03_r8,     0.2026519E+03_r8,&  
         0.2013899E+03_r8,     0.2001327E+03_r8,     0.1988808E+03_r8,     0.1976345E+03_r8,     0.1963940E+03_r8,&  
         0.1951596E+03_r8,     0.1939315E+03_r8,     0.1927098E+03_r8,     0.1914948E+03_r8,     0.1902864E+03_r8/)  

    psaditmk(1:150,  107)= (/ &
         0.3160123E+03_r8,     0.3153843E+03_r8,     0.3147570E+03_r8,     0.3141306E+03_r8,     0.3135045E+03_r8,&  
         0.3128797E+03_r8,     0.3122556E+03_r8,     0.3116322E+03_r8,     0.3110094E+03_r8,     0.3103876E+03_r8,&  
         0.3097663E+03_r8,     0.3091454E+03_r8,     0.3085257E+03_r8,     0.3079062E+03_r8,     0.3072867E+03_r8,&  
         0.3066685E+03_r8,     0.3060503E+03_r8,     0.3054329E+03_r8,     0.3048157E+03_r8,     0.3041989E+03_r8,&  
         0.3035821E+03_r8,     0.3029664E+03_r8,     0.3023499E+03_r8,     0.3017346E+03_r8,     0.3011190E+03_r8,&  
         0.3005036E+03_r8,     0.2998883E+03_r8,     0.2992732E+03_r8,     0.2986578E+03_r8,     0.2980425E+03_r8,&  
         0.2974273E+03_r8,     0.2968117E+03_r8,     0.2961964E+03_r8,     0.2955803E+03_r8,     0.2949644E+03_r8,&  
         0.2943481E+03_r8,     0.2937313E+03_r8,     0.2931143E+03_r8,     0.2924967E+03_r8,     0.2918784E+03_r8,&  
         0.2912599E+03_r8,     0.2906405E+03_r8,     0.2900204E+03_r8,     0.2893999E+03_r8,     0.2887782E+03_r8,&  
         0.2881554E+03_r8,     0.2875325E+03_r8,     0.2869079E+03_r8,     0.2862822E+03_r8,     0.2856557E+03_r8,&  
         0.2850273E+03_r8,     0.2843982E+03_r8,     0.2837674E+03_r8,     0.2831351E+03_r8,     0.2825017E+03_r8,&  
         0.2818661E+03_r8,     0.2812288E+03_r8,     0.2805897E+03_r8,     0.2799486E+03_r8,     0.2793055E+03_r8,&  
         0.2786603E+03_r8,     0.2780127E+03_r8,     0.2773627E+03_r8,     0.2767102E+03_r8,     0.2760550E+03_r8,&  
         0.2753972E+03_r8,     0.2747365E+03_r8,     0.2740726E+03_r8,     0.2734057E+03_r8,     0.2727354E+03_r8,&  
         0.2720614E+03_r8,     0.2713839E+03_r8,     0.2707026E+03_r8,     0.2700173E+03_r8,     0.2693279E+03_r8,&  
         0.2686339E+03_r8,     0.2679355E+03_r8,     0.2672324E+03_r8,     0.2665242E+03_r8,     0.2658106E+03_r8,&  
         0.2650917E+03_r8,     0.2643676E+03_r8,     0.2636367E+03_r8,     0.2629001E+03_r8,     0.2621570E+03_r8,&  
         0.2614069E+03_r8,     0.2606502E+03_r8,     0.2598860E+03_r8,     0.2591140E+03_r8,     0.2583341E+03_r8,&  
         0.2575455E+03_r8,     0.2567487E+03_r8,     0.2559427E+03_r8,     0.2551273E+03_r8,     0.2543020E+03_r8,&  
         0.2534667E+03_r8,     0.2526208E+03_r8,     0.2517638E+03_r8,     0.2508954E+03_r8,     0.2500151E+03_r8,&  
         0.2491227E+03_r8,     0.2482177E+03_r8,     0.2472994E+03_r8,     0.2463678E+03_r8,     0.2454220E+03_r8,&  
         0.2444625E+03_r8,     0.2434880E+03_r8,     0.2424988E+03_r8,     0.2414942E+03_r8,     0.2404743E+03_r8,&  
         0.2394388E+03_r8,     0.2383874E+03_r8,     0.2373199E+03_r8,     0.2362368E+03_r8,     0.2351378E+03_r8,&  
         0.2340231E+03_r8,     0.2328931E+03_r8,     0.2317478E+03_r8,     0.2305880E+03_r8,     0.2294141E+03_r8,&  
         0.2282267E+03_r8,     0.2270264E+03_r8,     0.2258141E+03_r8,     0.2245906E+03_r8,     0.2233569E+03_r8,&  
         0.2221140E+03_r8,     0.2208627E+03_r8,     0.2196043E+03_r8,     0.2183396E+03_r8,     0.2170699E+03_r8,&  
         0.2157960E+03_r8,     0.2145190E+03_r8,     0.2132400E+03_r8,     0.2119597E+03_r8,     0.2106792E+03_r8,&  
         0.2093991E+03_r8,     0.2081204E+03_r8,     0.2068436E+03_r8,     0.2055695E+03_r8,     0.2042987E+03_r8,&  
         0.2030316E+03_r8,     0.2017688E+03_r8,     0.2005107E+03_r8,     0.1992577E+03_r8,     0.1980101E+03_r8,&  
         0.1967682E+03_r8,     0.1955323E+03_r8,     0.1943025E+03_r8,     0.1930791E+03_r8,     0.1918622E+03_r8/)  

    psaditmk(1:150,  108)= (/ &
         0.3163077E+03_r8,     0.3156807E+03_r8,     0.3150552E+03_r8,     0.3144307E+03_r8,     0.3138067E+03_r8,&  
         0.3131835E+03_r8,     0.3125612E+03_r8,     0.3119397E+03_r8,     0.3113191E+03_r8,     0.3106985E+03_r8,&  
         0.3100795E+03_r8,     0.3094609E+03_r8,     0.3088427E+03_r8,     0.3082252E+03_r8,     0.3076084E+03_r8,&  
         0.3069920E+03_r8,     0.3063760E+03_r8,     0.3057608E+03_r8,     0.3051458E+03_r8,     0.3045313E+03_r8,&  
         0.3039175E+03_r8,     0.3033031E+03_r8,     0.3026901E+03_r8,     0.3020770E+03_r8,     0.3014637E+03_r8,&  
         0.3008512E+03_r8,     0.3002383E+03_r8,     0.2996259E+03_r8,     0.2990135E+03_r8,     0.2984009E+03_r8,&  
         0.2977886E+03_r8,     0.2971762E+03_r8,     0.2965635E+03_r8,     0.2959512E+03_r8,     0.2953382E+03_r8,&  
         0.2947248E+03_r8,     0.2941117E+03_r8,     0.2934978E+03_r8,     0.2928837E+03_r8,     0.2922693E+03_r8,&  
         0.2916541E+03_r8,     0.2910384E+03_r8,     0.2904225E+03_r8,     0.2898054E+03_r8,     0.2891876E+03_r8,&  
         0.2885696E+03_r8,     0.2879503E+03_r8,     0.2873296E+03_r8,     0.2867089E+03_r8,     0.2860867E+03_r8,&  
         0.2854632E+03_r8,     0.2848386E+03_r8,     0.2842126E+03_r8,     0.2835854E+03_r8,     0.2829572E+03_r8,&  
         0.2823264E+03_r8,     0.2816948E+03_r8,     0.2810613E+03_r8,     0.2804259E+03_r8,     0.2797885E+03_r8,&  
         0.2791497E+03_r8,     0.2785082E+03_r8,     0.2778645E+03_r8,     0.2772186E+03_r8,     0.2765702E+03_r8,&  
         0.2759193E+03_r8,     0.2752656E+03_r8,     0.2746091E+03_r8,     0.2739497E+03_r8,     0.2732873E+03_r8,&  
         0.2726213E+03_r8,     0.2719521E+03_r8,     0.2712796E+03_r8,     0.2706029E+03_r8,     0.2699225E+03_r8,&  
         0.2692381E+03_r8,     0.2685496E+03_r8,     0.2678561E+03_r8,     0.2671585E+03_r8,     0.2664557E+03_r8,&  
         0.2657478E+03_r8,     0.2650345E+03_r8,     0.2643160E+03_r8,     0.2635914E+03_r8,     0.2628611E+03_r8,&  
         0.2621242E+03_r8,     0.2613806E+03_r8,     0.2606305E+03_r8,     0.2598730E+03_r8,     0.2591081E+03_r8,&  
         0.2583354E+03_r8,     0.2575544E+03_r8,     0.2567647E+03_r8,     0.2559668E+03_r8,     0.2551592E+03_r8,&  
         0.2543421E+03_r8,     0.2535153E+03_r8,     0.2526777E+03_r8,     0.2518295E+03_r8,     0.2509701E+03_r8,&  
         0.2500991E+03_r8,     0.2492160E+03_r8,     0.2483207E+03_r8,     0.2474125E+03_r8,     0.2464907E+03_r8,&  
         0.2455555E+03_r8,     0.2446059E+03_r8,     0.2436424E+03_r8,     0.2426637E+03_r8,     0.2416701E+03_r8,&  
         0.2406611E+03_r8,     0.2396364E+03_r8,     0.2385961E+03_r8,     0.2375400E+03_r8,     0.2364677E+03_r8,&  
         0.2353796E+03_r8,     0.2342757E+03_r8,     0.2331562E+03_r8,     0.2320213E+03_r8,     0.2308716E+03_r8,&  
         0.2297073E+03_r8,     0.2285292E+03_r8,     0.2273378E+03_r8,     0.2261338E+03_r8,     0.2249183E+03_r8,&  
         0.2236918E+03_r8,     0.2224556E+03_r8,     0.2212106E+03_r8,     0.2199577E+03_r8,     0.2186980E+03_r8,&  
         0.2174327E+03_r8,     0.2161627E+03_r8,     0.2148890E+03_r8,     0.2136127E+03_r8,     0.2123347E+03_r8,&  
         0.2110559E+03_r8,     0.2097773E+03_r8,     0.2084995E+03_r8,     0.2072233E+03_r8,     0.2059494E+03_r8,&  
         0.2046785E+03_r8,     0.2034111E+03_r8,     0.2021477E+03_r8,     0.2008887E+03_r8,     0.1996347E+03_r8,&  
         0.1983859E+03_r8,     0.1971427E+03_r8,     0.1959053E+03_r8,     0.1946739E+03_r8,     0.1934488E+03_r8/)  

    psaditmk(1:150,  109)= (/ &
         0.3165989E+03_r8,     0.3159736E+03_r8,     0.3153493E+03_r8,     0.3147262E+03_r8,     0.3141041E+03_r8,&  
         0.3134830E+03_r8,     0.3128623E+03_r8,     0.3122425E+03_r8,     0.3116233E+03_r8,     0.3110055E+03_r8,&  
         0.3103882E+03_r8,     0.3097708E+03_r8,     0.3091550E+03_r8,     0.3085395E+03_r8,     0.3079247E+03_r8,&  
         0.3073105E+03_r8,     0.3066969E+03_r8,     0.3060834E+03_r8,     0.3054710E+03_r8,     0.3048587E+03_r8,&  
         0.3042465E+03_r8,     0.3036353E+03_r8,     0.3030241E+03_r8,     0.3024133E+03_r8,     0.3018029E+03_r8,&  
         0.3011926E+03_r8,     0.3005826E+03_r8,     0.2999726E+03_r8,     0.2993629E+03_r8,     0.2987535E+03_r8,&  
         0.2981440E+03_r8,     0.2975340E+03_r8,     0.2969245E+03_r8,     0.2963148E+03_r8,     0.2957047E+03_r8,&  
         0.2950950E+03_r8,     0.2944847E+03_r8,     0.2938741E+03_r8,     0.2932637E+03_r8,     0.2926523E+03_r8,&  
         0.2920406E+03_r8,     0.2914289E+03_r8,     0.2908161E+03_r8,     0.2902028E+03_r8,     0.2895892E+03_r8,&  
         0.2889748E+03_r8,     0.2883593E+03_r8,     0.2877437E+03_r8,     0.2871267E+03_r8,     0.2865085E+03_r8,&  
         0.2858894E+03_r8,     0.2852697E+03_r8,     0.2846482E+03_r8,     0.2840257E+03_r8,     0.2834025E+03_r8,&  
         0.2827771E+03_r8,     0.2821504E+03_r8,     0.2815222E+03_r8,     0.2808922E+03_r8,     0.2802603E+03_r8,&  
         0.2796271E+03_r8,     0.2789917E+03_r8,     0.2783541E+03_r8,     0.2777145E+03_r8,     0.2770725E+03_r8,&  
         0.2764282E+03_r8,     0.2757813E+03_r8,     0.2751319E+03_r8,     0.2744795E+03_r8,     0.2738244E+03_r8,&  
         0.2731665E+03_r8,     0.2725052E+03_r8,     0.2718405E+03_r8,     0.2711725E+03_r8,     0.2705007E+03_r8,&  
         0.2698250E+03_r8,     0.2691456E+03_r8,     0.2684620E+03_r8,     0.2677742E+03_r8,     0.2670813E+03_r8,&  
         0.2663843E+03_r8,     0.2656818E+03_r8,     0.2649742E+03_r8,     0.2642614E+03_r8,     0.2635430E+03_r8,&  
         0.2628185E+03_r8,     0.2620879E+03_r8,     0.2613507E+03_r8,     0.2606070E+03_r8,     0.2598563E+03_r8,&  
         0.2590982E+03_r8,     0.2583325E+03_r8,     0.2575589E+03_r8,     0.2567767E+03_r8,     0.2559860E+03_r8,&  
         0.2551864E+03_r8,     0.2543773E+03_r8,     0.2535586E+03_r8,     0.2527296E+03_r8,     0.2518901E+03_r8,&  
         0.2510393E+03_r8,     0.2501776E+03_r8,     0.2493039E+03_r8,     0.2484182E+03_r8,     0.2475194E+03_r8,&  
         0.2466078E+03_r8,     0.2456826E+03_r8,     0.2447436E+03_r8,     0.2437905E+03_r8,     0.2428226E+03_r8,&  
         0.2418396E+03_r8,     0.2408415E+03_r8,     0.2398279E+03_r8,     0.2387986E+03_r8,     0.2377534E+03_r8,&  
         0.2366922E+03_r8,     0.2356151E+03_r8,     0.2345220E+03_r8,     0.2334132E+03_r8,     0.2322889E+03_r8,&  
         0.2311492E+03_r8,     0.2299949E+03_r8,     0.2288262E+03_r8,     0.2276439E+03_r8,     0.2264486E+03_r8,&  
         0.2252412E+03_r8,     0.2240224E+03_r8,     0.2227932E+03_r8,     0.2215547E+03_r8,     0.2203077E+03_r8,&  
         0.2190534E+03_r8,     0.2177928E+03_r8,     0.2165269E+03_r8,     0.2152569E+03_r8,     0.2139836E+03_r8,&  
         0.2127081E+03_r8,     0.2114314E+03_r8,     0.2101543E+03_r8,     0.2088776E+03_r8,     0.2076022E+03_r8,&  
         0.2063288E+03_r8,     0.2050579E+03_r8,     0.2037903E+03_r8,     0.2025264E+03_r8,     0.2012668E+03_r8,&  
         0.2000118E+03_r8,     0.1987619E+03_r8,     0.1975174E+03_r8,     0.1962786E+03_r8,     0.1950457E+03_r8/)  

    psaditmk(1:150,  110)= (/ &
         0.3168859E+03_r8,     0.3162622E+03_r8,     0.3156398E+03_r8,     0.3150183E+03_r8,     0.3143974E+03_r8,&  
         0.3137777E+03_r8,     0.3131589E+03_r8,     0.3125410E+03_r8,     0.3119236E+03_r8,     0.3113070E+03_r8,&  
         0.3106920E+03_r8,     0.3100770E+03_r8,     0.3094625E+03_r8,     0.3088487E+03_r8,     0.3082362E+03_r8,&  
         0.3076240E+03_r8,     0.3070121E+03_r8,     0.3064012E+03_r8,     0.3057904E+03_r8,     0.3051803E+03_r8,&  
         0.3045711E+03_r8,     0.3039615E+03_r8,     0.3033528E+03_r8,     0.3027445E+03_r8,     0.3021361E+03_r8,&  
         0.3015287E+03_r8,     0.3009211E+03_r8,     0.3003135E+03_r8,     0.2997066E+03_r8,     0.2990996E+03_r8,&  
         0.2984925E+03_r8,     0.2978856E+03_r8,     0.2972786E+03_r8,     0.2966721E+03_r8,     0.2960652E+03_r8,&  
         0.2954580E+03_r8,     0.2948509E+03_r8,     0.2942438E+03_r8,     0.2936362E+03_r8,     0.2930281E+03_r8,&  
         0.2924203E+03_r8,     0.2918114E+03_r8,     0.2912024E+03_r8,     0.2905932E+03_r8,     0.2899828E+03_r8,&  
         0.2893718E+03_r8,     0.2887609E+03_r8,     0.2881486E+03_r8,     0.2875354E+03_r8,     0.2869218E+03_r8,&  
         0.2863074E+03_r8,     0.2856916E+03_r8,     0.2850750E+03_r8,     0.2844570E+03_r8,     0.2838379E+03_r8,&  
         0.2832174E+03_r8,     0.2825958E+03_r8,     0.2819729E+03_r8,     0.2813481E+03_r8,     0.2807217E+03_r8,&  
         0.2800937E+03_r8,     0.2794638E+03_r8,     0.2788322E+03_r8,     0.2781986E+03_r8,     0.2775627E+03_r8,&  
         0.2769245E+03_r8,     0.2762842E+03_r8,     0.2756415E+03_r8,     0.2749962E+03_r8,     0.2743481E+03_r8,&  
         0.2736973E+03_r8,     0.2730435E+03_r8,     0.2723868E+03_r8,     0.2717267E+03_r8,     0.2710629E+03_r8,&  
         0.2703958E+03_r8,     0.2697253E+03_r8,     0.2690504E+03_r8,     0.2683719E+03_r8,     0.2676891E+03_r8,&  
         0.2670016E+03_r8,     0.2663095E+03_r8,     0.2656128E+03_r8,     0.2649109E+03_r8,     0.2642034E+03_r8,&  
         0.2634911E+03_r8,     0.2627725E+03_r8,     0.2620482E+03_r8,     0.2613172E+03_r8,     0.2605799E+03_r8,&  
         0.2598357E+03_r8,     0.2590844E+03_r8,     0.2583256E+03_r8,     0.2575591E+03_r8,     0.2567842E+03_r8,&  
         0.2560010E+03_r8,     0.2552090E+03_r8,     0.2544078E+03_r8,     0.2535972E+03_r8,     0.2527765E+03_r8,&  
         0.2519455E+03_r8,     0.2511038E+03_r8,     0.2502506E+03_r8,     0.2493862E+03_r8,     0.2485097E+03_r8,&  
         0.2476209E+03_r8,     0.2467189E+03_r8,     0.2458038E+03_r8,     0.2448750E+03_r8,     0.2439321E+03_r8,&  
         0.2429749E+03_r8,     0.2420027E+03_r8,     0.2410155E+03_r8,     0.2400128E+03_r8,     0.2389944E+03_r8,&  
         0.2379604E+03_r8,     0.2369103E+03_r8,     0.2358441E+03_r8,     0.2347620E+03_r8,     0.2336640E+03_r8,&  
         0.2325502E+03_r8,     0.2314209E+03_r8,     0.2302766E+03_r8,     0.2291177E+03_r8,     0.2279447E+03_r8,&  
         0.2267584E+03_r8,     0.2255594E+03_r8,     0.2243485E+03_r8,     0.2231268E+03_r8,     0.2218950E+03_r8,&  
         0.2206543E+03_r8,     0.2194055E+03_r8,     0.2181500E+03_r8,     0.2168885E+03_r8,     0.2156224E+03_r8,&  
         0.2143524E+03_r8,     0.2130798E+03_r8,     0.2118053E+03_r8,     0.2105300E+03_r8,     0.2092548E+03_r8,&  
         0.2079803E+03_r8,     0.2067075E+03_r8,     0.2054369E+03_r8,     0.2041692E+03_r8,     0.2029050E+03_r8,&  
         0.2016447E+03_r8,     0.2003889E+03_r8,     0.1991380E+03_r8,     0.1978923E+03_r8,     0.1966521E+03_r8/)  

    psaditmk(1:150,  111)= (/ &
         0.3171682E+03_r8,     0.3165469E+03_r8,     0.3159259E+03_r8,     0.3153057E+03_r8,     0.3146868E+03_r8,&  
         0.3140685E+03_r8,     0.3134512E+03_r8,     0.3128349E+03_r8,     0.3122197E+03_r8,     0.3116051E+03_r8,&  
         0.3109906E+03_r8,     0.3103778E+03_r8,     0.3097657E+03_r8,     0.3091538E+03_r8,     0.3085429E+03_r8,&  
         0.3079328E+03_r8,     0.3073231E+03_r8,     0.3067139E+03_r8,     0.3061053E+03_r8,     0.3054974E+03_r8,&  
         0.3048898E+03_r8,     0.3042831E+03_r8,     0.3036761E+03_r8,     0.3030698E+03_r8,     0.3024645E+03_r8,&  
         0.3018586E+03_r8,     0.3012536E+03_r8,     0.3006491E+03_r8,     0.3000442E+03_r8,     0.2994397E+03_r8,&  
         0.2988357E+03_r8,     0.2982311E+03_r8,     0.2976269E+03_r8,     0.2970231E+03_r8,     0.2964190E+03_r8,&  
         0.2958149E+03_r8,     0.2952108E+03_r8,     0.2946064E+03_r8,     0.2940020E+03_r8,     0.2933973E+03_r8,&  
         0.2927923E+03_r8,     0.2921869E+03_r8,     0.2915817E+03_r8,     0.2909752E+03_r8,     0.2903686E+03_r8,&  
         0.2897618E+03_r8,     0.2891539E+03_r8,     0.2885454E+03_r8,     0.2879367E+03_r8,     0.2873269E+03_r8,&  
         0.2867163E+03_r8,     0.2861048E+03_r8,     0.2854924E+03_r8,     0.2848788E+03_r8,     0.2842644E+03_r8,&  
         0.2836483E+03_r8,     0.2830313E+03_r8,     0.2824134E+03_r8,     0.2817935E+03_r8,     0.2811723E+03_r8,&  
         0.2805497E+03_r8,     0.2799250E+03_r8,     0.2792994E+03_r8,     0.2786713E+03_r8,     0.2780411E+03_r8,&  
         0.2774091E+03_r8,     0.2767749E+03_r8,     0.2761385E+03_r8,     0.2754998E+03_r8,     0.2748585E+03_r8,&  
         0.2742147E+03_r8,     0.2735680E+03_r8,     0.2729185E+03_r8,     0.2722659E+03_r8,     0.2716104E+03_r8,&  
         0.2709513E+03_r8,     0.2702887E+03_r8,     0.2696228E+03_r8,     0.2689531E+03_r8,     0.2682792E+03_r8,&  
         0.2676013E+03_r8,     0.2669192E+03_r8,     0.2662322E+03_r8,     0.2655406E+03_r8,     0.2648446E+03_r8,&  
         0.2641426E+03_r8,     0.2634358E+03_r8,     0.2627233E+03_r8,     0.2620047E+03_r8,     0.2612803E+03_r8,&  
         0.2605491E+03_r8,     0.2598113E+03_r8,     0.2590666E+03_r8,     0.2583147E+03_r8,     0.2575552E+03_r8,&  
         0.2567874E+03_r8,     0.2560116E+03_r8,     0.2552272E+03_r8,     0.2544340E+03_r8,     0.2536310E+03_r8,&  
         0.2528186E+03_r8,     0.2519959E+03_r8,     0.2511627E+03_r8,     0.2503187E+03_r8,     0.2494632E+03_r8,&  
         0.2485959E+03_r8,     0.2477164E+03_r8,     0.2468243E+03_r8,     0.2459191E+03_r8,     0.2450003E+03_r8,&  
         0.2440679E+03_r8,     0.2431209E+03_r8,     0.2421595E+03_r8,     0.2411830E+03_r8,     0.2401912E+03_r8,&  
         0.2391840E+03_r8,     0.2381608E+03_r8,     0.2371218E+03_r8,     0.2360667E+03_r8,     0.2349956E+03_r8,&  
         0.2339084E+03_r8,     0.2328053E+03_r8,     0.2316866E+03_r8,     0.2305525E+03_r8,     0.2294036E+03_r8,&  
         0.2282402E+03_r8,     0.2270629E+03_r8,     0.2258727E+03_r8,     0.2246700E+03_r8,     0.2234559E+03_r8,&  
         0.2222312E+03_r8,     0.2209971E+03_r8,     0.2197543E+03_r8,     0.2185041E+03_r8,     0.2172474E+03_r8,&  
         0.2159854E+03_r8,     0.2147192E+03_r8,     0.2134496E+03_r8,     0.2121777E+03_r8,     0.2109044E+03_r8,&  
         0.2096308E+03_r8,     0.2083575E+03_r8,     0.2070854E+03_r8,     0.2058153E+03_r8,     0.2045476E+03_r8,&  
         0.2032832E+03_r8,     0.2020225E+03_r8,     0.2007660E+03_r8,     0.1995142E+03_r8,     0.1982673E+03_r8/)  

    psaditmk(1:150,  112)= (/ &
         0.3174477E+03_r8,     0.3168271E+03_r8,     0.3162078E+03_r8,     0.3155891E+03_r8,     0.3149721E+03_r8,&  
         0.3143554E+03_r8,     0.3137397E+03_r8,     0.3131248E+03_r8,     0.3125110E+03_r8,     0.3118983E+03_r8,&  
         0.3112863E+03_r8,     0.3106743E+03_r8,     0.3100641E+03_r8,     0.3094543E+03_r8,     0.3088452E+03_r8,&  
         0.3082367E+03_r8,     0.3076290E+03_r8,     0.3070221E+03_r8,     0.3064149E+03_r8,     0.3058093E+03_r8,&  
         0.3052039E+03_r8,     0.3045989E+03_r8,     0.3039946E+03_r8,     0.3033908E+03_r8,     0.3027870E+03_r8,&  
         0.3021839E+03_r8,     0.3015811E+03_r8,     0.3009783E+03_r8,     0.3003763E+03_r8,     0.2997744E+03_r8,&  
         0.2991724E+03_r8,     0.2985709E+03_r8,     0.2979697E+03_r8,     0.2973681E+03_r8,     0.2967668E+03_r8,&  
         0.2961657E+03_r8,     0.2955642E+03_r8,     0.2949627E+03_r8,     0.2943612E+03_r8,     0.2937595E+03_r8,&  
         0.2931577E+03_r8,     0.2925557E+03_r8,     0.2919530E+03_r8,     0.2913504E+03_r8,     0.2907477E+03_r8,&  
         0.2901437E+03_r8,     0.2895394E+03_r8,     0.2889350E+03_r8,     0.2883296E+03_r8,     0.2877237E+03_r8,&  
         0.2871168E+03_r8,     0.2865095E+03_r8,     0.2859012E+03_r8,     0.2852918E+03_r8,     0.2846818E+03_r8,&  
         0.2840701E+03_r8,     0.2834576E+03_r8,     0.2828446E+03_r8,     0.2822292E+03_r8,     0.2816130E+03_r8,&  
         0.2809954E+03_r8,     0.2803761E+03_r8,     0.2797553E+03_r8,     0.2791329E+03_r8,     0.2785083E+03_r8,&  
         0.2778821E+03_r8,     0.2772539E+03_r8,     0.2766236E+03_r8,     0.2759910E+03_r8,     0.2753563E+03_r8,&  
         0.2747190E+03_r8,     0.2740791E+03_r8,     0.2734366E+03_r8,     0.2727914E+03_r8,     0.2721433E+03_r8,&  
         0.2714918E+03_r8,     0.2708371E+03_r8,     0.2701794E+03_r8,     0.2695178E+03_r8,     0.2688528E+03_r8,&  
         0.2681840E+03_r8,     0.2675110E+03_r8,     0.2668336E+03_r8,     0.2661522E+03_r8,     0.2654658E+03_r8,&  
         0.2647749E+03_r8,     0.2640789E+03_r8,     0.2633774E+03_r8,     0.2626706E+03_r8,     0.2619580E+03_r8,&  
         0.2612394E+03_r8,     0.2605148E+03_r8,     0.2597832E+03_r8,     0.2590450E+03_r8,     0.2582998E+03_r8,&  
         0.2575472E+03_r8,     0.2567866E+03_r8,     0.2560180E+03_r8,     0.2552410E+03_r8,     0.2544552E+03_r8,&  
         0.2536602E+03_r8,     0.2528558E+03_r8,     0.2520415E+03_r8,     0.2512167E+03_r8,     0.2503813E+03_r8,&  
         0.2495347E+03_r8,     0.2486766E+03_r8,     0.2478064E+03_r8,     0.2469239E+03_r8,     0.2460284E+03_r8,&  
         0.2451200E+03_r8,     0.2441973E+03_r8,     0.2432610E+03_r8,     0.2423100E+03_r8,     0.2413443E+03_r8,&  
         0.2403634E+03_r8,     0.2393669E+03_r8,     0.2383550E+03_r8,     0.2373269E+03_r8,     0.2362829E+03_r8,&  
         0.2352228E+03_r8,     0.2341465E+03_r8,     0.2330543E+03_r8,     0.2319462E+03_r8,     0.2308226E+03_r8,&  
         0.2296836E+03_r8,     0.2285300E+03_r8,     0.2273622E+03_r8,     0.2261809E+03_r8,     0.2249867E+03_r8,&  
         0.2237806E+03_r8,     0.2225633E+03_r8,     0.2213360E+03_r8,     0.2200996E+03_r8,     0.2188550E+03_r8,&  
         0.2176035E+03_r8,     0.2163459E+03_r8,     0.2150836E+03_r8,     0.2138173E+03_r8,     0.2125483E+03_r8,&  
         0.2112774E+03_r8,     0.2100056E+03_r8,     0.2087337E+03_r8,     0.2074626E+03_r8,     0.2061930E+03_r8,&  
         0.2049257E+03_r8,     0.2036612E+03_r8,     0.2024001E+03_r8,     0.2011431E+03_r8,     0.1998904E+03_r8/)  

    psaditmk(1:150,  113)= (/ &
         0.3177231E+03_r8,     0.3171039E+03_r8,     0.3164858E+03_r8,     0.3158689E+03_r8,     0.3152529E+03_r8,&  
         0.3146383E+03_r8,     0.3140239E+03_r8,     0.3134109E+03_r8,     0.3127985E+03_r8,     0.3121872E+03_r8,&  
         0.3115768E+03_r8,     0.3109674E+03_r8,     0.3103580E+03_r8,     0.3097498E+03_r8,     0.3091430E+03_r8,&  
         0.3085363E+03_r8,     0.3079302E+03_r8,     0.3073249E+03_r8,     0.3067206E+03_r8,     0.3061164E+03_r8,&  
         0.3055128E+03_r8,     0.3049101E+03_r8,     0.3043076E+03_r8,     0.3037058E+03_r8,     0.3031046E+03_r8,&  
         0.3025036E+03_r8,     0.3019026E+03_r8,     0.3013027E+03_r8,     0.3007028E+03_r8,     0.3001031E+03_r8,&  
         0.2995040E+03_r8,     0.2989050E+03_r8,     0.2983059E+03_r8,     0.2977071E+03_r8,     0.2971085E+03_r8,&  
         0.2965099E+03_r8,     0.2959113E+03_r8,     0.2953129E+03_r8,     0.2947139E+03_r8,     0.2941154E+03_r8,&  
         0.2935165E+03_r8,     0.2929173E+03_r8,     0.2923183E+03_r8,     0.2917186E+03_r8,     0.2911186E+03_r8,&  
         0.2905183E+03_r8,     0.2899178E+03_r8,     0.2893167E+03_r8,     0.2887149E+03_r8,     0.2881127E+03_r8,&  
         0.2875098E+03_r8,     0.2869060E+03_r8,     0.2863017E+03_r8,     0.2856966E+03_r8,     0.2850905E+03_r8,&  
         0.2844833E+03_r8,     0.2838752E+03_r8,     0.2832661E+03_r8,     0.2826556E+03_r8,     0.2820441E+03_r8,&  
         0.2814312E+03_r8,     0.2808169E+03_r8,     0.2802012E+03_r8,     0.2795840E+03_r8,     0.2789648E+03_r8,&  
         0.2783440E+03_r8,     0.2777214E+03_r8,     0.2770970E+03_r8,     0.2764704E+03_r8,     0.2758416E+03_r8,&  
         0.2752108E+03_r8,     0.2745775E+03_r8,     0.2739417E+03_r8,     0.2733035E+03_r8,     0.2726621E+03_r8,&  
         0.2720182E+03_r8,     0.2713712E+03_r8,     0.2707211E+03_r8,     0.2700677E+03_r8,     0.2694108E+03_r8,&  
         0.2687503E+03_r8,     0.2680861E+03_r8,     0.2674181E+03_r8,     0.2667457E+03_r8,     0.2660692E+03_r8,&  
         0.2653883E+03_r8,     0.2647021E+03_r8,     0.2640118E+03_r8,     0.2633158E+03_r8,     0.2626146E+03_r8,&  
         0.2619080E+03_r8,     0.2611954E+03_r8,     0.2604769E+03_r8,     0.2597516E+03_r8,     0.2590199E+03_r8,&  
         0.2582812E+03_r8,     0.2575352E+03_r8,     0.2567816E+03_r8,     0.2560201E+03_r8,     0.2552505E+03_r8,&  
         0.2544722E+03_r8,     0.2536849E+03_r8,     0.2528884E+03_r8,     0.2520821E+03_r8,     0.2512659E+03_r8,&  
         0.2504389E+03_r8,     0.2496011E+03_r8,     0.2487519E+03_r8,     0.2478911E+03_r8,     0.2470180E+03_r8,&  
         0.2461321E+03_r8,     0.2452331E+03_r8,     0.2443211E+03_r8,     0.2433947E+03_r8,     0.2424543E+03_r8,&  
         0.2414993E+03_r8,     0.2405291E+03_r8,     0.2395435E+03_r8,     0.2385425E+03_r8,     0.2375255E+03_r8,&  
         0.2364926E+03_r8,     0.2354435E+03_r8,     0.2343783E+03_r8,     0.2332970E+03_r8,     0.2321997E+03_r8,&  
         0.2310864E+03_r8,     0.2299579E+03_r8,     0.2288143E+03_r8,     0.2276561E+03_r8,     0.2264839E+03_r8,&  
         0.2252986E+03_r8,     0.2241006E+03_r8,     0.2228912E+03_r8,     0.2216709E+03_r8,     0.2204411E+03_r8,&  
         0.2192026E+03_r8,     0.2179564E+03_r8,     0.2167037E+03_r8,     0.2154455E+03_r8,     0.2141830E+03_r8,&  
         0.2129170E+03_r8,     0.2116487E+03_r8,     0.2103790E+03_r8,     0.2091087E+03_r8,     0.2078388E+03_r8,&  
         0.2065700E+03_r8,     0.2053031E+03_r8,     0.2040387E+03_r8,     0.2027775E+03_r8,     0.2015199E+03_r8/)  

    psaditmk(1:150,  114)= (/ &
         0.3179943E+03_r8,     0.3173770E+03_r8,     0.3167603E+03_r8,     0.3161448E+03_r8,     0.3155302E+03_r8,&  
         0.3149166E+03_r8,     0.3143040E+03_r8,     0.3136926E+03_r8,     0.3130822E+03_r8,     0.3124723E+03_r8,&  
         0.3118633E+03_r8,     0.3112552E+03_r8,     0.3106482E+03_r8,     0.3100416E+03_r8,     0.3094358E+03_r8,&  
         0.3088313E+03_r8,     0.3082275E+03_r8,     0.3076238E+03_r8,     0.3070206E+03_r8,     0.3064189E+03_r8,&  
         0.3058176E+03_r8,     0.3052163E+03_r8,     0.3046164E+03_r8,     0.3040167E+03_r8,     0.3034169E+03_r8,&  
         0.3028186E+03_r8,     0.3022201E+03_r8,     0.3016216E+03_r8,     0.3010242E+03_r8,     0.3004270E+03_r8,&  
         0.2998297E+03_r8,     0.2992331E+03_r8,     0.2986370E+03_r8,     0.2980403E+03_r8,     0.2974441E+03_r8,&  
         0.2968485E+03_r8,     0.2962524E+03_r8,     0.2956566E+03_r8,     0.2950606E+03_r8,     0.2944646E+03_r8,&  
         0.2938686E+03_r8,     0.2932730E+03_r8,     0.2926766E+03_r8,     0.2920798E+03_r8,     0.2914832E+03_r8,&  
         0.2908863E+03_r8,     0.2902887E+03_r8,     0.2896909E+03_r8,     0.2890928E+03_r8,     0.2884942E+03_r8,&  
         0.2878949E+03_r8,     0.2872947E+03_r8,     0.2866944E+03_r8,     0.2860928E+03_r8,     0.2854907E+03_r8,&  
         0.2848881E+03_r8,     0.2842839E+03_r8,     0.2836790E+03_r8,     0.2830730E+03_r8,     0.2824658E+03_r8,&  
         0.2818578E+03_r8,     0.2812480E+03_r8,     0.2806374E+03_r8,     0.2800250E+03_r8,     0.2794109E+03_r8,&  
         0.2787955E+03_r8,     0.2781782E+03_r8,     0.2775596E+03_r8,     0.2769384E+03_r8,     0.2763155E+03_r8,&  
         0.2756908E+03_r8,     0.2750636E+03_r8,     0.2744342E+03_r8,     0.2738024E+03_r8,     0.2731680E+03_r8,&  
         0.2725311E+03_r8,     0.2718914E+03_r8,     0.2712484E+03_r8,     0.2706027E+03_r8,     0.2699538E+03_r8,&  
         0.2693014E+03_r8,     0.2686456E+03_r8,     0.2679859E+03_r8,     0.2673225E+03_r8,     0.2666553E+03_r8,&  
         0.2659835E+03_r8,     0.2653078E+03_r8,     0.2646271E+03_r8,     0.2639418E+03_r8,     0.2632516E+03_r8,&  
         0.2625558E+03_r8,     0.2618549E+03_r8,     0.2611481E+03_r8,     0.2604354E+03_r8,     0.2597165E+03_r8,&  
         0.2589910E+03_r8,     0.2582589E+03_r8,     0.2575197E+03_r8,     0.2567727E+03_r8,     0.2560182E+03_r8,&  
         0.2552557E+03_r8,     0.2544849E+03_r8,     0.2537052E+03_r8,     0.2529164E+03_r8,     0.2521181E+03_r8,&  
         0.2513099E+03_r8,     0.2504916E+03_r8,     0.2496624E+03_r8,     0.2488222E+03_r8,     0.2479702E+03_r8,&  
         0.2471065E+03_r8,     0.2462301E+03_r8,     0.2453413E+03_r8,     0.2444388E+03_r8,     0.2435227E+03_r8,&  
         0.2425927E+03_r8,     0.2416480E+03_r8,     0.2406885E+03_r8,     0.2397139E+03_r8,     0.2387237E+03_r8,&  
         0.2377178E+03_r8,     0.2366959E+03_r8,     0.2356579E+03_r8,     0.2346037E+03_r8,     0.2335334E+03_r8,&  
         0.2324468E+03_r8,     0.2313444E+03_r8,     0.2302262E+03_r8,     0.2290928E+03_r8,     0.2279444E+03_r8,&  
         0.2267817E+03_r8,     0.2256053E+03_r8,     0.2244159E+03_r8,     0.2232145E+03_r8,     0.2220018E+03_r8,&  
         0.2207787E+03_r8,     0.2195466E+03_r8,     0.2183062E+03_r8,     0.2170586E+03_r8,     0.2158050E+03_r8,&  
         0.2145463E+03_r8,     0.2132838E+03_r8,     0.2120184E+03_r8,     0.2107510E+03_r8,     0.2094826E+03_r8,&  
         0.2082141E+03_r8,     0.2069463E+03_r8,     0.2056800E+03_r8,     0.2044159E+03_r8,     0.2031546E+03_r8/)  

    psaditmk(1:150,  115)= (/ &
         0.3182625E+03_r8,     0.3176460E+03_r8,     0.3170309E+03_r8,     0.3164169E+03_r8,     0.3158037E+03_r8,&  
         0.3151916E+03_r8,     0.3145804E+03_r8,     0.3139705E+03_r8,     0.3133615E+03_r8,     0.3127534E+03_r8,&  
         0.3121458E+03_r8,     0.3115393E+03_r8,     0.3109337E+03_r8,     0.3103293E+03_r8,     0.3097254E+03_r8,&  
         0.3091220E+03_r8,     0.3085197E+03_r8,     0.3079180E+03_r8,     0.3073174E+03_r8,     0.3067166E+03_r8,&  
         0.3061174E+03_r8,     0.3055185E+03_r8,     0.3049198E+03_r8,     0.3043221E+03_r8,     0.3037248E+03_r8,&  
         0.3031279E+03_r8,     0.3025316E+03_r8,     0.3019357E+03_r8,     0.3013401E+03_r8,     0.3007452E+03_r8,&  
         0.3001506E+03_r8,     0.2995562E+03_r8,     0.2989620E+03_r8,     0.2983683E+03_r8,     0.2977749E+03_r8,&  
         0.2971811E+03_r8,     0.2965879E+03_r8,     0.2959946E+03_r8,     0.2954012E+03_r8,     0.2948079E+03_r8,&  
         0.2942150E+03_r8,     0.2936214E+03_r8,     0.2930280E+03_r8,     0.2924348E+03_r8,     0.2918411E+03_r8,&  
         0.2912469E+03_r8,     0.2906529E+03_r8,     0.2900585E+03_r8,     0.2894633E+03_r8,     0.2888681E+03_r8,&  
         0.2882723E+03_r8,     0.2876761E+03_r8,     0.2870789E+03_r8,     0.2864814E+03_r8,     0.2858831E+03_r8,&  
         0.2852842E+03_r8,     0.2846841E+03_r8,     0.2840833E+03_r8,     0.2834812E+03_r8,     0.2828784E+03_r8,&  
         0.2822751E+03_r8,     0.2816699E+03_r8,     0.2810639E+03_r8,     0.2804561E+03_r8,     0.2798471E+03_r8,&  
         0.2792367E+03_r8,     0.2786245E+03_r8,     0.2780109E+03_r8,     0.2773955E+03_r8,     0.2767783E+03_r8,&  
         0.2761593E+03_r8,     0.2755379E+03_r8,     0.2749147E+03_r8,     0.2742891E+03_r8,     0.2736612E+03_r8,&  
         0.2730308E+03_r8,     0.2723978E+03_r8,     0.2717622E+03_r8,     0.2711237E+03_r8,     0.2704818E+03_r8,&  
         0.2698375E+03_r8,     0.2691897E+03_r8,     0.2685382E+03_r8,     0.2678835E+03_r8,     0.2672246E+03_r8,&  
         0.2665622E+03_r8,     0.2658955E+03_r8,     0.2652242E+03_r8,     0.2645491E+03_r8,     0.2638690E+03_r8,&  
         0.2631840E+03_r8,     0.2624940E+03_r8,     0.2617985E+03_r8,     0.2610975E+03_r8,     0.2603906E+03_r8,&  
         0.2596778E+03_r8,     0.2589586E+03_r8,     0.2582328E+03_r8,     0.2574999E+03_r8,     0.2567598E+03_r8,&  
         0.2560123E+03_r8,     0.2552568E+03_r8,     0.2544932E+03_r8,     0.2537211E+03_r8,     0.2529399E+03_r8,&  
         0.2521495E+03_r8,     0.2513495E+03_r8,     0.2505392E+03_r8,     0.2497184E+03_r8,     0.2488871E+03_r8,&  
         0.2480440E+03_r8,     0.2471894E+03_r8,     0.2463226E+03_r8,     0.2454431E+03_r8,     0.2445507E+03_r8,&  
         0.2436448E+03_r8,     0.2427246E+03_r8,     0.2417906E+03_r8,     0.2408417E+03_r8,     0.2398777E+03_r8,&  
         0.2388984E+03_r8,     0.2379035E+03_r8,     0.2368927E+03_r8,     0.2358657E+03_r8,     0.2348227E+03_r8,&  
         0.2337632E+03_r8,     0.2326877E+03_r8,     0.2315961E+03_r8,     0.2304886E+03_r8,     0.2293653E+03_r8,&  
         0.2282270E+03_r8,     0.2270739E+03_r8,     0.2259067E+03_r8,     0.2247262E+03_r8,     0.2235331E+03_r8,&  
         0.2223281E+03_r8,     0.2211125E+03_r8,     0.2198869E+03_r8,     0.2186525E+03_r8,     0.2174103E+03_r8,&  
         0.2161616E+03_r8,     0.2149072E+03_r8,     0.2136484E+03_r8,     0.2123860E+03_r8,     0.2111213E+03_r8,&  
         0.2098550E+03_r8,     0.2085881E+03_r8,     0.2073215E+03_r8,     0.2060560E+03_r8,     0.2047924E+03_r8/)  

    psaditmk(1:150,  116)= (/ &
         0.3185270E+03_r8,     0.3179120E+03_r8,     0.3172978E+03_r8,     0.3166854E+03_r8,     0.3160738E+03_r8,&  
         0.3154631E+03_r8,     0.3148531E+03_r8,     0.3142444E+03_r8,     0.3136367E+03_r8,     0.3130298E+03_r8,&  
         0.3124246E+03_r8,     0.3118196E+03_r8,     0.3112155E+03_r8,     0.3106120E+03_r8,     0.3100097E+03_r8,&  
         0.3094087E+03_r8,     0.3088081E+03_r8,     0.3082080E+03_r8,     0.3076091E+03_r8,     0.3070104E+03_r8,&  
         0.3064126E+03_r8,     0.3058153E+03_r8,     0.3052188E+03_r8,     0.3046230E+03_r8,     0.3040275E+03_r8,&  
         0.3034329E+03_r8,     0.3028388E+03_r8,     0.3022449E+03_r8,     0.3016516E+03_r8,     0.3010586E+03_r8,&  
         0.3004659E+03_r8,     0.2998740E+03_r8,     0.2992821E+03_r8,     0.2986906E+03_r8,     0.2980991E+03_r8,&  
         0.2975082E+03_r8,     0.2969174E+03_r8,     0.2963266E+03_r8,     0.2957361E+03_r8,     0.2951455E+03_r8,&  
         0.2945546E+03_r8,     0.2939645E+03_r8,     0.2933739E+03_r8,     0.2927831E+03_r8,     0.2921923E+03_r8,&  
         0.2916016E+03_r8,     0.2910105E+03_r8,     0.2904188E+03_r8,     0.2898272E+03_r8,     0.2892352E+03_r8,&  
         0.2886427E+03_r8,     0.2880497E+03_r8,     0.2874563E+03_r8,     0.2868624E+03_r8,     0.2862676E+03_r8,&  
         0.2856723E+03_r8,     0.2850766E+03_r8,     0.2844794E+03_r8,     0.2838816E+03_r8,     0.2832831E+03_r8,&  
         0.2826834E+03_r8,     0.2820827E+03_r8,     0.2814810E+03_r8,     0.2808779E+03_r8,     0.2802736E+03_r8,&  
         0.2796680E+03_r8,     0.2790607E+03_r8,     0.2784523E+03_r8,     0.2778423E+03_r8,     0.2772303E+03_r8,&  
         0.2766169E+03_r8,     0.2760011E+03_r8,     0.2753836E+03_r8,     0.2747639E+03_r8,     0.2741423E+03_r8,&  
         0.2735184E+03_r8,     0.2728917E+03_r8,     0.2722627E+03_r8,     0.2716313E+03_r8,     0.2709965E+03_r8,&  
         0.2703594E+03_r8,     0.2697191E+03_r8,     0.2690757E+03_r8,     0.2684286E+03_r8,     0.2677782E+03_r8,&  
         0.2671244E+03_r8,     0.2664663E+03_r8,     0.2658047E+03_r8,     0.2651385E+03_r8,     0.2644684E+03_r8,&  
         0.2637935E+03_r8,     0.2631134E+03_r8,     0.2624290E+03_r8,     0.2617391E+03_r8,     0.2610437E+03_r8,&  
         0.2603426E+03_r8,     0.2596358E+03_r8,     0.2589225E+03_r8,     0.2582032E+03_r8,     0.2574767E+03_r8,&  
         0.2567432E+03_r8,     0.2560024E+03_r8,     0.2552539E+03_r8,     0.2544973E+03_r8,     0.2537325E+03_r8,&  
         0.2529590E+03_r8,     0.2521763E+03_r8,     0.2513840E+03_r8,     0.2505819E+03_r8,     0.2497697E+03_r8,&  
         0.2489467E+03_r8,     0.2481127E+03_r8,     0.2472671E+03_r8,     0.2464095E+03_r8,     0.2455395E+03_r8,&  
         0.2446567E+03_r8,     0.2437606E+03_r8,     0.2428509E+03_r8,     0.2419271E+03_r8,     0.2409886E+03_r8,&  
         0.2400352E+03_r8,     0.2390669E+03_r8,     0.2380828E+03_r8,     0.2370830E+03_r8,     0.2360672E+03_r8,&  
         0.2350351E+03_r8,     0.2339868E+03_r8,     0.2329223E+03_r8,     0.2318416E+03_r8,     0.2307445E+03_r8,&  
         0.2296320E+03_r8,     0.2285037E+03_r8,     0.2273605E+03_r8,     0.2262029E+03_r8,     0.2250314E+03_r8,&  
         0.2238468E+03_r8,     0.2226500E+03_r8,     0.2214418E+03_r8,     0.2202232E+03_r8,     0.2189952E+03_r8,&  
         0.2177588E+03_r8,     0.2165152E+03_r8,     0.2152655E+03_r8,     0.2140106E+03_r8,     0.2127517E+03_r8,&  
         0.2114898E+03_r8,     0.2102259E+03_r8,     0.2089609E+03_r8,     0.2076958E+03_r8,     0.2064312E+03_r8/)  

    psaditmk(1:150,  117)= (/ &
         0.3187878E+03_r8,     0.3181745E+03_r8,     0.3175619E+03_r8,     0.3169501E+03_r8,     0.3163396E+03_r8,&  
         0.3157304E+03_r8,     0.3151222E+03_r8,     0.3145148E+03_r8,     0.3139087E+03_r8,     0.3133034E+03_r8,&  
         0.3126990E+03_r8,     0.3120958E+03_r8,     0.3114934E+03_r8,     0.3108917E+03_r8,     0.3102911E+03_r8,&  
         0.3096909E+03_r8,     0.3090922E+03_r8,     0.3084939E+03_r8,     0.3078962E+03_r8,     0.3072992E+03_r8,&  
         0.3067038E+03_r8,     0.3061082E+03_r8,     0.3055135E+03_r8,     0.3049196E+03_r8,     0.3043263E+03_r8,&  
         0.3037331E+03_r8,     0.3031406E+03_r8,     0.3025493E+03_r8,     0.3019579E+03_r8,     0.3013667E+03_r8,&  
         0.3007765E+03_r8,     0.3001868E+03_r8,     0.2995969E+03_r8,     0.2990078E+03_r8,     0.2984185E+03_r8,&  
         0.2978299E+03_r8,     0.2972413E+03_r8,     0.2966533E+03_r8,     0.2960650E+03_r8,     0.2954769E+03_r8,&  
         0.2948891E+03_r8,     0.2943015E+03_r8,     0.2937131E+03_r8,     0.2931256E+03_r8,     0.2925376E+03_r8,&  
         0.2919493E+03_r8,     0.2913611E+03_r8,     0.2907729E+03_r8,     0.2901845E+03_r8,     0.2895953E+03_r8,&  
         0.2890060E+03_r8,     0.2884164E+03_r8,     0.2878266E+03_r8,     0.2872356E+03_r8,     0.2866443E+03_r8,&  
         0.2860530E+03_r8,     0.2854608E+03_r8,     0.2848673E+03_r8,     0.2842737E+03_r8,     0.2836788E+03_r8,&  
         0.2830834E+03_r8,     0.2824870E+03_r8,     0.2818893E+03_r8,     0.2812910E+03_r8,     0.2806909E+03_r8,&  
         0.2800901E+03_r8,     0.2794873E+03_r8,     0.2788839E+03_r8,     0.2782789E+03_r8,     0.2776721E+03_r8,&  
         0.2770635E+03_r8,     0.2764534E+03_r8,     0.2758414E+03_r8,     0.2752275E+03_r8,     0.2746118E+03_r8,&  
         0.2739938E+03_r8,     0.2733734E+03_r8,     0.2727509E+03_r8,     0.2721257E+03_r8,     0.2714983E+03_r8,&  
         0.2708682E+03_r8,     0.2702346E+03_r8,     0.2695987E+03_r8,     0.2689593E+03_r8,     0.2683171E+03_r8,&  
         0.2676710E+03_r8,     0.2670214E+03_r8,     0.2663685E+03_r8,     0.2657113E+03_r8,     0.2650502E+03_r8,&  
         0.2643845E+03_r8,     0.2637151E+03_r8,     0.2630404E+03_r8,     0.2623610E+03_r8,     0.2616768E+03_r8,&  
         0.2609868E+03_r8,     0.2602915E+03_r8,     0.2595905E+03_r8,     0.2588833E+03_r8,     0.2581698E+03_r8,&  
         0.2574499E+03_r8,     0.2567228E+03_r8,     0.2559886E+03_r8,     0.2552470E+03_r8,     0.2544974E+03_r8,&  
         0.2537397E+03_r8,     0.2529736E+03_r8,     0.2521986E+03_r8,     0.2514142E+03_r8,     0.2506202E+03_r8,&  
         0.2498161E+03_r8,     0.2490016E+03_r8,     0.2481761E+03_r8,     0.2473394E+03_r8,     0.2464910E+03_r8,&  
         0.2456303E+03_r8,     0.2447572E+03_r8,     0.2438710E+03_r8,     0.2429709E+03_r8,     0.2420573E+03_r8,&  
         0.2411294E+03_r8,     0.2401868E+03_r8,     0.2392290E+03_r8,     0.2382558E+03_r8,     0.2372670E+03_r8,&  
         0.2362621E+03_r8,     0.2352412E+03_r8,     0.2342040E+03_r8,     0.2331505E+03_r8,     0.2320804E+03_r8,&  
         0.2309944E+03_r8,     0.2298923E+03_r8,     0.2287747E+03_r8,     0.2276415E+03_r8,     0.2264935E+03_r8,&  
         0.2253314E+03_r8,     0.2241557E+03_r8,     0.2229672E+03_r8,     0.2217668E+03_r8,     0.2205555E+03_r8,&  
         0.2193342E+03_r8,     0.2181039E+03_r8,     0.2168658E+03_r8,     0.2156209E+03_r8,     0.2143703E+03_r8,&  
         0.2131152E+03_r8,     0.2118564E+03_r8,     0.2105951E+03_r8,     0.2093323E+03_r8,     0.2080688E+03_r8/)  

    psaditmk(1:150,  118)= (/ &
         0.3190461E+03_r8,     0.3184332E+03_r8,     0.3178218E+03_r8,     0.3172119E+03_r8,     0.3166027E+03_r8,&  
         0.3159946E+03_r8,     0.3153873E+03_r8,     0.3147815E+03_r8,     0.3141768E+03_r8,     0.3135729E+03_r8,&  
         0.3129699E+03_r8,     0.3123680E+03_r8,     0.3117667E+03_r8,     0.3111672E+03_r8,     0.3105679E+03_r8,&  
         0.3099693E+03_r8,     0.3093719E+03_r8,     0.3087757E+03_r8,     0.3081794E+03_r8,     0.3075843E+03_r8,&  
         0.3069904E+03_r8,     0.3063967E+03_r8,     0.3058038E+03_r8,     0.3052112E+03_r8,     0.3046196E+03_r8,&  
         0.3040291E+03_r8,     0.3034381E+03_r8,     0.3028485E+03_r8,     0.3022594E+03_r8,     0.3016707E+03_r8,&  
         0.3010822E+03_r8,     0.3004943E+03_r8,     0.2999067E+03_r8,     0.2993194E+03_r8,     0.2987329E+03_r8,&  
         0.2981466E+03_r8,     0.2975602E+03_r8,     0.2969742E+03_r8,     0.2963887E+03_r8,     0.2958031E+03_r8,&  
         0.2952177E+03_r8,     0.2946321E+03_r8,     0.2940472E+03_r8,     0.2934615E+03_r8,     0.2928764E+03_r8,&  
         0.2922915E+03_r8,     0.2917061E+03_r8,     0.2911204E+03_r8,     0.2905345E+03_r8,     0.2899490E+03_r8,&  
         0.2893628E+03_r8,     0.2887762E+03_r8,     0.2881894E+03_r8,     0.2876021E+03_r8,     0.2870145E+03_r8,&  
         0.2864262E+03_r8,     0.2858373E+03_r8,     0.2852480E+03_r8,     0.2846578E+03_r8,     0.2840668E+03_r8,&  
         0.2834753E+03_r8,     0.2828830E+03_r8,     0.2822892E+03_r8,     0.2816949E+03_r8,     0.2810995E+03_r8,&  
         0.2805031E+03_r8,     0.2799051E+03_r8,     0.2793059E+03_r8,     0.2787056E+03_r8,     0.2781036E+03_r8,&  
         0.2775004E+03_r8,     0.2768952E+03_r8,     0.2762887E+03_r8,     0.2756803E+03_r8,     0.2750702E+03_r8,&  
         0.2744578E+03_r8,     0.2738435E+03_r8,     0.2732270E+03_r8,     0.2726080E+03_r8,     0.2719871E+03_r8,&  
         0.2713633E+03_r8,     0.2707368E+03_r8,     0.2701081E+03_r8,     0.2694761E+03_r8,     0.2688411E+03_r8,&  
         0.2682029E+03_r8,     0.2675617E+03_r8,     0.2669164E+03_r8,     0.2662679E+03_r8,     0.2656154E+03_r8,&  
         0.2649591E+03_r8,     0.2642989E+03_r8,     0.2636337E+03_r8,     0.2629645E+03_r8,     0.2622905E+03_r8,&  
         0.2616111E+03_r8,     0.2609269E+03_r8,     0.2602374E+03_r8,     0.2595417E+03_r8,     0.2588406E+03_r8,&  
         0.2581332E+03_r8,     0.2574194E+03_r8,     0.2566988E+03_r8,     0.2559711E+03_r8,     0.2552361E+03_r8,&  
         0.2544935E+03_r8,     0.2537429E+03_r8,     0.2529840E+03_r8,     0.2522163E+03_r8,     0.2514396E+03_r8,&  
         0.2506537E+03_r8,     0.2498577E+03_r8,     0.2490515E+03_r8,     0.2482345E+03_r8,     0.2474066E+03_r8,&  
         0.2465671E+03_r8,     0.2457157E+03_r8,     0.2448519E+03_r8,     0.2439752E+03_r8,     0.2430854E+03_r8,&  
         0.2421819E+03_r8,     0.2412641E+03_r8,     0.2403319E+03_r8,     0.2393847E+03_r8,     0.2384224E+03_r8,&  
         0.2374445E+03_r8,     0.2364507E+03_r8,     0.2354408E+03_r8,     0.2344148E+03_r8,     0.2333721E+03_r8,&  
         0.2323133E+03_r8,     0.2312382E+03_r8,     0.2301468E+03_r8,     0.2290395E+03_r8,     0.2279167E+03_r8,&  
         0.2267786E+03_r8,     0.2256260E+03_r8,     0.2244594E+03_r8,     0.2232795E+03_r8,     0.2220873E+03_r8,&  
         0.2208835E+03_r8,     0.2196692E+03_r8,     0.2184454E+03_r8,     0.2172131E+03_r8,     0.2159734E+03_r8,&  
         0.2147275E+03_r8,     0.2134763E+03_r8,     0.2122210E+03_r8,     0.2109627E+03_r8,     0.2097022E+03_r8/)  

    psaditmk(1:150,  119)= (/ &
         0.3193003E+03_r8,     0.3186888E+03_r8,     0.3180786E+03_r8,     0.3174695E+03_r8,     0.3168615E+03_r8,&  
         0.3162549E+03_r8,     0.3156494E+03_r8,     0.3150446E+03_r8,     0.3144410E+03_r8,     0.3138389E+03_r8,&  
         0.3132371E+03_r8,     0.3126365E+03_r8,     0.3120369E+03_r8,     0.3114385E+03_r8,     0.3108409E+03_r8,&  
         0.3102442E+03_r8,     0.3096482E+03_r8,     0.3090529E+03_r8,     0.3084586E+03_r8,     0.3078654E+03_r8,&  
         0.3072729E+03_r8,     0.3066806E+03_r8,     0.3060892E+03_r8,     0.3054993E+03_r8,     0.3049089E+03_r8,&  
         0.3043201E+03_r8,     0.3037317E+03_r8,     0.3031434E+03_r8,     0.3025561E+03_r8,     0.3019691E+03_r8,&  
         0.3013830E+03_r8,     0.3007970E+03_r8,     0.3002115E+03_r8,     0.2996267E+03_r8,     0.2990420E+03_r8,&  
         0.2984576E+03_r8,     0.2978739E+03_r8,     0.2972900E+03_r8,     0.2967067E+03_r8,     0.2961234E+03_r8,&  
         0.2955406E+03_r8,     0.2949577E+03_r8,     0.2943749E+03_r8,     0.2937923E+03_r8,     0.2932098E+03_r8,&  
         0.2926268E+03_r8,     0.2920445E+03_r8,     0.2914619E+03_r8,     0.2908790E+03_r8,     0.2902960E+03_r8,&  
         0.2897128E+03_r8,     0.2891297E+03_r8,     0.2885459E+03_r8,     0.2879616E+03_r8,     0.2873770E+03_r8,&  
         0.2867922E+03_r8,     0.2862070E+03_r8,     0.2856207E+03_r8,     0.2850342E+03_r8,     0.2844471E+03_r8,&  
         0.2838593E+03_r8,     0.2832707E+03_r8,     0.2826810E+03_r8,     0.2820907E+03_r8,     0.2814993E+03_r8,&  
         0.2809070E+03_r8,     0.2803135E+03_r8,     0.2797190E+03_r8,     0.2791231E+03_r8,     0.2785256E+03_r8,&  
         0.2779274E+03_r8,     0.2773272E+03_r8,     0.2767259E+03_r8,     0.2761227E+03_r8,     0.2755180E+03_r8,&  
         0.2749110E+03_r8,     0.2743023E+03_r8,     0.2736918E+03_r8,     0.2730788E+03_r8,     0.2724638E+03_r8,&  
         0.2718463E+03_r8,     0.2712267E+03_r8,     0.2706043E+03_r8,     0.2699794E+03_r8,     0.2693515E+03_r8,&  
         0.2687207E+03_r8,     0.2680869E+03_r8,     0.2674496E+03_r8,     0.2668094E+03_r8,     0.2661649E+03_r8,&  
         0.2655174E+03_r8,     0.2648654E+03_r8,     0.2642101E+03_r8,     0.2635503E+03_r8,     0.2628857E+03_r8,&  
         0.2622170E+03_r8,     0.2615430E+03_r8,     0.2608640E+03_r8,     0.2601798E+03_r8,     0.2594901E+03_r8,&  
         0.2587947E+03_r8,     0.2580930E+03_r8,     0.2573855E+03_r8,     0.2566712E+03_r8,     0.2559499E+03_r8,&  
         0.2552215E+03_r8,     0.2544856E+03_r8,     0.2537420E+03_r8,     0.2529902E+03_r8,     0.2522300E+03_r8,&  
         0.2514608E+03_r8,     0.2506823E+03_r8,     0.2498944E+03_r8,     0.2490963E+03_r8,     0.2482880E+03_r8,&  
         0.2474687E+03_r8,     0.2466380E+03_r8,     0.2457957E+03_r8,     0.2449414E+03_r8,     0.2440741E+03_r8,&  
         0.2431940E+03_r8,     0.2423003E+03_r8,     0.2413928E+03_r8,     0.2404710E+03_r8,     0.2395344E+03_r8,&  
         0.2385828E+03_r8,     0.2376157E+03_r8,     0.2366329E+03_r8,     0.2356340E+03_r8,     0.2346189E+03_r8,&  
         0.2335875E+03_r8,     0.2325398E+03_r8,     0.2314757E+03_r8,     0.2303951E+03_r8,     0.2292983E+03_r8,&  
         0.2281859E+03_r8,     0.2270580E+03_r8,     0.2259151E+03_r8,     0.2247578E+03_r8,     0.2235869E+03_r8,&  
         0.2224030E+03_r8,     0.2212072E+03_r8,     0.2200002E+03_r8,     0.2187831E+03_r8,     0.2175569E+03_r8,&  
         0.2163228E+03_r8,     0.2150818E+03_r8,     0.2138350E+03_r8,     0.2125834E+03_r8,     0.2113282E+03_r8/)  

    psaditmk(1:150,  120)= (/ &
         0.3195515E+03_r8,     0.3189412E+03_r8,     0.3183321E+03_r8,     0.3177245E+03_r8,     0.3171178E+03_r8,&  
         0.3165120E+03_r8,     0.3159077E+03_r8,     0.3153046E+03_r8,     0.3147021E+03_r8,     0.3141010E+03_r8,&  
         0.3135009E+03_r8,     0.3129019E+03_r8,     0.3123036E+03_r8,     0.3117065E+03_r8,     0.3111101E+03_r8,&  
         0.3105146E+03_r8,     0.3099203E+03_r8,     0.3093270E+03_r8,     0.3087339E+03_r8,     0.3081422E+03_r8,&  
         0.3075513E+03_r8,     0.3069609E+03_r8,     0.3063715E+03_r8,     0.3057821E+03_r8,     0.3051945E+03_r8,&  
         0.3046072E+03_r8,     0.3040201E+03_r8,     0.3034340E+03_r8,     0.3028485E+03_r8,     0.3022633E+03_r8,&  
         0.3016790E+03_r8,     0.3010952E+03_r8,     0.3005120E+03_r8,     0.2999285E+03_r8,     0.2993460E+03_r8,&  
         0.2987641E+03_r8,     0.2981820E+03_r8,     0.2976008E+03_r8,     0.2970198E+03_r8,     0.2964388E+03_r8,&  
         0.2958580E+03_r8,     0.2952775E+03_r8,     0.2946973E+03_r8,     0.2941172E+03_r8,     0.2935371E+03_r8,&  
         0.2929572E+03_r8,     0.2923773E+03_r8,     0.2917969E+03_r8,     0.2912169E+03_r8,     0.2906372E+03_r8,&  
         0.2900567E+03_r8,     0.2894759E+03_r8,     0.2888953E+03_r8,     0.2883146E+03_r8,     0.2877332E+03_r8,&  
         0.2871513E+03_r8,     0.2865693E+03_r8,     0.2859866E+03_r8,     0.2854036E+03_r8,     0.2848200E+03_r8,&  
         0.2842356E+03_r8,     0.2836508E+03_r8,     0.2830651E+03_r8,     0.2824783E+03_r8,     0.2818912E+03_r8,&  
         0.2813026E+03_r8,     0.2807136E+03_r8,     0.2801232E+03_r8,     0.2795316E+03_r8,     0.2789389E+03_r8,&  
         0.2783452E+03_r8,     0.2777497E+03_r8,     0.2771530E+03_r8,     0.2765549E+03_r8,     0.2759553E+03_r8,&  
         0.2753538E+03_r8,     0.2747505E+03_r8,     0.2741454E+03_r8,     0.2735381E+03_r8,     0.2729290E+03_r8,&  
         0.2723177E+03_r8,     0.2717043E+03_r8,     0.2710880E+03_r8,     0.2704698E+03_r8,     0.2698488E+03_r8,&  
         0.2692249E+03_r8,     0.2685984E+03_r8,     0.2679685E+03_r8,     0.2673358E+03_r8,     0.2666993E+03_r8,&  
         0.2660602E+03_r8,     0.2654167E+03_r8,     0.2647698E+03_r8,     0.2641187E+03_r8,     0.2634639E+03_r8,&  
         0.2628044E+03_r8,     0.2621402E+03_r8,     0.2614719E+03_r8,     0.2607981E+03_r8,     0.2601193E+03_r8,&  
         0.2594354E+03_r8,     0.2587456E+03_r8,     0.2580499E+03_r8,     0.2573482E+03_r8,     0.2566400E+03_r8,&  
         0.2559251E+03_r8,     0.2552032E+03_r8,     0.2544740E+03_r8,     0.2537372E+03_r8,     0.2529924E+03_r8,&  
         0.2522393E+03_r8,     0.2514776E+03_r8,     0.2507068E+03_r8,     0.2499266E+03_r8,     0.2491367E+03_r8,&  
         0.2483364E+03_r8,     0.2475256E+03_r8,     0.2467037E+03_r8,     0.2458704E+03_r8,     0.2450251E+03_r8,&  
         0.2441672E+03_r8,     0.2432969E+03_r8,     0.2424132E+03_r8,     0.2415157E+03_r8,     0.2406040E+03_r8,&  
         0.2396779E+03_r8,     0.2387369E+03_r8,     0.2377806E+03_r8,     0.2368087E+03_r8,     0.2358208E+03_r8,&  
         0.2348168E+03_r8,     0.2337966E+03_r8,     0.2327599E+03_r8,     0.2317066E+03_r8,     0.2306370E+03_r8,&  
         0.2295511E+03_r8,     0.2284493E+03_r8,     0.2273315E+03_r8,     0.2261987E+03_r8,     0.2250509E+03_r8,&  
         0.2238892E+03_r8,     0.2227140E+03_r8,     0.2215263E+03_r8,     0.2203269E+03_r8,     0.2191169E+03_r8,&  
         0.2178972E+03_r8,     0.2166689E+03_r8,     0.2154331E+03_r8,     0.2141910E+03_r8,     0.2129435E+03_r8/)  

    psaditmk(1:150,  121)= (/ &
         0.3197993E+03_r8,     0.3191903E+03_r8,     0.3185823E+03_r8,     0.3179757E+03_r8,     0.3173701E+03_r8,&  
         0.3167660E+03_r8,     0.3161627E+03_r8,     0.3155605E+03_r8,     0.3149597E+03_r8,     0.3143598E+03_r8,&  
         0.3137607E+03_r8,     0.3131629E+03_r8,     0.3125663E+03_r8,     0.3119706E+03_r8,     0.3113759E+03_r8,&  
         0.3107817E+03_r8,     0.3101886E+03_r8,     0.3095967E+03_r8,     0.3090057E+03_r8,     0.3084149E+03_r8,&  
         0.3078254E+03_r8,     0.3072368E+03_r8,     0.3066491E+03_r8,     0.3060614E+03_r8,     0.3054754E+03_r8,&  
         0.3048897E+03_r8,     0.3043047E+03_r8,     0.3037200E+03_r8,     0.3031365E+03_r8,     0.3025534E+03_r8,&  
         0.3019707E+03_r8,     0.3013887E+03_r8,     0.3008071E+03_r8,     0.3002264E+03_r8,     0.2996454E+03_r8,&  
         0.2990653E+03_r8,     0.2984859E+03_r8,     0.2979066E+03_r8,     0.2973275E+03_r8,     0.2967488E+03_r8,&  
         0.2961705E+03_r8,     0.2955923E+03_r8,     0.2950144E+03_r8,     0.2944368E+03_r8,     0.2938591E+03_r8,&  
         0.2932815E+03_r8,     0.2927037E+03_r8,     0.2921268E+03_r8,     0.2915492E+03_r8,     0.2909717E+03_r8,&  
         0.2903942E+03_r8,     0.2898168E+03_r8,     0.2892388E+03_r8,     0.2886607E+03_r8,     0.2880823E+03_r8,&  
         0.2875038E+03_r8,     0.2869248E+03_r8,     0.2863456E+03_r8,     0.2857657E+03_r8,     0.2851854E+03_r8,&  
         0.2846047E+03_r8,     0.2840230E+03_r8,     0.2834411E+03_r8,     0.2828582E+03_r8,     0.2822747E+03_r8,&  
         0.2816902E+03_r8,     0.2811049E+03_r8,     0.2805188E+03_r8,     0.2799316E+03_r8,     0.2793430E+03_r8,&  
         0.2787540E+03_r8,     0.2781631E+03_r8,     0.2775710E+03_r8,     0.2769777E+03_r8,     0.2763827E+03_r8,&  
         0.2757864E+03_r8,     0.2751881E+03_r8,     0.2745881E+03_r8,     0.2739866E+03_r8,     0.2733832E+03_r8,&  
         0.2727775E+03_r8,     0.2721700E+03_r8,     0.2715599E+03_r8,     0.2709481E+03_r8,     0.2703331E+03_r8,&  
         0.2697164E+03_r8,     0.2690963E+03_r8,     0.2684739E+03_r8,     0.2678482E+03_r8,     0.2672196E+03_r8,&  
         0.2665881E+03_r8,     0.2659526E+03_r8,     0.2653140E+03_r8,     0.2646714E+03_r8,     0.2640254E+03_r8,&  
         0.2633746E+03_r8,     0.2627204E+03_r8,     0.2620616E+03_r8,     0.2613978E+03_r8,     0.2607295E+03_r8,&  
         0.2600560E+03_r8,     0.2593774E+03_r8,     0.2586933E+03_r8,     0.2580035E+03_r8,     0.2573074E+03_r8,&  
         0.2566055E+03_r8,     0.2558968E+03_r8,     0.2551812E+03_r8,     0.2544585E+03_r8,     0.2537284E+03_r8,&  
         0.2529906E+03_r8,     0.2522447E+03_r8,     0.2514900E+03_r8,     0.2507266E+03_r8,     0.2499542E+03_r8,&  
         0.2491722E+03_r8,     0.2483802E+03_r8,     0.2475779E+03_r8,     0.2467644E+03_r8,     0.2459398E+03_r8,&  
         0.2451035E+03_r8,     0.2442551E+03_r8,     0.2433941E+03_r8,     0.2425199E+03_r8,     0.2416325E+03_r8,&  
         0.2407309E+03_r8,     0.2398154E+03_r8,     0.2388848E+03_r8,     0.2379392E+03_r8,     0.2369781E+03_r8,&  
         0.2360012E+03_r8,     0.2350083E+03_r8,     0.2339991E+03_r8,     0.2329734E+03_r8,     0.2319312E+03_r8,&  
         0.2308726E+03_r8,     0.2297977E+03_r8,     0.2287064E+03_r8,     0.2275993E+03_r8,     0.2264764E+03_r8,&  
         0.2253385E+03_r8,     0.2241861E+03_r8,     0.2230199E+03_r8,     0.2218407E+03_r8,     0.2206493E+03_r8,&  
         0.2194466E+03_r8,     0.2182337E+03_r8,     0.2170116E+03_r8,     0.2157814E+03_r8,     0.2145442E+03_r8/)  

    psaditmk(1:150,  122)= (/ &
         0.3200440E+03_r8,     0.3194363E+03_r8,     0.3188294E+03_r8,     0.3182239E+03_r8,     0.3176196E+03_r8,&  
         0.3170165E+03_r8,     0.3164142E+03_r8,     0.3158132E+03_r8,     0.3152135E+03_r8,     0.3146149E+03_r8,&  
         0.3140175E+03_r8,     0.3134208E+03_r8,     0.3128253E+03_r8,     0.3122308E+03_r8,     0.3116376E+03_r8,&  
         0.3110449E+03_r8,     0.3104532E+03_r8,     0.3098625E+03_r8,     0.3092728E+03_r8,     0.3086845E+03_r8,&  
         0.3080963E+03_r8,     0.3075090E+03_r8,     0.3069227E+03_r8,     0.3063373E+03_r8,     0.3057524E+03_r8,&  
         0.3051679E+03_r8,     0.3045846E+03_r8,     0.3040023E+03_r8,     0.3034203E+03_r8,     0.3028387E+03_r8,&  
         0.3022580E+03_r8,     0.3016780E+03_r8,     0.3010983E+03_r8,     0.3005187E+03_r8,     0.2999407E+03_r8,&  
         0.2993625E+03_r8,     0.2987845E+03_r8,     0.2982075E+03_r8,     0.2976305E+03_r8,     0.2970539E+03_r8,&  
         0.2964777E+03_r8,     0.2959018E+03_r8,     0.2953263E+03_r8,     0.2947506E+03_r8,     0.2941754E+03_r8,&  
         0.2936006E+03_r8,     0.2930254E+03_r8,     0.2924505E+03_r8,     0.2918756E+03_r8,     0.2913011E+03_r8,&  
         0.2907262E+03_r8,     0.2901510E+03_r8,     0.2895759E+03_r8,     0.2890007E+03_r8,     0.2884254E+03_r8,&  
         0.2878501E+03_r8,     0.2872738E+03_r8,     0.2866976E+03_r8,     0.2861211E+03_r8,     0.2855441E+03_r8,&  
         0.2849664E+03_r8,     0.2843884E+03_r8,     0.2838102E+03_r8,     0.2832304E+03_r8,     0.2826509E+03_r8,&  
         0.2820701E+03_r8,     0.2814885E+03_r8,     0.2809063E+03_r8,     0.2803232E+03_r8,     0.2797390E+03_r8,&  
         0.2791537E+03_r8,     0.2785674E+03_r8,     0.2779795E+03_r8,     0.2773909E+03_r8,     0.2768006E+03_r8,&  
         0.2762089E+03_r8,     0.2756158E+03_r8,     0.2750211E+03_r8,     0.2744247E+03_r8,     0.2738265E+03_r8,&  
         0.2732265E+03_r8,     0.2726247E+03_r8,     0.2720204E+03_r8,     0.2714145E+03_r8,     0.2708061E+03_r8,&  
         0.2701954E+03_r8,     0.2695820E+03_r8,     0.2689663E+03_r8,     0.2683476E+03_r8,     0.2677261E+03_r8,&  
         0.2671015E+03_r8,     0.2664740E+03_r8,     0.2658430E+03_r8,     0.2652088E+03_r8,     0.2645706E+03_r8,&  
         0.2639292E+03_r8,     0.2632837E+03_r8,     0.2626336E+03_r8,     0.2619798E+03_r8,     0.2613212E+03_r8,&  
         0.2606580E+03_r8,     0.2599901E+03_r8,     0.2593166E+03_r8,     0.2586378E+03_r8,     0.2579539E+03_r8,&  
         0.2572637E+03_r8,     0.2565674E+03_r8,     0.2558649E+03_r8,     0.2551556E+03_r8,     0.2544393E+03_r8,&  
         0.2537159E+03_r8,     0.2529851E+03_r8,     0.2522458E+03_r8,     0.2514986E+03_r8,     0.2507427E+03_r8,&  
         0.2499778E+03_r8,     0.2492034E+03_r8,     0.2484193E+03_r8,     0.2476251E+03_r8,     0.2468201E+03_r8,&  
         0.2460042E+03_r8,     0.2451767E+03_r8,     0.2443376E+03_r8,     0.2434859E+03_r8,     0.2426215E+03_r8,&  
         0.2417438E+03_r8,     0.2408523E+03_r8,     0.2399467E+03_r8,     0.2390268E+03_r8,     0.2380915E+03_r8,&  
         0.2371413E+03_r8,     0.2361753E+03_r8,     0.2351934E+03_r8,     0.2341951E+03_r8,     0.2331806E+03_r8,&  
         0.2321496E+03_r8,     0.2311021E+03_r8,     0.2300381E+03_r8,     0.2289575E+03_r8,     0.2278610E+03_r8,&  
         0.2267485E+03_r8,     0.2256205E+03_r8,     0.2244778E+03_r8,     0.2233208E+03_r8,     0.2221503E+03_r8,&  
         0.2209671E+03_r8,     0.2197720E+03_r8,     0.2185663E+03_r8,     0.2173507E+03_r8,     0.2161263E+03_r8/)  

    psaditmk(1:150,  123)= (/ &
         0.3202856E+03_r8,     0.3196789E+03_r8,     0.3190735E+03_r8,     0.3184689E+03_r8,     0.3178657E+03_r8,&  
         0.3172638E+03_r8,     0.3166631E+03_r8,     0.3160629E+03_r8,     0.3154644E+03_r8,     0.3148669E+03_r8,&  
         0.3142706E+03_r8,     0.3136752E+03_r8,     0.3130813E+03_r8,     0.3124879E+03_r8,     0.3118958E+03_r8,&  
         0.3113049E+03_r8,     0.3107148E+03_r8,     0.3101251E+03_r8,     0.3095367E+03_r8,     0.3089495E+03_r8,&  
         0.3083633E+03_r8,     0.3077775E+03_r8,     0.3071926E+03_r8,     0.3066083E+03_r8,     0.3060253E+03_r8,&  
         0.3054427E+03_r8,     0.3048611E+03_r8,     0.3042798E+03_r8,     0.3037000E+03_r8,     0.3031199E+03_r8,&  
         0.3025408E+03_r8,     0.3019626E+03_r8,     0.3013848E+03_r8,     0.3008074E+03_r8,     0.3002305E+03_r8,&  
         0.2996545E+03_r8,     0.2990788E+03_r8,     0.2985036E+03_r8,     0.2979285E+03_r8,     0.2973542E+03_r8,&  
         0.2967804E+03_r8,     0.2962062E+03_r8,     0.2956327E+03_r8,     0.2950598E+03_r8,     0.2944867E+03_r8,&  
         0.2939139E+03_r8,     0.2933413E+03_r8,     0.2927689E+03_r8,     0.2921967E+03_r8,     0.2916241E+03_r8,&  
         0.2910518E+03_r8,     0.2904798E+03_r8,     0.2899074E+03_r8,     0.2893348E+03_r8,     0.2887622E+03_r8,&  
         0.2881893E+03_r8,     0.2876165E+03_r8,     0.2870435E+03_r8,     0.2864697E+03_r8,     0.2858957E+03_r8,&  
         0.2853216E+03_r8,     0.2847472E+03_r8,     0.2841714E+03_r8,     0.2835958E+03_r8,     0.2830193E+03_r8,&  
         0.2824424E+03_r8,     0.2818643E+03_r8,     0.2812861E+03_r8,     0.2807065E+03_r8,     0.2801266E+03_r8,&  
         0.2795453E+03_r8,     0.2789631E+03_r8,     0.2783796E+03_r8,     0.2777952E+03_r8,     0.2772093E+03_r8,&  
         0.2766226E+03_r8,     0.2760340E+03_r8,     0.2754442E+03_r8,     0.2748528E+03_r8,     0.2742599E+03_r8,&  
         0.2736650E+03_r8,     0.2730686E+03_r8,     0.2724700E+03_r8,     0.2718698E+03_r8,     0.2712672E+03_r8,&  
         0.2706625E+03_r8,     0.2700555E+03_r8,     0.2694461E+03_r8,     0.2688339E+03_r8,     0.2682195E+03_r8,&  
         0.2676018E+03_r8,     0.2669815E+03_r8,     0.2663579E+03_r8,     0.2657315E+03_r8,     0.2651014E+03_r8,&  
         0.2644680E+03_r8,     0.2638307E+03_r8,     0.2631897E+03_r8,     0.2625448E+03_r8,     0.2618955E+03_r8,&  
         0.2612421E+03_r8,     0.2605835E+03_r8,     0.2599207E+03_r8,     0.2592529E+03_r8,     0.2585796E+03_r8,&  
         0.2579010E+03_r8,     0.2572166E+03_r8,     0.2565262E+03_r8,     0.2558297E+03_r8,     0.2551267E+03_r8,&  
         0.2544165E+03_r8,     0.2536996E+03_r8,     0.2529754E+03_r8,     0.2522431E+03_r8,     0.2515029E+03_r8,&  
         0.2507542E+03_r8,     0.2499968E+03_r8,     0.2492301E+03_r8,     0.2484538E+03_r8,     0.2476676E+03_r8,&  
         0.2468710E+03_r8,     0.2460636E+03_r8,     0.2452450E+03_r8,     0.2444146E+03_r8,     0.2435722E+03_r8,&  
         0.2427170E+03_r8,     0.2418490E+03_r8,     0.2409676E+03_r8,     0.2400722E+03_r8,     0.2391624E+03_r8,&  
         0.2382379E+03_r8,     0.2372983E+03_r8,     0.2363430E+03_r8,     0.2353720E+03_r8,     0.2343849E+03_r8,&  
         0.2333816E+03_r8,     0.2323617E+03_r8,     0.2313251E+03_r8,     0.2302720E+03_r8,     0.2292025E+03_r8,&  
         0.2281164E+03_r8,     0.2270145E+03_r8,     0.2258969E+03_r8,     0.2247639E+03_r8,     0.2236163E+03_r8,&  
         0.2224548E+03_r8,     0.2212801E+03_r8,     0.2200930E+03_r8,     0.2188947E+03_r8,     0.2176859E+03_r8/)  

    psaditmk(1:150,  124)= (/ &
         0.3205251E+03_r8,     0.3199188E+03_r8,     0.3193144E+03_r8,     0.3187111E+03_r8,     0.3181088E+03_r8,&  
         0.3175078E+03_r8,     0.3169078E+03_r8,     0.3163096E+03_r8,     0.3157119E+03_r8,     0.3151158E+03_r8,&  
         0.3145205E+03_r8,     0.3139266E+03_r8,     0.3133334E+03_r8,     0.3127418E+03_r8,     0.3121510E+03_r8,&  
         0.3115609E+03_r8,     0.3109721E+03_r8,     0.3103846E+03_r8,     0.3097975E+03_r8,     0.3092114E+03_r8,&  
         0.3086263E+03_r8,     0.3080421E+03_r8,     0.3074590E+03_r8,     0.3068762E+03_r8,     0.3062945E+03_r8,&  
         0.3057133E+03_r8,     0.3051334E+03_r8,     0.3045541E+03_r8,     0.3039749E+03_r8,     0.3033968E+03_r8,&  
         0.3028199E+03_r8,     0.3022431E+03_r8,     0.3016670E+03_r8,     0.3010917E+03_r8,     0.3005168E+03_r8,&  
         0.2999421E+03_r8,     0.2993681E+03_r8,     0.2987951E+03_r8,     0.2982223E+03_r8,     0.2976497E+03_r8,&  
         0.2970778E+03_r8,     0.2965062E+03_r8,     0.2959347E+03_r8,     0.2953636E+03_r8,     0.2947929E+03_r8,&  
         0.2942224E+03_r8,     0.2936522E+03_r8,     0.2930818E+03_r8,     0.2925119E+03_r8,     0.2919421E+03_r8,&  
         0.2913725E+03_r8,     0.2908026E+03_r8,     0.2902327E+03_r8,     0.2896631E+03_r8,     0.2890932E+03_r8,&  
         0.2885231E+03_r8,     0.2879529E+03_r8,     0.2873826E+03_r8,     0.2868122E+03_r8,     0.2862412E+03_r8,&  
         0.2856698E+03_r8,     0.2850983E+03_r8,     0.2845263E+03_r8,     0.2839539E+03_r8,     0.2833806E+03_r8,&  
         0.2828076E+03_r8,     0.2822328E+03_r8,     0.2816582E+03_r8,     0.2810824E+03_r8,     0.2805060E+03_r8,&  
         0.2799288E+03_r8,     0.2793505E+03_r8,     0.2787712E+03_r8,     0.2781912E+03_r8,     0.2776096E+03_r8,&  
         0.2770269E+03_r8,     0.2764432E+03_r8,     0.2758579E+03_r8,     0.2752713E+03_r8,     0.2746832E+03_r8,&  
         0.2740935E+03_r8,     0.2735023E+03_r8,     0.2729090E+03_r8,     0.2723142E+03_r8,     0.2717171E+03_r8,&  
         0.2711185E+03_r8,     0.2705172E+03_r8,     0.2699140E+03_r8,     0.2693081E+03_r8,     0.2687002E+03_r8,&  
         0.2680891E+03_r8,     0.2674758E+03_r8,     0.2668592E+03_r8,     0.2662402E+03_r8,     0.2656176E+03_r8,&  
         0.2649919E+03_r8,     0.2643628E+03_r8,     0.2637302E+03_r8,     0.2630937E+03_r8,     0.2624532E+03_r8,&  
         0.2618090E+03_r8,     0.2611600E+03_r8,     0.2605069E+03_r8,     0.2598488E+03_r8,     0.2591863E+03_r8,&  
         0.2585183E+03_r8,     0.2578452E+03_r8,     0.2571663E+03_r8,     0.2564818E+03_r8,     0.2557911E+03_r8,&  
         0.2550941E+03_r8,     0.2543903E+03_r8,     0.2536796E+03_r8,     0.2529620E+03_r8,     0.2522364E+03_r8,&  
         0.2515032E+03_r8,     0.2507617E+03_r8,     0.2500112E+03_r8,     0.2492521E+03_r8,     0.2484836E+03_r8,&  
         0.2477053E+03_r8,     0.2469169E+03_r8,     0.2461179E+03_r8,     0.2453079E+03_r8,     0.2444865E+03_r8,&  
         0.2436531E+03_r8,     0.2428076E+03_r8,     0.2419491E+03_r8,     0.2410772E+03_r8,     0.2401917E+03_r8,&  
         0.2392922E+03_r8,     0.2383781E+03_r8,     0.2374490E+03_r8,     0.2365045E+03_r8,     0.2355444E+03_r8,&  
         0.2345683E+03_r8,     0.2335761E+03_r8,     0.2325671E+03_r8,     0.2315417E+03_r8,     0.2304997E+03_r8,&  
         0.2294411E+03_r8,     0.2283660E+03_r8,     0.2272745E+03_r8,     0.2261672E+03_r8,     0.2250444E+03_r8,&  
         0.2239065E+03_r8,     0.2227542E+03_r8,     0.2215883E+03_r8,     0.2204095E+03_r8,     0.2192188E+03_r8/)  

    psaditmk(1:150,  125)= (/ &
         0.3207607E+03_r8,     0.3201557E+03_r8,     0.3195525E+03_r8,     0.3189500E+03_r8,     0.3183489E+03_r8,&  
         0.3177490E+03_r8,     0.3171505E+03_r8,     0.3165531E+03_r8,     0.3159566E+03_r8,     0.3153612E+03_r8,&  
         0.3147674E+03_r8,     0.3141744E+03_r8,     0.3135829E+03_r8,     0.3129921E+03_r8,     0.3124026E+03_r8,&  
         0.3118138E+03_r8,     0.3112263E+03_r8,     0.3106400E+03_r8,     0.3100540E+03_r8,     0.3094699E+03_r8,&  
         0.3088856E+03_r8,     0.3083029E+03_r8,     0.3077212E+03_r8,     0.3071402E+03_r8,     0.3065598E+03_r8,&  
         0.3059802E+03_r8,     0.3054013E+03_r8,     0.3048240E+03_r8,     0.3042468E+03_r8,     0.3036700E+03_r8,&  
         0.3030945E+03_r8,     0.3025199E+03_r8,     0.3019453E+03_r8,     0.3013715E+03_r8,     0.3007983E+03_r8,&  
         0.3002261E+03_r8,     0.2996541E+03_r8,     0.2990821E+03_r8,     0.2985112E+03_r8,     0.2979411E+03_r8,&  
         0.2973706E+03_r8,     0.2968008E+03_r8,     0.2962316E+03_r8,     0.2956628E+03_r8,     0.2950941E+03_r8,&  
         0.2945258E+03_r8,     0.2939576E+03_r8,     0.2933899E+03_r8,     0.2928222E+03_r8,     0.2922545E+03_r8,&  
         0.2916872E+03_r8,     0.2911198E+03_r8,     0.2905528E+03_r8,     0.2899855E+03_r8,     0.2894182E+03_r8,&  
         0.2888508E+03_r8,     0.2882836E+03_r8,     0.2877159E+03_r8,     0.2871480E+03_r8,     0.2865800E+03_r8,&  
         0.2860121E+03_r8,     0.2854436E+03_r8,     0.2848746E+03_r8,     0.2843050E+03_r8,     0.2837356E+03_r8,&  
         0.2831650E+03_r8,     0.2825944E+03_r8,     0.2820229E+03_r8,     0.2814507E+03_r8,     0.2808777E+03_r8,&  
         0.2803046E+03_r8,     0.2797299E+03_r8,     0.2791548E+03_r8,     0.2785786E+03_r8,     0.2780013E+03_r8,&  
         0.2774228E+03_r8,     0.2768435E+03_r8,     0.2762629E+03_r8,     0.2756806E+03_r8,     0.2750974E+03_r8,&  
         0.2745125E+03_r8,     0.2739262E+03_r8,     0.2733380E+03_r8,     0.2727484E+03_r8,     0.2721567E+03_r8,&  
         0.2715634E+03_r8,     0.2709677E+03_r8,     0.2703705E+03_r8,     0.2697706E+03_r8,     0.2691689E+03_r8,&  
         0.2685641E+03_r8,     0.2679576E+03_r8,     0.2673476E+03_r8,     0.2667354E+03_r8,     0.2661201E+03_r8,&  
         0.2655018E+03_r8,     0.2648803E+03_r8,     0.2642556E+03_r8,     0.2636274E+03_r8,     0.2629951E+03_r8,&  
         0.2623595E+03_r8,     0.2617195E+03_r8,     0.2610757E+03_r8,     0.2604273E+03_r8,     0.2597747E+03_r8,&  
         0.2591169E+03_r8,     0.2584542E+03_r8,     0.2577863E+03_r8,     0.2571131E+03_r8,     0.2564341E+03_r8,&  
         0.2557494E+03_r8,     0.2550582E+03_r8,     0.2543605E+03_r8,     0.2536561E+03_r8,     0.2529447E+03_r8,&  
         0.2522260E+03_r8,     0.2514995E+03_r8,     0.2507650E+03_r8,     0.2500217E+03_r8,     0.2492699E+03_r8,&  
         0.2485090E+03_r8,     0.2477389E+03_r8,     0.2469582E+03_r8,     0.2461673E+03_r8,     0.2453660E+03_r8,&  
         0.2445530E+03_r8,     0.2437287E+03_r8,     0.2428922E+03_r8,     0.2420432E+03_r8,     0.2411811E+03_r8,&  
         0.2403057E+03_r8,     0.2394160E+03_r8,     0.2385123E+03_r8,     0.2375936E+03_r8,     0.2366599E+03_r8,&  
         0.2357106E+03_r8,     0.2347455E+03_r8,     0.2337640E+03_r8,     0.2327663E+03_r8,     0.2317520E+03_r8,&  
         0.2307210E+03_r8,     0.2296735E+03_r8,     0.2286093E+03_r8,     0.2275287E+03_r8,     0.2264317E+03_r8,&  
         0.2253191E+03_r8,     0.2241911E+03_r8,     0.2230483E+03_r8,     0.2218915E+03_r8,     0.2207213E+03_r8/)  

    psaditmk(1:150,  126)= (/ &
         0.3209937E+03_r8,     0.3203900E+03_r8,     0.3197874E+03_r8,     0.3191863E+03_r8,     0.3185862E+03_r8,&  
         0.3179874E+03_r8,     0.3173896E+03_r8,     0.3167934E+03_r8,     0.3161978E+03_r8,     0.3156042E+03_r8,&  
         0.3150111E+03_r8,     0.3144196E+03_r8,     0.3138293E+03_r8,     0.3132395E+03_r8,     0.3126512E+03_r8,&  
         0.3120635E+03_r8,     0.3114771E+03_r8,     0.3108918E+03_r8,     0.3103073E+03_r8,     0.3097245E+03_r8,&  
         0.3091418E+03_r8,     0.3085606E+03_r8,     0.3079797E+03_r8,     0.3074003E+03_r8,     0.3068215E+03_r8,&  
         0.3062435E+03_r8,     0.3056663E+03_r8,     0.3050896E+03_r8,     0.3045139E+03_r8,     0.3039396E+03_r8,&  
         0.3033653E+03_r8,     0.3027919E+03_r8,     0.3022190E+03_r8,     0.3016474E+03_r8,     0.3010759E+03_r8,&  
         0.3005049E+03_r8,     0.2999347E+03_r8,     0.2993651E+03_r8,     0.2987961E+03_r8,     0.2982270E+03_r8,&  
         0.2976591E+03_r8,     0.2970912E+03_r8,     0.2965240E+03_r8,     0.2959569E+03_r8,     0.2953902E+03_r8,&  
         0.2948244E+03_r8,     0.2942582E+03_r8,     0.2936927E+03_r8,     0.2931270E+03_r8,     0.2925618E+03_r8,&  
         0.2919971E+03_r8,     0.2914319E+03_r8,     0.2908672E+03_r8,     0.2903023E+03_r8,     0.2897377E+03_r8,&  
         0.2891729E+03_r8,     0.2886080E+03_r8,     0.2880428E+03_r8,     0.2874785E+03_r8,     0.2869130E+03_r8,&  
         0.2863475E+03_r8,     0.2857819E+03_r8,     0.2852160E+03_r8,     0.2846500E+03_r8,     0.2840831E+03_r8,&  
         0.2835160E+03_r8,     0.2829489E+03_r8,     0.2823804E+03_r8,     0.2818120E+03_r8,     0.2812426E+03_r8,&  
         0.2806727E+03_r8,     0.2801018E+03_r8,     0.2795304E+03_r8,     0.2789580E+03_r8,     0.2783849E+03_r8,&  
         0.2778105E+03_r8,     0.2772352E+03_r8,     0.2766584E+03_r8,     0.2760811E+03_r8,     0.2755020E+03_r8,&  
         0.2749219E+03_r8,     0.2743401E+03_r8,     0.2737570E+03_r8,     0.2731722E+03_r8,     0.2725858E+03_r8,&  
         0.2719976E+03_r8,     0.2714079E+03_r8,     0.2708159E+03_r8,     0.2702218E+03_r8,     0.2696260E+03_r8,&  
         0.2690274E+03_r8,     0.2684271E+03_r8,     0.2678235E+03_r8,     0.2672182E+03_r8,     0.2666095E+03_r8,&  
         0.2659983E+03_r8,     0.2653840E+03_r8,     0.2647666E+03_r8,     0.2641462E+03_r8,     0.2635220E+03_r8,&  
         0.2628945E+03_r8,     0.2622632E+03_r8,     0.2616280E+03_r8,     0.2609890E+03_r8,     0.2603452E+03_r8,&  
         0.2596974E+03_r8,     0.2590447E+03_r8,     0.2583870E+03_r8,     0.2577249E+03_r8,     0.2570569E+03_r8,&  
         0.2563835E+03_r8,     0.2557042E+03_r8,     0.2550190E+03_r8,     0.2543275E+03_r8,     0.2536291E+03_r8,&  
         0.2529241E+03_r8,     0.2522118E+03_r8,     0.2514921E+03_r8,     0.2507640E+03_r8,     0.2500281E+03_r8,&  
         0.2492837E+03_r8,     0.2485302E+03_r8,     0.2477675E+03_r8,     0.2469950E+03_r8,     0.2462124E+03_r8,&  
         0.2454192E+03_r8,     0.2446150E+03_r8,     0.2437993E+03_r8,     0.2429719E+03_r8,     0.2421322E+03_r8,&  
         0.2412796E+03_r8,     0.2404137E+03_r8,     0.2395341E+03_r8,     0.2386404E+03_r8,     0.2377321E+03_r8,&  
         0.2368089E+03_r8,     0.2358703E+03_r8,     0.2349160E+03_r8,     0.2339456E+03_r8,     0.2329590E+03_r8,&  
         0.2319559E+03_r8,     0.2309361E+03_r8,     0.2298997E+03_r8,     0.2288464E+03_r8,     0.2277764E+03_r8,&  
         0.2266904E+03_r8,     0.2255880E+03_r8,     0.2244701E+03_r8,     0.2233371E+03_r8,     0.2221895E+03_r8/)  

    psaditmk(1:150,  127)= (/ &
         0.3212245E+03_r8,     0.3206209E+03_r8,     0.3200195E+03_r8,     0.3194191E+03_r8,     0.3188203E+03_r8,&  
         0.3182223E+03_r8,     0.3176261E+03_r8,     0.3170307E+03_r8,     0.3164363E+03_r8,     0.3158435E+03_r8,&  
         0.3152517E+03_r8,     0.3146616E+03_r8,     0.3140718E+03_r8,     0.3134836E+03_r8,     0.3128961E+03_r8,&  
         0.3123103E+03_r8,     0.3117250E+03_r8,     0.3111405E+03_r8,     0.3105576E+03_r8,     0.3099755E+03_r8,&  
         0.3093945E+03_r8,     0.3088144E+03_r8,     0.3082354E+03_r8,     0.3076567E+03_r8,     0.3070795E+03_r8,&  
         0.3065031E+03_r8,     0.3059272E+03_r8,     0.3053525E+03_r8,     0.3047778E+03_r8,     0.3042047E+03_r8,&  
         0.3036325E+03_r8,     0.3030608E+03_r8,     0.3024893E+03_r8,     0.3019189E+03_r8,     0.3013494E+03_r8,&  
         0.3007802E+03_r8,     0.3002115E+03_r8,     0.2996437E+03_r8,     0.2990763E+03_r8,     0.2985096E+03_r8,&  
         0.2979430E+03_r8,     0.2973769E+03_r8,     0.2968118E+03_r8,     0.2962468E+03_r8,     0.2956822E+03_r8,&  
         0.2951180E+03_r8,     0.2945542E+03_r8,     0.2939908E+03_r8,     0.2934274E+03_r8,     0.2928642E+03_r8,&  
         0.2923012E+03_r8,     0.2917389E+03_r8,     0.2911766E+03_r8,     0.2906139E+03_r8,     0.2900518E+03_r8,&  
         0.2894894E+03_r8,     0.2889272E+03_r8,     0.2883647E+03_r8,     0.2878025E+03_r8,     0.2872399E+03_r8,&  
         0.2866775E+03_r8,     0.2861149E+03_r8,     0.2855516E+03_r8,     0.2849881E+03_r8,     0.2844244E+03_r8,&  
         0.2838607E+03_r8,     0.2832962E+03_r8,     0.2827315E+03_r8,     0.2821661E+03_r8,     0.2816002E+03_r8,&  
         0.2810334E+03_r8,     0.2804666E+03_r8,     0.2798983E+03_r8,     0.2793299E+03_r8,     0.2787603E+03_r8,&  
         0.2781900E+03_r8,     0.2776186E+03_r8,     0.2770463E+03_r8,     0.2764730E+03_r8,     0.2758979E+03_r8,&  
         0.2753224E+03_r8,     0.2747451E+03_r8,     0.2741668E+03_r8,     0.2735866E+03_r8,     0.2730054E+03_r8,&  
         0.2724219E+03_r8,     0.2718374E+03_r8,     0.2712508E+03_r8,     0.2706624E+03_r8,     0.2700722E+03_r8,&  
         0.2694793E+03_r8,     0.2688849E+03_r8,     0.2682875E+03_r8,     0.2676885E+03_r8,     0.2670861E+03_r8,&  
         0.2664819E+03_r8,     0.2658744E+03_r8,     0.2652641E+03_r8,     0.2646510E+03_r8,     0.2640345E+03_r8,&  
         0.2634149E+03_r8,     0.2627915E+03_r8,     0.2621648E+03_r8,     0.2615341E+03_r8,     0.2608994E+03_r8,&  
         0.2602607E+03_r8,     0.2596175E+03_r8,     0.2589699E+03_r8,     0.2583175E+03_r8,     0.2576602E+03_r8,&  
         0.2569980E+03_r8,     0.2563297E+03_r8,     0.2556560E+03_r8,     0.2549769E+03_r8,     0.2542910E+03_r8,&  
         0.2535986E+03_r8,     0.2528998E+03_r8,     0.2521940E+03_r8,     0.2514806E+03_r8,     0.2507594E+03_r8,&  
         0.2500305E+03_r8,     0.2492933E+03_r8,     0.2485472E+03_r8,     0.2477920E+03_r8,     0.2470273E+03_r8,&  
         0.2462524E+03_r8,     0.2454674E+03_r8,     0.2446717E+03_r8,     0.2438650E+03_r8,     0.2430463E+03_r8,&  
         0.2422156E+03_r8,     0.2413725E+03_r8,     0.2405160E+03_r8,     0.2396465E+03_r8,     0.2387627E+03_r8,&  
         0.2378647E+03_r8,     0.2369520E+03_r8,     0.2360240E+03_r8,     0.2350805E+03_r8,     0.2341212E+03_r8,&  
         0.2331456E+03_r8,     0.2321534E+03_r8,     0.2311447E+03_r8,     0.2301193E+03_r8,     0.2290771E+03_r8,&  
         0.2280183E+03_r8,     0.2269427E+03_r8,     0.2258511E+03_r8,     0.2247434E+03_r8,     0.2236202E+03_r8/)  

    psaditmk(1:150,  128)= (/ &
         0.3214514E+03_r8,     0.3208495E+03_r8,     0.3202490E+03_r8,     0.3196495E+03_r8,     0.3190518E+03_r8,&  
         0.3184547E+03_r8,     0.3178591E+03_r8,     0.3172652E+03_r8,     0.3166719E+03_r8,     0.3160800E+03_r8,&  
         0.3154895E+03_r8,     0.3149002E+03_r8,     0.3143118E+03_r8,     0.3137245E+03_r8,     0.3131385E+03_r8,&  
         0.3125535E+03_r8,     0.3119696E+03_r8,     0.3113864E+03_r8,     0.3108044E+03_r8,     0.3102235E+03_r8,&  
         0.3096438E+03_r8,     0.3090647E+03_r8,     0.3084873E+03_r8,     0.3079100E+03_r8,     0.3073340E+03_r8,&  
         0.3067586E+03_r8,     0.3061844E+03_r8,     0.3056113E+03_r8,     0.3050386E+03_r8,     0.3044662E+03_r8,&  
         0.3038952E+03_r8,     0.3033253E+03_r8,     0.3027558E+03_r8,     0.3021869E+03_r8,     0.3016184E+03_r8,&  
         0.3010509E+03_r8,     0.3004844E+03_r8,     0.2999181E+03_r8,     0.2993523E+03_r8,     0.2987871E+03_r8,&  
         0.2982229E+03_r8,     0.2976585E+03_r8,     0.2970951E+03_r8,     0.2965322E+03_r8,     0.2959695E+03_r8,&  
         0.2954073E+03_r8,     0.2948451E+03_r8,     0.2942838E+03_r8,     0.2937226E+03_r8,     0.2931618E+03_r8,&  
         0.2926010E+03_r8,     0.2920405E+03_r8,     0.2914801E+03_r8,     0.2909204E+03_r8,     0.2903606E+03_r8,&  
         0.2898004E+03_r8,     0.2892406E+03_r8,     0.2886811E+03_r8,     0.2881213E+03_r8,     0.2875612E+03_r8,&  
         0.2870013E+03_r8,     0.2864410E+03_r8,     0.2858811E+03_r8,     0.2853207E+03_r8,     0.2847598E+03_r8,&  
         0.2841984E+03_r8,     0.2836375E+03_r8,     0.2830756E+03_r8,     0.2825135E+03_r8,     0.2819511E+03_r8,&  
         0.2813875E+03_r8,     0.2808238E+03_r8,     0.2802592E+03_r8,     0.2796943E+03_r8,     0.2791284E+03_r8,&  
         0.2785620E+03_r8,     0.2779942E+03_r8,     0.2774258E+03_r8,     0.2768562E+03_r8,     0.2762859E+03_r8,&  
         0.2757143E+03_r8,     0.2751417E+03_r8,     0.2745676E+03_r8,     0.2739919E+03_r8,     0.2734154E+03_r8,&  
         0.2728367E+03_r8,     0.2722572E+03_r8,     0.2716757E+03_r8,     0.2710925E+03_r8,     0.2705074E+03_r8,&  
         0.2699206E+03_r8,     0.2693317E+03_r8,     0.2687401E+03_r8,     0.2681472E+03_r8,     0.2675509E+03_r8,&  
         0.2669531E+03_r8,     0.2663521E+03_r8,     0.2657488E+03_r8,     0.2651425E+03_r8,     0.2645333E+03_r8,&  
         0.2639211E+03_r8,     0.2633053E+03_r8,     0.2626867E+03_r8,     0.2620641E+03_r8,     0.2614378E+03_r8,&  
         0.2608079E+03_r8,     0.2601737E+03_r8,     0.2595356E+03_r8,     0.2588925E+03_r8,     0.2582451E+03_r8,&  
         0.2575929E+03_r8,     0.2569356E+03_r8,     0.2562729E+03_r8,     0.2556051E+03_r8,     0.2549310E+03_r8,&  
         0.2542512E+03_r8,     0.2535649E+03_r8,     0.2528721E+03_r8,     0.2521725E+03_r8,     0.2514659E+03_r8,&  
         0.2507512E+03_r8,     0.2500290E+03_r8,     0.2492988E+03_r8,     0.2485600E+03_r8,     0.2478122E+03_r8,&  
         0.2470548E+03_r8,     0.2462880E+03_r8,     0.2455112E+03_r8,     0.2447234E+03_r8,     0.2439252E+03_r8,&  
         0.2431154E+03_r8,     0.2422937E+03_r8,     0.2414595E+03_r8,     0.2406130E+03_r8,     0.2397529E+03_r8,&  
         0.2388792E+03_r8,     0.2379914E+03_r8,     0.2370890E+03_r8,     0.2361716E+03_r8,     0.2352387E+03_r8,&  
         0.2342902E+03_r8,     0.2333255E+03_r8,     0.2323446E+03_r8,     0.2313469E+03_r8,     0.2303326E+03_r8,&  
         0.2293015E+03_r8,     0.2282537E+03_r8,     0.2271892E+03_r8,     0.2261080E+03_r8,     0.2250108E+03_r8/)  

    psaditmk(1:150,  129)= (/ &
         0.3216765E+03_r8,     0.3210758E+03_r8,     0.3204755E+03_r8,     0.3198774E+03_r8,     0.3192802E+03_r8,&  
         0.3186845E+03_r8,     0.3180895E+03_r8,     0.3174969E+03_r8,     0.3169049E+03_r8,     0.3163135E+03_r8,&  
         0.3157242E+03_r8,     0.3151358E+03_r8,     0.3145483E+03_r8,     0.3139623E+03_r8,     0.3133773E+03_r8,&  
         0.3127934E+03_r8,     0.3122105E+03_r8,     0.3116291E+03_r8,     0.3110482E+03_r8,     0.3104684E+03_r8,&  
         0.3098897E+03_r8,     0.3093121E+03_r8,     0.3087353E+03_r8,     0.3081602E+03_r8,     0.3075853E+03_r8,&  
         0.3070115E+03_r8,     0.3064383E+03_r8,     0.3058663E+03_r8,     0.3052953E+03_r8,     0.3047247E+03_r8,&  
         0.3041551E+03_r8,     0.3035860E+03_r8,     0.3030181E+03_r8,     0.3024510E+03_r8,     0.3018843E+03_r8,&  
         0.3013182E+03_r8,     0.3007532E+03_r8,     0.3001885E+03_r8,     0.2996246E+03_r8,     0.2990610E+03_r8,&  
         0.2984980E+03_r8,     0.2979362E+03_r8,     0.2973742E+03_r8,     0.2968130E+03_r8,     0.2962522E+03_r8,&  
         0.2956917E+03_r8,     0.2951318E+03_r8,     0.2945725E+03_r8,     0.2940132E+03_r8,     0.2934542E+03_r8,&  
         0.2928959E+03_r8,     0.2923378E+03_r8,     0.2917794E+03_r8,     0.2912215E+03_r8,     0.2906636E+03_r8,&  
         0.2901065E+03_r8,     0.2895490E+03_r8,     0.2889914E+03_r8,     0.2884341E+03_r8,     0.2878767E+03_r8,&  
         0.2873196E+03_r8,     0.2867620E+03_r8,     0.2862043E+03_r8,     0.2856466E+03_r8,     0.2850888E+03_r8,&  
         0.2845309E+03_r8,     0.2839722E+03_r8,     0.2834134E+03_r8,     0.2828546E+03_r8,     0.2822947E+03_r8,&  
         0.2817349E+03_r8,     0.2811742E+03_r8,     0.2806132E+03_r8,     0.2800513E+03_r8,     0.2794893E+03_r8,&  
         0.2789262E+03_r8,     0.2783621E+03_r8,     0.2777978E+03_r8,     0.2772321E+03_r8,     0.2766657E+03_r8,&  
         0.2760980E+03_r8,     0.2755295E+03_r8,     0.2749597E+03_r8,     0.2743887E+03_r8,     0.2738164E+03_r8,&  
         0.2732427E+03_r8,     0.2726677E+03_r8,     0.2720910E+03_r8,     0.2715128E+03_r8,     0.2709326E+03_r8,&  
         0.2703511E+03_r8,     0.2697674E+03_r8,     0.2691817E+03_r8,     0.2685943E+03_r8,     0.2680044E+03_r8,&  
         0.2674124E+03_r8,     0.2668180E+03_r8,     0.2662210E+03_r8,     0.2656214E+03_r8,     0.2650190E+03_r8,&  
         0.2644138E+03_r8,     0.2638053E+03_r8,     0.2631941E+03_r8,     0.2625793E+03_r8,     0.2619610E+03_r8,&  
         0.2613394E+03_r8,     0.2607137E+03_r8,     0.2600843E+03_r8,     0.2594505E+03_r8,     0.2588128E+03_r8,&  
         0.2581701E+03_r8,     0.2575230E+03_r8,     0.2568708E+03_r8,     0.2562135E+03_r8,     0.2555509E+03_r8,&  
         0.2548824E+03_r8,     0.2542085E+03_r8,     0.2535281E+03_r8,     0.2528410E+03_r8,     0.2521475E+03_r8,&  
         0.2514474E+03_r8,     0.2507391E+03_r8,     0.2500237E+03_r8,     0.2493003E+03_r8,     0.2485687E+03_r8,&  
         0.2478278E+03_r8,     0.2470782E+03_r8,     0.2463190E+03_r8,     0.2455502E+03_r8,     0.2447709E+03_r8,&  
         0.2439808E+03_r8,     0.2431796E+03_r8,     0.2423670E+03_r8,     0.2415419E+03_r8,     0.2407044E+03_r8,&  
         0.2398542E+03_r8,     0.2389901E+03_r8,     0.2381122E+03_r8,     0.2372200E+03_r8,     0.2363130E+03_r8,&  
         0.2353908E+03_r8,     0.2344530E+03_r8,     0.2334992E+03_r8,     0.2325294E+03_r8,     0.2315429E+03_r8,&  
         0.2305397E+03_r8,     0.2295197E+03_r8,     0.2284830E+03_r8,     0.2274294E+03_r8,     0.2263591E+03_r8/)  

    psaditmk(1:150,  130)= (/ &
         0.3218990E+03_r8,     0.3212989E+03_r8,     0.3207000E+03_r8,     0.3201022E+03_r8,     0.3195060E+03_r8,&  
         0.3189113E+03_r8,     0.3183179E+03_r8,     0.3177255E+03_r8,     0.3171341E+03_r8,     0.3165443E+03_r8,&  
         0.3159558E+03_r8,     0.3153684E+03_r8,     0.3147823E+03_r8,     0.3141973E+03_r8,     0.3136135E+03_r8,&  
         0.3130304E+03_r8,     0.3124487E+03_r8,     0.3118680E+03_r8,     0.3112885E+03_r8,     0.3107103E+03_r8,&  
         0.3101331E+03_r8,     0.3095563E+03_r8,     0.3089808E+03_r8,     0.3084063E+03_r8,     0.3078330E+03_r8,&  
         0.3072605E+03_r8,     0.3066886E+03_r8,     0.3061179E+03_r8,     0.3055482E+03_r8,     0.3049788E+03_r8,&  
         0.3044109E+03_r8,     0.3038435E+03_r8,     0.3032767E+03_r8,     0.3027108E+03_r8,     0.3021460E+03_r8,&  
         0.3015817E+03_r8,     0.3010179E+03_r8,     0.3004548E+03_r8,     0.2998923E+03_r8,     0.2993309E+03_r8,&  
         0.2987696E+03_r8,     0.2982090E+03_r8,     0.2976491E+03_r8,     0.2970896E+03_r8,     0.2965308E+03_r8,&  
         0.2959722E+03_r8,     0.2954141E+03_r8,     0.2948565E+03_r8,     0.2942995E+03_r8,     0.2937425E+03_r8,&  
         0.2931858E+03_r8,     0.2926297E+03_r8,     0.2920736E+03_r8,     0.2915182E+03_r8,     0.2909623E+03_r8,&  
         0.2904070E+03_r8,     0.2898520E+03_r8,     0.2892971E+03_r8,     0.2887421E+03_r8,     0.2881870E+03_r8,&  
         0.2876322E+03_r8,     0.2870773E+03_r8,     0.2865223E+03_r8,     0.2859673E+03_r8,     0.2854121E+03_r8,&  
         0.2848563E+03_r8,     0.2843008E+03_r8,     0.2837451E+03_r8,     0.2831887E+03_r8,     0.2826324E+03_r8,&  
         0.2820754E+03_r8,     0.2815181E+03_r8,     0.2809604E+03_r8,     0.2804017E+03_r8,     0.2798429E+03_r8,&  
         0.2792831E+03_r8,     0.2787231E+03_r8,     0.2781618E+03_r8,     0.2776000E+03_r8,     0.2770372E+03_r8,&  
         0.2764738E+03_r8,     0.2759091E+03_r8,     0.2753431E+03_r8,     0.2747768E+03_r8,     0.2742083E+03_r8,&  
         0.2736393E+03_r8,     0.2730685E+03_r8,     0.2724969E+03_r8,     0.2719236E+03_r8,     0.2713484E+03_r8,&  
         0.2707718E+03_r8,     0.2701933E+03_r8,     0.2696131E+03_r8,     0.2690308E+03_r8,     0.2684467E+03_r8,&  
         0.2678605E+03_r8,     0.2672717E+03_r8,     0.2666812E+03_r8,     0.2660879E+03_r8,     0.2654918E+03_r8,&  
         0.2648936E+03_r8,     0.2642921E+03_r8,     0.2636880E+03_r8,     0.2630804E+03_r8,     0.2624700E+03_r8,&  
         0.2618561E+03_r8,     0.2612384E+03_r8,     0.2606175E+03_r8,     0.2599925E+03_r8,     0.2593634E+03_r8,&  
         0.2587302E+03_r8,     0.2580925E+03_r8,     0.2574502E+03_r8,     0.2568029E+03_r8,     0.2561510E+03_r8,&  
         0.2554939E+03_r8,     0.2548310E+03_r8,     0.2541622E+03_r8,     0.2534878E+03_r8,     0.2528065E+03_r8,&  
         0.2521193E+03_r8,     0.2514252E+03_r8,     0.2507234E+03_r8,     0.2500145E+03_r8,     0.2492980E+03_r8,&  
         0.2485730E+03_r8,     0.2478395E+03_r8,     0.2470974E+03_r8,     0.2463458E+03_r8,     0.2455847E+03_r8,&  
         0.2448133E+03_r8,     0.2440315E+03_r8,     0.2432388E+03_r8,     0.2424346E+03_r8,     0.2416186E+03_r8,&  
         0.2407903E+03_r8,     0.2399494E+03_r8,     0.2390952E+03_r8,     0.2382273E+03_r8,     0.2373452E+03_r8,&  
         0.2364485E+03_r8,     0.2355368E+03_r8,     0.2346097E+03_r8,     0.2336668E+03_r8,     0.2327078E+03_r8,&  
         0.2317324E+03_r8,     0.2307403E+03_r8,     0.2297314E+03_r8,     0.2287058E+03_r8,     0.2276632E+03_r8/)  

    psaditmk(1:150,  131)= (/ &
         0.3221185E+03_r8,     0.3215196E+03_r8,     0.3209214E+03_r8,     0.3203247E+03_r8,     0.3197292E+03_r8,&  
         0.3191355E+03_r8,     0.3185429E+03_r8,     0.3179514E+03_r8,     0.3173614E+03_r8,     0.3167725E+03_r8,&  
         0.3161847E+03_r8,     0.3155985E+03_r8,     0.3150135E+03_r8,     0.3144294E+03_r8,     0.3138463E+03_r8,&  
         0.3132643E+03_r8,     0.3126837E+03_r8,     0.3121044E+03_r8,     0.3115259E+03_r8,     0.3109488E+03_r8,&  
         0.3103725E+03_r8,     0.3097974E+03_r8,     0.3092233E+03_r8,     0.3086495E+03_r8,     0.3080773E+03_r8,&  
         0.3075059E+03_r8,     0.3069358E+03_r8,     0.3063665E+03_r8,     0.3057979E+03_r8,     0.3052302E+03_r8,&  
         0.3046634E+03_r8,     0.3040975E+03_r8,     0.3035322E+03_r8,     0.3029680E+03_r8,     0.3024041E+03_r8,&  
         0.3018411E+03_r8,     0.3012791E+03_r8,     0.3007178E+03_r8,     0.3001568E+03_r8,     0.2995963E+03_r8,&  
         0.2990375E+03_r8,     0.2984783E+03_r8,     0.2979197E+03_r8,     0.2973618E+03_r8,     0.2968051E+03_r8,&  
         0.2962482E+03_r8,     0.2956920E+03_r8,     0.2951363E+03_r8,     0.2945809E+03_r8,     0.2940260E+03_r8,&  
         0.2934716E+03_r8,     0.2929171E+03_r8,     0.2923633E+03_r8,     0.2918095E+03_r8,     0.2912563E+03_r8,&  
         0.2907034E+03_r8,     0.2901501E+03_r8,     0.2895974E+03_r8,     0.2890445E+03_r8,     0.2884922E+03_r8,&  
         0.2879397E+03_r8,     0.2873872E+03_r8,     0.2868344E+03_r8,     0.2862820E+03_r8,     0.2857296E+03_r8,&  
         0.2851767E+03_r8,     0.2846240E+03_r8,     0.2840704E+03_r8,     0.2835175E+03_r8,     0.2829638E+03_r8,&  
         0.2824096E+03_r8,     0.2818555E+03_r8,     0.2813006E+03_r8,     0.2807454E+03_r8,     0.2801894E+03_r8,&  
         0.2796336E+03_r8,     0.2790764E+03_r8,     0.2785189E+03_r8,     0.2779610E+03_r8,     0.2774016E+03_r8,&  
         0.2768417E+03_r8,     0.2762808E+03_r8,     0.2757191E+03_r8,     0.2751563E+03_r8,     0.2745925E+03_r8,&  
         0.2740276E+03_r8,     0.2734614E+03_r8,     0.2728939E+03_r8,     0.2723250E+03_r8,     0.2717546E+03_r8,&  
         0.2711830E+03_r8,     0.2706093E+03_r8,     0.2700343E+03_r8,     0.2694572E+03_r8,     0.2688786E+03_r8,&  
         0.2682978E+03_r8,     0.2677148E+03_r8,     0.2671300E+03_r8,     0.2665427E+03_r8,     0.2659532E+03_r8,&  
         0.2653611E+03_r8,     0.2647662E+03_r8,     0.2641689E+03_r8,     0.2635684E+03_r8,     0.2629654E+03_r8,&  
         0.2623587E+03_r8,     0.2617488E+03_r8,     0.2611358E+03_r8,     0.2605188E+03_r8,     0.2598985E+03_r8,&  
         0.2592739E+03_r8,     0.2586453E+03_r8,     0.2580125E+03_r8,     0.2573751E+03_r8,     0.2567329E+03_r8,&  
         0.2560860E+03_r8,     0.2554337E+03_r8,     0.2547764E+03_r8,     0.2541131E+03_r8,     0.2534443E+03_r8,&  
         0.2527691E+03_r8,     0.2520875E+03_r8,     0.2513997E+03_r8,     0.2507041E+03_r8,     0.2500018E+03_r8,&  
         0.2492920E+03_r8,     0.2485735E+03_r8,     0.2478473E+03_r8,     0.2471124E+03_r8,     0.2463684E+03_r8,&  
         0.2456146E+03_r8,     0.2448512E+03_r8,     0.2440774E+03_r8,     0.2432932E+03_r8,     0.2424976E+03_r8,&  
         0.2416906E+03_r8,     0.2408713E+03_r8,     0.2400395E+03_r8,     0.2391949E+03_r8,     0.2383366E+03_r8,&  
         0.2374645E+03_r8,     0.2365780E+03_r8,     0.2356766E+03_r8,     0.2347602E+03_r8,     0.2338281E+03_r8,&  
         0.2328800E+03_r8,     0.2319156E+03_r8,     0.2309347E+03_r8,     0.2299370E+03_r8,     0.2289223E+03_r8/)  

    psaditmk(1:150,  132)= (/ &
         0.3223362E+03_r8,     0.3217372E+03_r8,     0.3211405E+03_r8,     0.3205443E+03_r8,     0.3199502E+03_r8,&  
         0.3193571E+03_r8,     0.3187653E+03_r8,     0.3181750E+03_r8,     0.3175855E+03_r8,     0.3169979E+03_r8,&  
         0.3164109E+03_r8,     0.3158253E+03_r8,     0.3152415E+03_r8,     0.3146585E+03_r8,     0.3140765E+03_r8,&  
         0.3134954E+03_r8,     0.3129160E+03_r8,     0.3123375E+03_r8,     0.3117603E+03_r8,     0.3111839E+03_r8,&  
         0.3106090E+03_r8,     0.3100350E+03_r8,     0.3094619E+03_r8,     0.3088900E+03_r8,     0.3083192E+03_r8,&  
         0.3077487E+03_r8,     0.3071794E+03_r8,     0.3066113E+03_r8,     0.3060443E+03_r8,     0.3054778E+03_r8,&  
         0.3049119E+03_r8,     0.3043474E+03_r8,     0.3037838E+03_r8,     0.3032206E+03_r8,     0.3026589E+03_r8,&  
         0.3020971E+03_r8,     0.3015365E+03_r8,     0.3009766E+03_r8,     0.3004174E+03_r8,     0.2998589E+03_r8,&  
         0.2993008E+03_r8,     0.2987432E+03_r8,     0.2981867E+03_r8,     0.2976306E+03_r8,     0.2970753E+03_r8,&  
         0.2965202E+03_r8,     0.2959657E+03_r8,     0.2954117E+03_r8,     0.2948584E+03_r8,     0.2943051E+03_r8,&  
         0.2937524E+03_r8,     0.2932003E+03_r8,     0.2926483E+03_r8,     0.2920965E+03_r8,     0.2915453E+03_r8,&  
         0.2909941E+03_r8,     0.2904434E+03_r8,     0.2898929E+03_r8,     0.2893422E+03_r8,     0.2887920E+03_r8,&  
         0.2882417E+03_r8,     0.2876918E+03_r8,     0.2871417E+03_r8,     0.2865916E+03_r8,     0.2860413E+03_r8,&  
         0.2854912E+03_r8,     0.2849410E+03_r8,     0.2843908E+03_r8,     0.2838399E+03_r8,     0.2832889E+03_r8,&  
         0.2827383E+03_r8,     0.2821865E+03_r8,     0.2816349E+03_r8,     0.2810826E+03_r8,     0.2805299E+03_r8,&  
         0.2799769E+03_r8,     0.2794231E+03_r8,     0.2788690E+03_r8,     0.2783141E+03_r8,     0.2777587E+03_r8,&  
         0.2772020E+03_r8,     0.2766451E+03_r8,     0.2760873E+03_r8,     0.2755282E+03_r8,     0.2749684E+03_r8,&  
         0.2744073E+03_r8,     0.2738457E+03_r8,     0.2732825E+03_r8,     0.2727177E+03_r8,     0.2721519E+03_r8,&  
         0.2715847E+03_r8,     0.2710160E+03_r8,     0.2704458E+03_r8,     0.2698737E+03_r8,     0.2693004E+03_r8,&  
         0.2687246E+03_r8,     0.2681475E+03_r8,     0.2675680E+03_r8,     0.2669865E+03_r8,     0.2664029E+03_r8,&  
         0.2658167E+03_r8,     0.2652283E+03_r8,     0.2646375E+03_r8,     0.2640436E+03_r8,     0.2634474E+03_r8,&  
         0.2628477E+03_r8,     0.2622454E+03_r8,     0.2616399E+03_r8,     0.2610307E+03_r8,     0.2604182E+03_r8,&  
         0.2598021E+03_r8,     0.2591820E+03_r8,     0.2585580E+03_r8,     0.2579299E+03_r8,     0.2572974E+03_r8,&  
         0.2566601E+03_r8,     0.2560183E+03_r8,     0.2553712E+03_r8,     0.2547188E+03_r8,     0.2540611E+03_r8,&  
         0.2533976E+03_r8,     0.2527286E+03_r8,     0.2520525E+03_r8,     0.2513706E+03_r8,     0.2506813E+03_r8,&  
         0.2499854E+03_r8,     0.2492821E+03_r8,     0.2485704E+03_r8,     0.2478510E+03_r8,     0.2471233E+03_r8,&  
         0.2463863E+03_r8,     0.2456403E+03_r8,     0.2448847E+03_r8,     0.2441190E+03_r8,     0.2433427E+03_r8,&  
         0.2425558E+03_r8,     0.2417572E+03_r8,     0.2409468E+03_r8,     0.2401242E+03_r8,     0.2392889E+03_r8,&  
         0.2384403E+03_r8,     0.2375781E+03_r8,     0.2367017E+03_r8,     0.2358107E+03_r8,     0.2349046E+03_r8,&  
         0.2339832E+03_r8,     0.2330459E+03_r8,     0.2320925E+03_r8,     0.2311226E+03_r8,     0.2301361E+03_r8/)  

    psaditmk(1:150,  133)= (/ &
         0.3225504E+03_r8,     0.3219531E+03_r8,     0.3213569E+03_r8,     0.3207617E+03_r8,     0.3201681E+03_r8,&  
         0.3195757E+03_r8,     0.3189850E+03_r8,     0.3183953E+03_r8,     0.3178072E+03_r8,     0.3172203E+03_r8,&  
         0.3166344E+03_r8,     0.3160499E+03_r8,     0.3154668E+03_r8,     0.3148845E+03_r8,     0.3143039E+03_r8,&  
         0.3137240E+03_r8,     0.3131454E+03_r8,     0.3125680E+03_r8,     0.3119919E+03_r8,     0.3114165E+03_r8,&  
         0.3108429E+03_r8,     0.3102697E+03_r8,     0.3096978E+03_r8,     0.3091273E+03_r8,     0.3085569E+03_r8,&  
         0.3079882E+03_r8,     0.3074200E+03_r8,     0.3068529E+03_r8,     0.3062871E+03_r8,     0.3057220E+03_r8,&  
         0.3051579E+03_r8,     0.3045946E+03_r8,     0.3040319E+03_r8,     0.3034702E+03_r8,     0.3029097E+03_r8,&  
         0.3023496E+03_r8,     0.3017906E+03_r8,     0.3012318E+03_r8,     0.3006742E+03_r8,     0.3001172E+03_r8,&  
         0.2995608E+03_r8,     0.2990045E+03_r8,     0.2984497E+03_r8,     0.2978953E+03_r8,     0.2973416E+03_r8,&  
         0.2967881E+03_r8,     0.2962354E+03_r8,     0.2956831E+03_r8,     0.2951314E+03_r8,     0.2945801E+03_r8,&  
         0.2940292E+03_r8,     0.2934785E+03_r8,     0.2929286E+03_r8,     0.2923792E+03_r8,     0.2918297E+03_r8,&  
         0.2912804E+03_r8,     0.2907319E+03_r8,     0.2901833E+03_r8,     0.2896352E+03_r8,     0.2890872E+03_r8,&  
         0.2885392E+03_r8,     0.2879911E+03_r8,     0.2874433E+03_r8,     0.2868960E+03_r8,     0.2863481E+03_r8,&  
         0.2858003E+03_r8,     0.2852524E+03_r8,     0.2847043E+03_r8,     0.2841565E+03_r8,     0.2836086E+03_r8,&  
         0.2830601E+03_r8,     0.2825117E+03_r8,     0.2819627E+03_r8,     0.2814134E+03_r8,     0.2808639E+03_r8,&  
         0.2803135E+03_r8,     0.2797636E+03_r8,     0.2792122E+03_r8,     0.2786610E+03_r8,     0.2781086E+03_r8,&  
         0.2775556E+03_r8,     0.2770021E+03_r8,     0.2764477E+03_r8,     0.2758925E+03_r8,     0.2753367E+03_r8,&  
         0.2747794E+03_r8,     0.2742213E+03_r8,     0.2736621E+03_r8,     0.2731020E+03_r8,     0.2725407E+03_r8,&  
         0.2719775E+03_r8,     0.2714136E+03_r8,     0.2708478E+03_r8,     0.2702808E+03_r8,     0.2697122E+03_r8,&  
         0.2691416E+03_r8,     0.2685697E+03_r8,     0.2679951E+03_r8,     0.2674194E+03_r8,     0.2668415E+03_r8,&  
         0.2662612E+03_r8,     0.2656789E+03_r8,     0.2650942E+03_r8,     0.2645066E+03_r8,     0.2639168E+03_r8,&  
         0.2633242E+03_r8,     0.2627287E+03_r8,     0.2621301E+03_r8,     0.2615286E+03_r8,     0.2609240E+03_r8,&  
         0.2603154E+03_r8,     0.2597034E+03_r8,     0.2590882E+03_r8,     0.2584684E+03_r8,     0.2578449E+03_r8,&  
         0.2572169E+03_r8,     0.2565845E+03_r8,     0.2559474E+03_r8,     0.2553056E+03_r8,     0.2546584E+03_r8,&  
         0.2540063E+03_r8,     0.2533481E+03_r8,     0.2526845E+03_r8,     0.2520145E+03_r8,     0.2513384E+03_r8,&  
         0.2506554E+03_r8,     0.2499654E+03_r8,     0.2492686E+03_r8,     0.2485634E+03_r8,     0.2478510E+03_r8,&  
         0.2471302E+03_r8,     0.2464003E+03_r8,     0.2456618E+03_r8,     0.2449137E+03_r8,     0.2441557E+03_r8,&  
         0.2433876E+03_r8,     0.2426087E+03_r8,     0.2418189E+03_r8,     0.2410174E+03_r8,     0.2402037E+03_r8,&  
         0.2393776E+03_r8,     0.2385386E+03_r8,     0.2376860E+03_r8,     0.2368195E+03_r8,     0.2359389E+03_r8,&  
         0.2350433E+03_r8,     0.2341323E+03_r8,     0.2332057E+03_r8,     0.2322630E+03_r8,     0.2313042E+03_r8/)  

    psaditmk(1:150,  134)= (/ &
         0.3227627E+03_r8,     0.3221658E+03_r8,     0.3215702E+03_r8,     0.3209764E+03_r8,     0.3203833E+03_r8,&  
         0.3197926E+03_r8,     0.3192024E+03_r8,     0.3186134E+03_r8,     0.3180262E+03_r8,     0.3174403E+03_r8,&  
         0.3168553E+03_r8,     0.3162713E+03_r8,     0.3156891E+03_r8,     0.3151082E+03_r8,     0.3145285E+03_r8,&  
         0.3139496E+03_r8,     0.3133720E+03_r8,     0.3127957E+03_r8,     0.3122205E+03_r8,     0.3116463E+03_r8,&  
         0.3110732E+03_r8,     0.3105013E+03_r8,     0.3099308E+03_r8,     0.3093607E+03_r8,     0.3087922E+03_r8,&  
         0.3082245E+03_r8,     0.3076580E+03_r8,     0.3070920E+03_r8,     0.3065268E+03_r8,     0.3059630E+03_r8,&  
         0.3054000E+03_r8,     0.3048383E+03_r8,     0.3042768E+03_r8,     0.3037165E+03_r8,     0.3031570E+03_r8,&  
         0.3025987E+03_r8,     0.3020406E+03_r8,     0.3014839E+03_r8,     0.3009271E+03_r8,     0.3003714E+03_r8,&  
         0.2998168E+03_r8,     0.2992627E+03_r8,     0.2987091E+03_r8,     0.2981560E+03_r8,     0.2976037E+03_r8,&  
         0.2970522E+03_r8,     0.2965011E+03_r8,     0.2959503E+03_r8,     0.2954003E+03_r8,     0.2948506E+03_r8,&  
         0.2943019E+03_r8,     0.2937531E+03_r8,     0.2932049E+03_r8,     0.2926573E+03_r8,     0.2921097E+03_r8,&  
         0.2915629E+03_r8,     0.2910159E+03_r8,     0.2904692E+03_r8,     0.2899230E+03_r8,     0.2893770E+03_r8,&  
         0.2888315E+03_r8,     0.2882856E+03_r8,     0.2877402E+03_r8,     0.2871945E+03_r8,     0.2866493E+03_r8,&  
         0.2861039E+03_r8,     0.2855590E+03_r8,     0.2850133E+03_r8,     0.2844680E+03_r8,     0.2839220E+03_r8,&  
         0.2833768E+03_r8,     0.2828308E+03_r8,     0.2822845E+03_r8,     0.2817384E+03_r8,     0.2811916E+03_r8,&  
         0.2806445E+03_r8,     0.2800968E+03_r8,     0.2795490E+03_r8,     0.2790007E+03_r8,     0.2784517E+03_r8,&  
         0.2779023E+03_r8,     0.2773520E+03_r8,     0.2768011E+03_r8,     0.2762498E+03_r8,     0.2756969E+03_r8,&  
         0.2751437E+03_r8,     0.2745898E+03_r8,     0.2740342E+03_r8,     0.2734781E+03_r8,     0.2729206E+03_r8,&  
         0.2723622E+03_r8,     0.2718023E+03_r8,     0.2712411E+03_r8,     0.2706787E+03_r8,     0.2701146E+03_r8,&  
         0.2695492E+03_r8,     0.2689819E+03_r8,     0.2684128E+03_r8,     0.2678423E+03_r8,     0.2672697E+03_r8,&  
         0.2666947E+03_r8,     0.2661182E+03_r8,     0.2655392E+03_r8,     0.2649579E+03_r8,     0.2643745E+03_r8,&  
         0.2637880E+03_r8,     0.2631992E+03_r8,     0.2626079E+03_r8,     0.2620132E+03_r8,     0.2614155E+03_r8,&  
         0.2608145E+03_r8,     0.2602106E+03_r8,     0.2596029E+03_r8,     0.2589916E+03_r8,     0.2583766E+03_r8,&  
         0.2577576E+03_r8,     0.2571341E+03_r8,     0.2565065E+03_r8,     0.2558745E+03_r8,     0.2552374E+03_r8,&  
         0.2545955E+03_r8,     0.2539483E+03_r8,     0.2532957E+03_r8,     0.2526375E+03_r8,     0.2519732E+03_r8,&  
         0.2513027E+03_r8,     0.2506258E+03_r8,     0.2499420E+03_r8,     0.2492515E+03_r8,     0.2485528E+03_r8,&  
         0.2478470E+03_r8,     0.2471331E+03_r8,     0.2464102E+03_r8,     0.2456790E+03_r8,     0.2449384E+03_r8,&  
         0.2441881E+03_r8,     0.2434277E+03_r8,     0.2426573E+03_r8,     0.2418756E+03_r8,     0.2410828E+03_r8,&  
         0.2402780E+03_r8,     0.2394611E+03_r8,     0.2386313E+03_r8,     0.2377884E+03_r8,     0.2369317E+03_r8,&  
         0.2360609E+03_r8,     0.2351756E+03_r8,     0.2342751E+03_r8,     0.2333593E+03_r8,     0.2324276E+03_r8/)  

    psaditmk(1:150,  135)= (/ &
         0.3229725E+03_r8,     0.3223766E+03_r8,     0.3217819E+03_r8,     0.3211890E+03_r8,     0.3205965E+03_r8,&  
         0.3200063E+03_r8,     0.3194170E+03_r8,     0.3188288E+03_r8,     0.3182428E+03_r8,     0.3176576E+03_r8,&  
         0.3170738E+03_r8,     0.3164905E+03_r8,     0.3159093E+03_r8,     0.3153295E+03_r8,     0.3147503E+03_r8,&  
         0.3141725E+03_r8,     0.3135957E+03_r8,     0.3130205E+03_r8,     0.3124462E+03_r8,     0.3118732E+03_r8,&  
         0.3113012E+03_r8,     0.3107301E+03_r8,     0.3101604E+03_r8,     0.3095919E+03_r8,     0.3090241E+03_r8,&  
         0.3084576E+03_r8,     0.3078922E+03_r8,     0.3073275E+03_r8,     0.3067638E+03_r8,     0.3062010E+03_r8,&  
         0.3056392E+03_r8,     0.3050783E+03_r8,     0.3045184E+03_r8,     0.3039594E+03_r8,     0.3034013E+03_r8,&  
         0.3028439E+03_r8,     0.3022875E+03_r8,     0.3017320E+03_r8,     0.3011772E+03_r8,     0.3006229E+03_r8,&  
         0.3000693E+03_r8,     0.2995164E+03_r8,     0.2989646E+03_r8,     0.2984132E+03_r8,     0.2978621E+03_r8,&  
         0.2973123E+03_r8,     0.2967627E+03_r8,     0.2962137E+03_r8,     0.2956656E+03_r8,     0.2951173E+03_r8,&  
         0.2945702E+03_r8,     0.2940235E+03_r8,     0.2934772E+03_r8,     0.2929312E+03_r8,     0.2923853E+03_r8,&  
         0.2918402E+03_r8,     0.2912954E+03_r8,     0.2907508E+03_r8,     0.2902067E+03_r8,     0.2896626E+03_r8,&  
         0.2891187E+03_r8,     0.2885753E+03_r8,     0.2880321E+03_r8,     0.2874892E+03_r8,     0.2869458E+03_r8,&  
         0.2864024E+03_r8,     0.2858597E+03_r8,     0.2853167E+03_r8,     0.2847738E+03_r8,     0.2842308E+03_r8,&  
         0.2836877E+03_r8,     0.2831441E+03_r8,     0.2826009E+03_r8,     0.2820570E+03_r8,     0.2815132E+03_r8,&  
         0.2809693E+03_r8,     0.2804241E+03_r8,     0.2798796E+03_r8,     0.2793338E+03_r8,     0.2787886E+03_r8,&  
         0.2782420E+03_r8,     0.2776950E+03_r8,     0.2771478E+03_r8,     0.2765993E+03_r8,     0.2760503E+03_r8,&  
         0.2755009E+03_r8,     0.2749500E+03_r8,     0.2743987E+03_r8,     0.2738462E+03_r8,     0.2732928E+03_r8,&  
         0.2727385E+03_r8,     0.2721827E+03_r8,     0.2716259E+03_r8,     0.2710678E+03_r8,     0.2705081E+03_r8,&  
         0.2699475E+03_r8,     0.2693848E+03_r8,     0.2688207E+03_r8,     0.2682550E+03_r8,     0.2676877E+03_r8,&  
         0.2671183E+03_r8,     0.2665471E+03_r8,     0.2659737E+03_r8,     0.2653984E+03_r8,     0.2648207E+03_r8,&  
         0.2642403E+03_r8,     0.2636579E+03_r8,     0.2630727E+03_r8,     0.2624847E+03_r8,     0.2618940E+03_r8,&  
         0.2613005E+03_r8,     0.2607039E+03_r8,     0.2601035E+03_r8,     0.2594999E+03_r8,     0.2588932E+03_r8,&  
         0.2582822E+03_r8,     0.2576676E+03_r8,     0.2570491E+03_r8,     0.2564263E+03_r8,     0.2557987E+03_r8,&  
         0.2551665E+03_r8,     0.2545298E+03_r8,     0.2538876E+03_r8,     0.2532407E+03_r8,     0.2525874E+03_r8,&  
         0.2519291E+03_r8,     0.2512639E+03_r8,     0.2505930E+03_r8,     0.2499151E+03_r8,     0.2492308E+03_r8,&  
         0.2485385E+03_r8,     0.2478393E+03_r8,     0.2471319E+03_r8,     0.2464162E+03_r8,     0.2456921E+03_r8,&  
         0.2449587E+03_r8,     0.2442162E+03_r8,     0.2434637E+03_r8,     0.2427010E+03_r8,     0.2419278E+03_r8,&  
         0.2411431E+03_r8,     0.2403473E+03_r8,     0.2395393E+03_r8,     0.2387186E+03_r8,     0.2378852E+03_r8,&  
         0.2370382E+03_r8,     0.2361775E+03_r8,     0.2353022E+03_r8,     0.2344120E+03_r8,     0.2335067E+03_r8/)  

    psaditmk(1:150,  136)= (/ &
         0.3231807E+03_r8,     0.3225849E+03_r8,     0.3219903E+03_r8,     0.3213983E+03_r8,     0.3208076E+03_r8,&  
         0.3202176E+03_r8,     0.3196290E+03_r8,     0.3190419E+03_r8,     0.3184569E+03_r8,     0.3178724E+03_r8,&  
         0.3172890E+03_r8,     0.3167075E+03_r8,     0.3161264E+03_r8,     0.3155473E+03_r8,     0.3149695E+03_r8,&  
         0.3143927E+03_r8,     0.3138169E+03_r8,     0.3132425E+03_r8,     0.3126695E+03_r8,     0.3120970E+03_r8,&  
         0.3115264E+03_r8,     0.3109565E+03_r8,     0.3103877E+03_r8,     0.3098198E+03_r8,     0.3092533E+03_r8,&  
         0.3086880E+03_r8,     0.3081235E+03_r8,     0.3075601E+03_r8,     0.3069975E+03_r8,     0.3064362E+03_r8,&  
         0.3058757E+03_r8,     0.3053157E+03_r8,     0.3047570E+03_r8,     0.3041991E+03_r8,     0.3036425E+03_r8,&  
         0.3030865E+03_r8,     0.3025311E+03_r8,     0.3019767E+03_r8,     0.3014229E+03_r8,     0.3008703E+03_r8,&  
         0.3003184E+03_r8,     0.2997671E+03_r8,     0.2992162E+03_r8,     0.2986665E+03_r8,     0.2981171E+03_r8,&  
         0.2975689E+03_r8,     0.2970205E+03_r8,     0.2964733E+03_r8,     0.2959268E+03_r8,     0.2953803E+03_r8,&  
         0.2948347E+03_r8,     0.2942892E+03_r8,     0.2937447E+03_r8,     0.2932007E+03_r8,     0.2926570E+03_r8,&  
         0.2921136E+03_r8,     0.2915706E+03_r8,     0.2910277E+03_r8,     0.2904857E+03_r8,     0.2899439E+03_r8,&  
         0.2894019E+03_r8,     0.2888605E+03_r8,     0.2883192E+03_r8,     0.2877780E+03_r8,     0.2872369E+03_r8,&  
         0.2866966E+03_r8,     0.2861557E+03_r8,     0.2856153E+03_r8,     0.2850742E+03_r8,     0.2845335E+03_r8,&  
         0.2839930E+03_r8,     0.2834525E+03_r8,     0.2829113E+03_r8,     0.2823705E+03_r8,     0.2818292E+03_r8,&  
         0.2812876E+03_r8,     0.2807460E+03_r8,     0.2802038E+03_r8,     0.2796617E+03_r8,     0.2791187E+03_r8,&  
         0.2785752E+03_r8,     0.2780318E+03_r8,     0.2774872E+03_r8,     0.2769426E+03_r8,     0.2763970E+03_r8,&  
         0.2758504E+03_r8,     0.2753035E+03_r8,     0.2747557E+03_r8,     0.2742069E+03_r8,     0.2736573E+03_r8,&  
         0.2731067E+03_r8,     0.2725551E+03_r8,     0.2720026E+03_r8,     0.2714482E+03_r8,     0.2708932E+03_r8,&  
         0.2703367E+03_r8,     0.2697788E+03_r8,     0.2692197E+03_r8,     0.2686584E+03_r8,     0.2680958E+03_r8,&  
         0.2675316E+03_r8,     0.2669656E+03_r8,     0.2663973E+03_r8,     0.2658278E+03_r8,     0.2652555E+03_r8,&  
         0.2646812E+03_r8,     0.2641049E+03_r8,     0.2635255E+03_r8,     0.2629444E+03_r8,     0.2623603E+03_r8,&  
         0.2617735E+03_r8,     0.2611834E+03_r8,     0.2605907E+03_r8,     0.2599948E+03_r8,     0.2593954E+03_r8,&  
         0.2587926E+03_r8,     0.2581861E+03_r8,     0.2575761E+03_r8,     0.2569617E+03_r8,     0.2563431E+03_r8,&  
         0.2557205E+03_r8,     0.2550935E+03_r8,     0.2544613E+03_r8,     0.2538243E+03_r8,     0.2531820E+03_r8,&  
         0.2525347E+03_r8,     0.2518816E+03_r8,     0.2512224E+03_r8,     0.2505571E+03_r8,     0.2498852E+03_r8,&  
         0.2492068E+03_r8,     0.2485207E+03_r8,     0.2478279E+03_r8,     0.2471270E+03_r8,     0.2464183E+03_r8,&  
         0.2457013E+03_r8,     0.2449750E+03_r8,     0.2442399E+03_r8,     0.2434950E+03_r8,     0.2427405E+03_r8,&  
         0.2419750E+03_r8,     0.2411991E+03_r8,     0.2404116E+03_r8,     0.2396123E+03_r8,     0.2388008E+03_r8,&  
         0.2379767E+03_r8,     0.2371395E+03_r8,     0.2362881E+03_r8,     0.2354230E+03_r8,     0.2345431E+03_r8/)  

    psaditmk(1:150,  137)= (/ &
         0.3233855E+03_r8,     0.3227918E+03_r8,     0.3221982E+03_r8,     0.3216064E+03_r8,     0.3210159E+03_r8,&  
         0.3204271E+03_r8,     0.3198392E+03_r8,     0.3192526E+03_r8,     0.3186677E+03_r8,     0.3180848E+03_r8,&  
         0.3175022E+03_r8,     0.3169214E+03_r8,     0.3163412E+03_r8,     0.3157629E+03_r8,     0.3151861E+03_r8,&  
         0.3146104E+03_r8,     0.3140353E+03_r8,     0.3134621E+03_r8,     0.3128900E+03_r8,     0.3123187E+03_r8,&  
         0.3117485E+03_r8,     0.3111797E+03_r8,     0.3106123E+03_r8,     0.3100453E+03_r8,     0.3094799E+03_r8,&  
         0.3089153E+03_r8,     0.3083520E+03_r8,     0.3077894E+03_r8,     0.3072279E+03_r8,     0.3066678E+03_r8,&  
         0.3061087E+03_r8,     0.3055501E+03_r8,     0.3049924E+03_r8,     0.3044359E+03_r8,     0.3038799E+03_r8,&  
         0.3033255E+03_r8,     0.3027713E+03_r8,     0.3022184E+03_r8,     0.3016656E+03_r8,     0.3011147E+03_r8,&  
         0.3005640E+03_r8,     0.3000141E+03_r8,     0.2994647E+03_r8,     0.2989163E+03_r8,     0.2983683E+03_r8,&  
         0.2978214E+03_r8,     0.2972750E+03_r8,     0.2967294E+03_r8,     0.2961839E+03_r8,     0.2956390E+03_r8,&  
         0.2950954E+03_r8,     0.2945518E+03_r8,     0.2940088E+03_r8,     0.2934661E+03_r8,     0.2929242E+03_r8,&  
         0.2923826E+03_r8,     0.2918414E+03_r8,     0.2913007E+03_r8,     0.2907602E+03_r8,     0.2902200E+03_r8,&  
         0.2896803E+03_r8,     0.2891412E+03_r8,     0.2886018E+03_r8,     0.2880629E+03_r8,     0.2875237E+03_r8,&  
         0.2869852E+03_r8,     0.2864467E+03_r8,     0.2859084E+03_r8,     0.2853700E+03_r8,     0.2848319E+03_r8,&  
         0.2842934E+03_r8,     0.2837548E+03_r8,     0.2832166E+03_r8,     0.2826781E+03_r8,     0.2821393E+03_r8,&  
         0.2816010E+03_r8,     0.2810616E+03_r8,     0.2805226E+03_r8,     0.2799829E+03_r8,     0.2794424E+03_r8,&  
         0.2789025E+03_r8,     0.2783615E+03_r8,     0.2778203E+03_r8,     0.2772786E+03_r8,     0.2767364E+03_r8,&  
         0.2761932E+03_r8,     0.2756497E+03_r8,     0.2751053E+03_r8,     0.2745603E+03_r8,     0.2740143E+03_r8,&  
         0.2734671E+03_r8,     0.2729198E+03_r8,     0.2723708E+03_r8,     0.2718208E+03_r8,     0.2712701E+03_r8,&  
         0.2707175E+03_r8,     0.2701641E+03_r8,     0.2696092E+03_r8,     0.2690526E+03_r8,     0.2684952E+03_r8,&  
         0.2679356E+03_r8,     0.2673741E+03_r8,     0.2668116E+03_r8,     0.2662471E+03_r8,     0.2656801E+03_r8,&  
         0.2651115E+03_r8,     0.2645407E+03_r8,     0.2639675E+03_r8,     0.2633922E+03_r8,     0.2628141E+03_r8,&  
         0.2622338E+03_r8,     0.2616506E+03_r8,     0.2610650E+03_r8,     0.2604758E+03_r8,     0.2598838E+03_r8,&  
         0.2592883E+03_r8,     0.2586899E+03_r8,     0.2580876E+03_r8,     0.2574818E+03_r8,     0.2568720E+03_r8,&  
         0.2562580E+03_r8,     0.2556402E+03_r8,     0.2550175E+03_r8,     0.2543901E+03_r8,     0.2537585E+03_r8,&  
         0.2531212E+03_r8,     0.2524790E+03_r8,     0.2518311E+03_r8,     0.2511774E+03_r8,     0.2505176E+03_r8,&  
         0.2498519E+03_r8,     0.2491793E+03_r8,     0.2484998E+03_r8,     0.2478130E+03_r8,     0.2471185E+03_r8,&  
         0.2464167E+03_r8,     0.2457065E+03_r8,     0.2449873E+03_r8,     0.2442594E+03_r8,     0.2435220E+03_r8,&  
         0.2427749E+03_r8,     0.2420178E+03_r8,     0.2412499E+03_r8,     0.2404708E+03_r8,     0.2396803E+03_r8,&  
         0.2388778E+03_r8,     0.2380628E+03_r8,     0.2372348E+03_r8,     0.2363935E+03_r8,     0.2355381E+03_r8/)  

    psaditmk(1:150,  138)= (/ &
         0.3235891E+03_r8,     0.3229945E+03_r8,     0.3224026E+03_r8,     0.3218117E+03_r8,     0.3212220E+03_r8,&  
         0.3206332E+03_r8,     0.3200469E+03_r8,     0.3194610E+03_r8,     0.3188771E+03_r8,     0.3182941E+03_r8,&  
         0.3177130E+03_r8,     0.3171330E+03_r8,     0.3165537E+03_r8,     0.3159766E+03_r8,     0.3154004E+03_r8,&  
         0.3148253E+03_r8,     0.3142512E+03_r8,     0.3136791E+03_r8,     0.3131078E+03_r8,     0.3125372E+03_r8,&  
         0.3119684E+03_r8,     0.3114005E+03_r8,     0.3108335E+03_r8,     0.3102682E+03_r8,     0.3097034E+03_r8,&  
         0.3091401E+03_r8,     0.3085774E+03_r8,     0.3080161E+03_r8,     0.3074559E+03_r8,     0.3068965E+03_r8,&  
         0.3063385E+03_r8,     0.3057814E+03_r8,     0.3052247E+03_r8,     0.3046695E+03_r8,     0.3041145E+03_r8,&  
         0.3035612E+03_r8,     0.3030086E+03_r8,     0.3024567E+03_r8,     0.3019058E+03_r8,     0.3013554E+03_r8,&  
         0.3008060E+03_r8,     0.3002575E+03_r8,     0.2997099E+03_r8,     0.2991625E+03_r8,     0.2986162E+03_r8,&  
         0.2980706E+03_r8,     0.2975253E+03_r8,     0.2969811E+03_r8,     0.2964378E+03_r8,     0.2958943E+03_r8,&  
         0.2953519E+03_r8,     0.2948102E+03_r8,     0.2942690E+03_r8,     0.2937280E+03_r8,     0.2931875E+03_r8,&  
         0.2926476E+03_r8,     0.2921082E+03_r8,     0.2915693E+03_r8,     0.2910308E+03_r8,     0.2904926E+03_r8,&  
         0.2899543E+03_r8,     0.2894168E+03_r8,     0.2888797E+03_r8,     0.2883429E+03_r8,     0.2878062E+03_r8,&  
         0.2872696E+03_r8,     0.2867329E+03_r8,     0.2861965E+03_r8,     0.2856606E+03_r8,     0.2851249E+03_r8,&  
         0.2845885E+03_r8,     0.2840527E+03_r8,     0.2835166E+03_r8,     0.2829804E+03_r8,     0.2824445E+03_r8,&  
         0.2819080E+03_r8,     0.2813717E+03_r8,     0.2808351E+03_r8,     0.2802981E+03_r8,     0.2797610E+03_r8,&  
         0.2792232E+03_r8,     0.2786855E+03_r8,     0.2781473E+03_r8,     0.2776086E+03_r8,     0.2770696E+03_r8,&  
         0.2765296E+03_r8,     0.2759893E+03_r8,     0.2754485E+03_r8,     0.2749064E+03_r8,     0.2743640E+03_r8,&  
         0.2738210E+03_r8,     0.2732767E+03_r8,     0.2727317E+03_r8,     0.2721857E+03_r8,     0.2716383E+03_r8,&  
         0.2710903E+03_r8,     0.2705408E+03_r8,     0.2699901E+03_r8,     0.2694385E+03_r8,     0.2688850E+03_r8,&  
         0.2683303E+03_r8,     0.2677740E+03_r8,     0.2672160E+03_r8,     0.2666562E+03_r8,     0.2660949E+03_r8,&  
         0.2655316E+03_r8,     0.2649659E+03_r8,     0.2643984E+03_r8,     0.2638292E+03_r8,     0.2632568E+03_r8,&  
         0.2626828E+03_r8,     0.2621061E+03_r8,     0.2615264E+03_r8,     0.2609439E+03_r8,     0.2603591E+03_r8,&  
         0.2597711E+03_r8,     0.2591798E+03_r8,     0.2585852E+03_r8,     0.2579871E+03_r8,     0.2573857E+03_r8,&  
         0.2567800E+03_r8,     0.2561707E+03_r8,     0.2555570E+03_r8,     0.2549390E+03_r8,     0.2543171E+03_r8,&  
         0.2536897E+03_r8,     0.2530577E+03_r8,     0.2524207E+03_r8,     0.2517780E+03_r8,     0.2511297E+03_r8,&  
         0.2504755E+03_r8,     0.2498152E+03_r8,     0.2491485E+03_r8,     0.2484752E+03_r8,     0.2477946E+03_r8,&  
         0.2471065E+03_r8,     0.2464109E+03_r8,     0.2457074E+03_r8,     0.2449957E+03_r8,     0.2442749E+03_r8,&  
         0.2435449E+03_r8,     0.2428054E+03_r8,     0.2420560E+03_r8,     0.2412961E+03_r8,     0.2405253E+03_r8,&  
         0.2397434E+03_r8,     0.2389498E+03_r8,     0.2381437E+03_r8,     0.2373250E+03_r8,     0.2364930E+03_r8/)  

    psaditmk(1:150,  139)= (/ &
         0.3237905E+03_r8,     0.3231973E+03_r8,     0.3226050E+03_r8,     0.3220146E+03_r8,     0.3214252E+03_r8,&  
         0.3208380E+03_r8,     0.3202515E+03_r8,     0.3196671E+03_r8,     0.3190841E+03_r8,     0.3185019E+03_r8,&  
         0.3179211E+03_r8,     0.3173418E+03_r8,     0.3167639E+03_r8,     0.3161870E+03_r8,     0.3156116E+03_r8,&  
         0.3150375E+03_r8,     0.3144648E+03_r8,     0.3138929E+03_r8,     0.3133227E+03_r8,     0.3127534E+03_r8,&  
         0.3121852E+03_r8,     0.3116182E+03_r8,     0.3110524E+03_r8,     0.3104879E+03_r8,     0.3099247E+03_r8,&  
         0.3093619E+03_r8,     0.3088006E+03_r8,     0.3082404E+03_r8,     0.3076808E+03_r8,     0.3071225E+03_r8,&  
         0.3065654E+03_r8,     0.3060092E+03_r8,     0.3054542E+03_r8,     0.3048996E+03_r8,     0.3043466E+03_r8,&  
         0.3037939E+03_r8,     0.3032422E+03_r8,     0.3026915E+03_r8,     0.3021423E+03_r8,     0.3015932E+03_r8,&  
         0.3010448E+03_r8,     0.3004978E+03_r8,     0.2999514E+03_r8,     0.2994056E+03_r8,     0.2988607E+03_r8,&  
         0.2983164E+03_r8,     0.2977723E+03_r8,     0.2972299E+03_r8,     0.2966876E+03_r8,     0.2961462E+03_r8,&  
         0.2956051E+03_r8,     0.2950647E+03_r8,     0.2945247E+03_r8,     0.2939858E+03_r8,     0.2934471E+03_r8,&  
         0.2929089E+03_r8,     0.2923711E+03_r8,     0.2918336E+03_r8,     0.2912970E+03_r8,     0.2907606E+03_r8,&  
         0.2902248E+03_r8,     0.2896887E+03_r8,     0.2891535E+03_r8,     0.2886180E+03_r8,     0.2880834E+03_r8,&  
         0.2875490E+03_r8,     0.2870146E+03_r8,     0.2864804E+03_r8,     0.2859466E+03_r8,     0.2854124E+03_r8,&  
         0.2848786E+03_r8,     0.2843448E+03_r8,     0.2838116E+03_r8,     0.2832779E+03_r8,     0.2827438E+03_r8,&  
         0.2822104E+03_r8,     0.2816761E+03_r8,     0.2811422E+03_r8,     0.2806081E+03_r8,     0.2800734E+03_r8,&  
         0.2795388E+03_r8,     0.2790038E+03_r8,     0.2784680E+03_r8,     0.2779326E+03_r8,     0.2773963E+03_r8,&  
         0.2768594E+03_r8,     0.2763225E+03_r8,     0.2757845E+03_r8,     0.2752460E+03_r8,     0.2747072E+03_r8,&  
         0.2741672E+03_r8,     0.2736262E+03_r8,     0.2730854E+03_r8,     0.2725427E+03_r8,     0.2719998E+03_r8,&  
         0.2714554E+03_r8,     0.2709097E+03_r8,     0.2703632E+03_r8,     0.2698155E+03_r8,     0.2692662E+03_r8,&  
         0.2687161E+03_r8,     0.2681644E+03_r8,     0.2676111E+03_r8,     0.2670562E+03_r8,     0.2664998E+03_r8,&  
         0.2659413E+03_r8,     0.2653811E+03_r8,     0.2648193E+03_r8,     0.2642549E+03_r8,     0.2636887E+03_r8,&  
         0.2631204E+03_r8,     0.2625495E+03_r8,     0.2619764E+03_r8,     0.2614005E+03_r8,     0.2608219E+03_r8,&  
         0.2602405E+03_r8,     0.2596563E+03_r8,     0.2590689E+03_r8,     0.2584785E+03_r8,     0.2578845E+03_r8,&  
         0.2572871E+03_r8,     0.2566859E+03_r8,     0.2560809E+03_r8,     0.2554717E+03_r8,     0.2548587E+03_r8,&  
         0.2542408E+03_r8,     0.2536186E+03_r8,     0.2529915E+03_r8,     0.2523592E+03_r8,     0.2517219E+03_r8,&  
         0.2510792E+03_r8,     0.2504304E+03_r8,     0.2497757E+03_r8,     0.2491146E+03_r8,     0.2484471E+03_r8,&  
         0.2477730E+03_r8,     0.2470912E+03_r8,     0.2464021E+03_r8,     0.2457049E+03_r8,     0.2449999E+03_r8,&  
         0.2442861E+03_r8,     0.2435634E+03_r8,     0.2428315E+03_r8,     0.2420898E+03_r8,     0.2413378E+03_r8,&  
         0.2405753E+03_r8,     0.2398016E+03_r8,     0.2390165E+03_r8,     0.2382194E+03_r8,     0.2374098E+03_r8/)  

    psaditmk(1:150,  140)= (/ &
         0.3239862E+03_r8,     0.3233954E+03_r8,     0.3228044E+03_r8,     0.3222159E+03_r8,     0.3216272E+03_r8,&  
         0.3210397E+03_r8,     0.3204544E+03_r8,     0.3198711E+03_r8,     0.3192885E+03_r8,     0.3187070E+03_r8,&  
         0.3181269E+03_r8,     0.3175490E+03_r8,     0.3169716E+03_r8,     0.3163958E+03_r8,     0.3158212E+03_r8,&  
         0.3152476E+03_r8,     0.3146757E+03_r8,     0.3141045E+03_r8,     0.3135350E+03_r8,     0.3129670E+03_r8,&  
         0.3123997E+03_r8,     0.3118340E+03_r8,     0.3112688E+03_r8,     0.3107053E+03_r8,     0.3101425E+03_r8,&  
         0.3095811E+03_r8,     0.3090207E+03_r8,     0.3084614E+03_r8,     0.3079033E+03_r8,     0.3073459E+03_r8,&  
         0.3067896E+03_r8,     0.3062345E+03_r8,     0.3056802E+03_r8,     0.3051274E+03_r8,     0.3045753E+03_r8,&  
         0.3040238E+03_r8,     0.3034734E+03_r8,     0.3029238E+03_r8,     0.3023754E+03_r8,     0.3018278E+03_r8,&  
         0.3012812E+03_r8,     0.3007350E+03_r8,     0.3001895E+03_r8,     0.2996451E+03_r8,     0.2991014E+03_r8,&  
         0.2985587E+03_r8,     0.2980164E+03_r8,     0.2974745E+03_r8,     0.2969339E+03_r8,     0.2963938E+03_r8,&  
         0.2958549E+03_r8,     0.2953160E+03_r8,     0.2947773E+03_r8,     0.2942396E+03_r8,     0.2937025E+03_r8,&  
         0.2931659E+03_r8,     0.2926302E+03_r8,     0.2920946E+03_r8,     0.2915591E+03_r8,     0.2910244E+03_r8,&  
         0.2904902E+03_r8,     0.2899564E+03_r8,     0.2894231E+03_r8,     0.2888898E+03_r8,     0.2883567E+03_r8,&  
         0.2878239E+03_r8,     0.2872916E+03_r8,     0.2867596E+03_r8,     0.2862280E+03_r8,     0.2856961E+03_r8,&  
         0.2851645E+03_r8,     0.2846324E+03_r8,     0.2841013E+03_r8,     0.2835695E+03_r8,     0.2830385E+03_r8,&  
         0.2825069E+03_r8,     0.2819755E+03_r8,     0.2814440E+03_r8,     0.2809120E+03_r8,     0.2803803E+03_r8,&  
         0.2798483E+03_r8,     0.2793159E+03_r8,     0.2787833E+03_r8,     0.2782502E+03_r8,     0.2777169E+03_r8,&  
         0.2771835E+03_r8,     0.2766489E+03_r8,     0.2761143E+03_r8,     0.2755793E+03_r8,     0.2750432E+03_r8,&  
         0.2745066E+03_r8,     0.2739697E+03_r8,     0.2734313E+03_r8,     0.2728929E+03_r8,     0.2723533E+03_r8,&  
         0.2718122E+03_r8,     0.2712708E+03_r8,     0.2707283E+03_r8,     0.2701844E+03_r8,     0.2696394E+03_r8,&  
         0.2690936E+03_r8,     0.2685459E+03_r8,     0.2679973E+03_r8,     0.2674471E+03_r8,     0.2668952E+03_r8,&  
         0.2663417E+03_r8,     0.2657867E+03_r8,     0.2652297E+03_r8,     0.2646708E+03_r8,     0.2641101E+03_r8,&  
         0.2635473E+03_r8,     0.2629821E+03_r8,     0.2624146E+03_r8,     0.2618449E+03_r8,     0.2612727E+03_r8,&  
         0.2606978E+03_r8,     0.2601204E+03_r8,     0.2595400E+03_r8,     0.2589565E+03_r8,     0.2583699E+03_r8,&  
         0.2577798E+03_r8,     0.2571868E+03_r8,     0.2565897E+03_r8,     0.2559891E+03_r8,     0.2553844E+03_r8,&  
         0.2547755E+03_r8,     0.2541622E+03_r8,     0.2535452E+03_r8,     0.2529226E+03_r8,     0.2522955E+03_r8,&  
         0.2516632E+03_r8,     0.2510254E+03_r8,     0.2503823E+03_r8,     0.2497332E+03_r8,     0.2490775E+03_r8,&  
         0.2484160E+03_r8,     0.2477475E+03_r8,     0.2470721E+03_r8,     0.2463894E+03_r8,     0.2456987E+03_r8,&  
         0.2450005E+03_r8,     0.2442936E+03_r8,     0.2435781E+03_r8,     0.2428534E+03_r8,     0.2421189E+03_r8,&  
         0.2413750E+03_r8,     0.2406205E+03_r8,     0.2398550E+03_r8,     0.2390783E+03_r8,     0.2382901E+03_r8/)  

    psaditmk(1:150,  141)= (/ &
         0.8999999E+10_r8,     0.3235944E+03_r8,     0.3230036E+03_r8,     0.3224138E+03_r8,     0.3218258E+03_r8,&  
         0.3212402E+03_r8,     0.3206551E+03_r8,     0.3200728E+03_r8,     0.3194904E+03_r8,     0.3189099E+03_r8,&  
         0.3183304E+03_r8,     0.3177533E+03_r8,     0.3171769E+03_r8,     0.3166019E+03_r8,     0.3160278E+03_r8,&  
         0.3154553E+03_r8,     0.3148843E+03_r8,     0.3143137E+03_r8,     0.3137451E+03_r8,     0.3131777E+03_r8,&  
         0.3126116E+03_r8,     0.3120465E+03_r8,     0.3114825E+03_r8,     0.3109200E+03_r8,     0.3103581E+03_r8,&  
         0.3097976E+03_r8,     0.3092384E+03_r8,     0.3086801E+03_r8,     0.3081228E+03_r8,     0.3075665E+03_r8,&  
         0.3070114E+03_r8,     0.3064568E+03_r8,     0.3059038E+03_r8,     0.3053517E+03_r8,     0.3048009E+03_r8,&  
         0.3042506E+03_r8,     0.3037017E+03_r8,     0.3031530E+03_r8,     0.3026057E+03_r8,     0.3020590E+03_r8,&  
         0.3015135E+03_r8,     0.3009688E+03_r8,     0.3004250E+03_r8,     0.2998817E+03_r8,     0.2993394E+03_r8,&  
         0.2987978E+03_r8,     0.2982567E+03_r8,     0.2977169E+03_r8,     0.2971775E+03_r8,     0.2966385E+03_r8,&  
         0.2961002E+03_r8,     0.2955630E+03_r8,     0.2950264E+03_r8,     0.2944904E+03_r8,     0.2939549E+03_r8,&  
         0.2934197E+03_r8,     0.2928851E+03_r8,     0.2923514E+03_r8,     0.2918179E+03_r8,     0.2912847E+03_r8,&  
         0.2907519E+03_r8,     0.2902196E+03_r8,     0.2896882E+03_r8,     0.2891567E+03_r8,     0.2886256E+03_r8,&  
         0.2880951E+03_r8,     0.2875647E+03_r8,     0.2870341E+03_r8,     0.2865043E+03_r8,     0.2859745E+03_r8,&  
         0.2854452E+03_r8,     0.2849158E+03_r8,     0.2843865E+03_r8,     0.2838575E+03_r8,     0.2833277E+03_r8,&  
         0.2827988E+03_r8,     0.2822700E+03_r8,     0.2817405E+03_r8,     0.2812114E+03_r8,     0.2806819E+03_r8,&  
         0.2801522E+03_r8,     0.2796225E+03_r8,     0.2790925E+03_r8,     0.2785622E+03_r8,     0.2780318E+03_r8,&  
         0.2775007E+03_r8,     0.2769696E+03_r8,     0.2764379E+03_r8,     0.2759055E+03_r8,     0.2753729E+03_r8,&  
         0.2748393E+03_r8,     0.2743055E+03_r8,     0.2737708E+03_r8,     0.2732358E+03_r8,     0.2726994E+03_r8,&  
         0.2721623E+03_r8,     0.2716245E+03_r8,     0.2710855E+03_r8,     0.2705460E+03_r8,     0.2700049E+03_r8,&  
         0.2694627E+03_r8,     0.2689194E+03_r8,     0.2683751E+03_r8,     0.2678292E+03_r8,     0.2672820E+03_r8,&  
         0.2667331E+03_r8,     0.2661830E+03_r8,     0.2656306E+03_r8,     0.2650768E+03_r8,     0.2645213E+03_r8,&  
         0.2639634E+03_r8,     0.2634040E+03_r8,     0.2628425E+03_r8,     0.2622787E+03_r8,     0.2617124E+03_r8,&  
         0.2611437E+03_r8,     0.2605720E+03_r8,     0.2599983E+03_r8,     0.2594216E+03_r8,     0.2588419E+03_r8,&  
         0.2582596E+03_r8,     0.2576734E+03_r8,     0.2570842E+03_r8,     0.2564913E+03_r8,     0.2558949E+03_r8,&  
         0.2552945E+03_r8,     0.2546901E+03_r8,     0.2540820E+03_r8,     0.2534687E+03_r8,     0.2528513E+03_r8,&  
         0.2522289E+03_r8,     0.2516016E+03_r8,     0.2509692E+03_r8,     0.2503311E+03_r8,     0.2496874E+03_r8,&  
         0.2490377E+03_r8,     0.2483819E+03_r8,     0.2477190E+03_r8,     0.2470497E+03_r8,     0.2463729E+03_r8,&  
         0.2456890E+03_r8,     0.2449973E+03_r8,     0.2442971E+03_r8,     0.2435887E+03_r8,     0.2428711E+03_r8,&  
         0.2421444E+03_r8,     0.2414080E+03_r8,     0.2406610E+03_r8,     0.2399039E+03_r8,     0.2391356E+03_r8/)  

    psaditmk(1:150,  142)= (/ &
         0.8999999E+10_r8,     0.3237857E+03_r8,     0.3231994E+03_r8,     0.3226109E+03_r8,     0.3220236E+03_r8,&  
         0.3214381E+03_r8,     0.3208542E+03_r8,     0.3202716E+03_r8,     0.3196906E+03_r8,     0.3191102E+03_r8,&  
         0.3185320E+03_r8,     0.3179554E+03_r8,     0.3173794E+03_r8,     0.3168052E+03_r8,     0.3162325E+03_r8,&  
         0.3156605E+03_r8,     0.3150904E+03_r8,     0.3145209E+03_r8,     0.3139529E+03_r8,     0.3133861E+03_r8,&  
         0.3128206E+03_r8,     0.3122565E+03_r8,     0.3116935E+03_r8,     0.3111319E+03_r8,     0.3105708E+03_r8,&  
         0.3100115E+03_r8,     0.3094526E+03_r8,     0.3088955E+03_r8,     0.3083392E+03_r8,     0.3077841E+03_r8,&  
         0.3072302E+03_r8,     0.3066770E+03_r8,     0.3061249E+03_r8,     0.3055735E+03_r8,     0.3050235E+03_r8,&  
         0.3044745E+03_r8,     0.3039267E+03_r8,     0.3033792E+03_r8,     0.3028329E+03_r8,     0.3022880E+03_r8,&  
         0.3017431E+03_r8,     0.3011995E+03_r8,     0.3006567E+03_r8,     0.3001151E+03_r8,     0.2995738E+03_r8,&  
         0.2990334E+03_r8,     0.2984940E+03_r8,     0.2979548E+03_r8,     0.2974172E+03_r8,     0.2968799E+03_r8,&  
         0.2963434E+03_r8,     0.2958068E+03_r8,     0.2952715E+03_r8,     0.2947367E+03_r8,     0.2942032E+03_r8,&  
         0.2936696E+03_r8,     0.2931364E+03_r8,     0.2926039E+03_r8,     0.2920719E+03_r8,     0.2915409E+03_r8,&  
         0.2910099E+03_r8,     0.2904796E+03_r8,     0.2899493E+03_r8,     0.2894198E+03_r8,     0.2888906E+03_r8,&  
         0.2883616E+03_r8,     0.2878332E+03_r8,     0.2873049E+03_r8,     0.2867766E+03_r8,     0.2862487E+03_r8,&  
         0.2857211E+03_r8,     0.2851938E+03_r8,     0.2846667E+03_r8,     0.2841396E+03_r8,     0.2836129E+03_r8,&  
         0.2830858E+03_r8,     0.2825588E+03_r8,     0.2820320E+03_r8,     0.2815049E+03_r8,     0.2809780E+03_r8,&  
         0.2804512E+03_r8,     0.2799237E+03_r8,     0.2793964E+03_r8,     0.2788689E+03_r8,     0.2783406E+03_r8,&  
         0.2778128E+03_r8,     0.2772842E+03_r8,     0.2767552E+03_r8,     0.2762260E+03_r8,     0.2756962E+03_r8,&  
         0.2751658E+03_r8,     0.2746351E+03_r8,     0.2741036E+03_r8,     0.2735716E+03_r8,     0.2730390E+03_r8,&  
         0.2725052E+03_r8,     0.2719709E+03_r8,     0.2714358E+03_r8,     0.2708995E+03_r8,     0.2703623E+03_r8,&  
         0.2698243E+03_r8,     0.2692850E+03_r8,     0.2687446E+03_r8,     0.2682027E+03_r8,     0.2676599E+03_r8,&  
         0.2671156E+03_r8,     0.2665698E+03_r8,     0.2660225E+03_r8,     0.2654737E+03_r8,     0.2649226E+03_r8,&  
         0.2643703E+03_r8,     0.2638162E+03_r8,     0.2632598E+03_r8,     0.2627011E+03_r8,     0.2621407E+03_r8,&  
         0.2615779E+03_r8,     0.2610129E+03_r8,     0.2604452E+03_r8,     0.2598749E+03_r8,     0.2593017E+03_r8,&  
         0.2587260E+03_r8,     0.2581470E+03_r8,     0.2575647E+03_r8,     0.2569798E+03_r8,     0.2563909E+03_r8,&  
         0.2557986E+03_r8,     0.2552029E+03_r8,     0.2546027E+03_r8,     0.2539985E+03_r8,     0.2533903E+03_r8,&  
         0.2527774E+03_r8,     0.2521600E+03_r8,     0.2515376E+03_r8,     0.2509100E+03_r8,     0.2502775E+03_r8,&  
         0.2496389E+03_r8,     0.2489945E+03_r8,     0.2483443E+03_r8,     0.2476872E+03_r8,     0.2470239E+03_r8,&  
         0.2463536E+03_r8,     0.2456756E+03_r8,     0.2449905E+03_r8,     0.2442969E+03_r8,     0.2435954E+03_r8,&  
         0.2428849E+03_r8,     0.2421654E+03_r8,     0.2414364E+03_r8,     0.2406974E+03_r8,     0.2399481E+03_r8/)  

    psaditmk(1:150,  143)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.3233897E+03_r8,     0.3228040E+03_r8,     0.3222179E+03_r8,&  
         0.3216338E+03_r8,     0.3210503E+03_r8,     0.3204682E+03_r8,     0.3198882E+03_r8,     0.3193093E+03_r8,&  
         0.3187312E+03_r8,     0.3181553E+03_r8,     0.3175802E+03_r8,     0.3170066E+03_r8,     0.3164343E+03_r8,&  
         0.3158631E+03_r8,     0.3152934E+03_r8,     0.3147251E+03_r8,     0.3141585E+03_r8,     0.3135926E+03_r8,&  
         0.3130279E+03_r8,     0.3124643E+03_r8,     0.3119025E+03_r8,     0.3113410E+03_r8,     0.3107811E+03_r8,&  
         0.3102225E+03_r8,     0.3096650E+03_r8,     0.3091086E+03_r8,     0.3085531E+03_r8,     0.3079992E+03_r8,&  
         0.3074461E+03_r8,     0.3068938E+03_r8,     0.3063430E+03_r8,     0.3057928E+03_r8,     0.3052437E+03_r8,&  
         0.3046957E+03_r8,     0.3041487E+03_r8,     0.3036024E+03_r8,     0.3030576E+03_r8,     0.3025135E+03_r8,&  
         0.3019702E+03_r8,     0.3014271E+03_r8,     0.3008857E+03_r8,     0.3003451E+03_r8,     0.2998052E+03_r8,&  
         0.2992663E+03_r8,     0.2987280E+03_r8,     0.2981906E+03_r8,     0.2976534E+03_r8,     0.2971175E+03_r8,&  
         0.2965822E+03_r8,     0.2960478E+03_r8,     0.2955135E+03_r8,     0.2949803E+03_r8,     0.2944477E+03_r8,&  
         0.2939160E+03_r8,     0.2933844E+03_r8,     0.2928534E+03_r8,     0.2923228E+03_r8,     0.2917934E+03_r8,&  
         0.2912639E+03_r8,     0.2907352E+03_r8,     0.2902066E+03_r8,     0.2896790E+03_r8,     0.2891515E+03_r8,&  
         0.2886239E+03_r8,     0.2880973E+03_r8,     0.2875709E+03_r8,     0.2870447E+03_r8,     0.2865190E+03_r8,&  
         0.2859935E+03_r8,     0.2854676E+03_r8,     0.2849423E+03_r8,     0.2844177E+03_r8,     0.2838929E+03_r8,&  
         0.2833680E+03_r8,     0.2828434E+03_r8,     0.2823185E+03_r8,     0.2817939E+03_r8,     0.2812695E+03_r8,&  
         0.2807444E+03_r8,     0.2802199E+03_r8,     0.2796946E+03_r8,     0.2791696E+03_r8,     0.2786444E+03_r8,&  
         0.2781187E+03_r8,     0.2775929E+03_r8,     0.2770670E+03_r8,     0.2765404E+03_r8,     0.2760136E+03_r8,&  
         0.2754865E+03_r8,     0.2749585E+03_r8,     0.2744300E+03_r8,     0.2739014E+03_r8,     0.2733717E+03_r8,&  
         0.2728413E+03_r8,     0.2723105E+03_r8,     0.2717789E+03_r8,     0.2712459E+03_r8,     0.2707126E+03_r8,&  
         0.2701783E+03_r8,     0.2696426E+03_r8,     0.2691061E+03_r8,     0.2685688E+03_r8,     0.2680297E+03_r8,&  
         0.2674898E+03_r8,     0.2669482E+03_r8,     0.2664055E+03_r8,     0.2658610E+03_r8,     0.2653151E+03_r8,&  
         0.2647676E+03_r8,     0.2642182E+03_r8,     0.2636667E+03_r8,     0.2631140E+03_r8,     0.2625584E+03_r8,&  
         0.2620013E+03_r8,     0.2614420E+03_r8,     0.2608803E+03_r8,     0.2603159E+03_r8,     0.2597494E+03_r8,&  
         0.2591798E+03_r8,     0.2586078E+03_r8,     0.2580328E+03_r8,     0.2574547E+03_r8,     0.2568734E+03_r8,&  
         0.2562887E+03_r8,     0.2557005E+03_r8,     0.2551085E+03_r8,     0.2545127E+03_r8,     0.2539133E+03_r8,&  
         0.2533094E+03_r8,     0.2527012E+03_r8,     0.2520885E+03_r8,     0.2514708E+03_r8,     0.2508484E+03_r8,&  
         0.2502205E+03_r8,     0.2495875E+03_r8,     0.2489485E+03_r8,     0.2483036E+03_r8,     0.2476526E+03_r8,&  
         0.2469949E+03_r8,     0.2463305E+03_r8,     0.2456588E+03_r8,     0.2449800E+03_r8,     0.2442931E+03_r8,&  
         0.2435982E+03_r8,     0.2428949E+03_r8,     0.2421823E+03_r8,     0.2414605E+03_r8,     0.2407294E+03_r8/)  

    psaditmk(1:150,  144)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.3230028E+03_r8,     0.3224122E+03_r8,&  
         0.3218277E+03_r8,     0.3212447E+03_r8,     0.3206636E+03_r8,     0.3200836E+03_r8,     0.3195050E+03_r8,&  
         0.3189283E+03_r8,     0.3183529E+03_r8,     0.3177787E+03_r8,     0.3172054E+03_r8,     0.3166344E+03_r8,&  
         0.3160642E+03_r8,     0.3154950E+03_r8,     0.3149272E+03_r8,     0.3143609E+03_r8,     0.3137963E+03_r8,&  
         0.3132323E+03_r8,     0.3126699E+03_r8,     0.3121087E+03_r8,     0.3115485E+03_r8,     0.3109892E+03_r8,&  
         0.3104314E+03_r8,     0.3098745E+03_r8,     0.3093192E+03_r8,     0.3087646E+03_r8,     0.3082115E+03_r8,&  
         0.3076593E+03_r8,     0.3071082E+03_r8,     0.3065584E+03_r8,     0.3060093E+03_r8,     0.3054613E+03_r8,&  
         0.3049145E+03_r8,     0.3043682E+03_r8,     0.3038229E+03_r8,     0.3032789E+03_r8,     0.3027357E+03_r8,&  
         0.3021938E+03_r8,     0.3016526E+03_r8,     0.3011121E+03_r8,     0.3005720E+03_r8,     0.3000333E+03_r8,&  
         0.2994954E+03_r8,     0.2989587E+03_r8,     0.2984222E+03_r8,     0.2978870E+03_r8,     0.2973521E+03_r8,&  
         0.2968178E+03_r8,     0.2962849E+03_r8,     0.2957524E+03_r8,     0.2952206E+03_r8,     0.2946889E+03_r8,&  
         0.2941585E+03_r8,     0.2936282E+03_r8,     0.2930994E+03_r8,     0.2925704E+03_r8,     0.2920423E+03_r8,&  
         0.2915141E+03_r8,     0.2909870E+03_r8,     0.2904604E+03_r8,     0.2899343E+03_r8,     0.2894084E+03_r8,&  
         0.2888829E+03_r8,     0.2883580E+03_r8,     0.2878331E+03_r8,     0.2873086E+03_r8,     0.2867843E+03_r8,&  
         0.2862611E+03_r8,     0.2857375E+03_r8,     0.2852144E+03_r8,     0.2846913E+03_r8,     0.2841681E+03_r8,&  
         0.2836456E+03_r8,     0.2831229E+03_r8,     0.2826006E+03_r8,     0.2820778E+03_r8,     0.2815556E+03_r8,&  
         0.2810333E+03_r8,     0.2805105E+03_r8,     0.2799879E+03_r8,     0.2794655E+03_r8,     0.2789425E+03_r8,&  
         0.2784196E+03_r8,     0.2778965E+03_r8,     0.2773729E+03_r8,     0.2768492E+03_r8,     0.2763250E+03_r8,&  
         0.2758006E+03_r8,     0.2752756E+03_r8,     0.2747506E+03_r8,     0.2742243E+03_r8,     0.2736980E+03_r8,&  
         0.2731712E+03_r8,     0.2726433E+03_r8,     0.2721147E+03_r8,     0.2715858E+03_r8,     0.2710557E+03_r8,&  
         0.2705247E+03_r8,     0.2699934E+03_r8,     0.2694605E+03_r8,     0.2689266E+03_r8,     0.2683918E+03_r8,&  
         0.2678558E+03_r8,     0.2673187E+03_r8,     0.2667798E+03_r8,     0.2662396E+03_r8,     0.2656986E+03_r8,&  
         0.2651553E+03_r8,     0.2646107E+03_r8,     0.2640648E+03_r8,     0.2635168E+03_r8,     0.2629666E+03_r8,&  
         0.2624149E+03_r8,     0.2618612E+03_r8,     0.2613049E+03_r8,     0.2607465E+03_r8,     0.2601858E+03_r8,&  
         0.2596227E+03_r8,     0.2590568E+03_r8,     0.2584883E+03_r8,     0.2579167E+03_r8,     0.2573423E+03_r8,&  
         0.2567653E+03_r8,     0.2561843E+03_r8,     0.2556001E+03_r8,     0.2550124E+03_r8,     0.2544214E+03_r8,&  
         0.2538258E+03_r8,     0.2532262E+03_r8,     0.2526227E+03_r8,     0.2520145E+03_r8,     0.2514017E+03_r8,&  
         0.2507840E+03_r8,     0.2501615E+03_r8,     0.2495332E+03_r8,     0.2488999E+03_r8,     0.2482601E+03_r8,&  
         0.2476148E+03_r8,     0.2469631E+03_r8,     0.2463042E+03_r8,     0.2456385E+03_r8,     0.2449661E+03_r8,&  
         0.2442856E+03_r8,     0.2435974E+03_r8,     0.2429006E+03_r8,     0.2421952E+03_r8,     0.2414808E+03_r8/)  

    psaditmk(1:150,  145)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.3226006E+03_r8,&  
         0.3220179E+03_r8,     0.3214376E+03_r8,     0.3208557E+03_r8,     0.3202767E+03_r8,     0.3196991E+03_r8,&  
         0.3191229E+03_r8,     0.3185481E+03_r8,     0.3179749E+03_r8,     0.3174029E+03_r8,     0.3168318E+03_r8,&  
         0.3162621E+03_r8,     0.3156940E+03_r8,     0.3151273E+03_r8,     0.3145618E+03_r8,     0.3139975E+03_r8,&  
         0.3134343E+03_r8,     0.3128728E+03_r8,     0.3123120E+03_r8,     0.3117529E+03_r8,     0.3111947E+03_r8,&  
         0.3106381E+03_r8,     0.3100822E+03_r8,     0.3095270E+03_r8,     0.3089738E+03_r8,     0.3084216E+03_r8,&  
         0.3078701E+03_r8,     0.3073200E+03_r8,     0.3067706E+03_r8,     0.3062231E+03_r8,     0.3056758E+03_r8,&  
         0.3051302E+03_r8,     0.3045849E+03_r8,     0.3040409E+03_r8,     0.3034977E+03_r8,     0.3029555E+03_r8,&  
         0.3024144E+03_r8,     0.3018744E+03_r8,     0.3013350E+03_r8,     0.3007968E+03_r8,     0.3002589E+03_r8,&  
         0.2997223E+03_r8,     0.2991865E+03_r8,     0.2986513E+03_r8,     0.2981169E+03_r8,     0.2975835E+03_r8,&  
         0.2970509E+03_r8,     0.2965189E+03_r8,     0.2959874E+03_r8,     0.2954569E+03_r8,     0.2949273E+03_r8,&  
         0.2943982E+03_r8,     0.2938691E+03_r8,     0.2933411E+03_r8,     0.2928138E+03_r8,     0.2922874E+03_r8,&  
         0.2917611E+03_r8,     0.2912354E+03_r8,     0.2907102E+03_r8,     0.2901853E+03_r8,     0.2896611E+03_r8,&  
         0.2891372E+03_r8,     0.2886144E+03_r8,     0.2880912E+03_r8,     0.2875688E+03_r8,     0.2870462E+03_r8,&  
         0.2865241E+03_r8,     0.2860025E+03_r8,     0.2854813E+03_r8,     0.2849604E+03_r8,     0.2844395E+03_r8,&  
         0.2839189E+03_r8,     0.2833982E+03_r8,     0.2828773E+03_r8,     0.2823575E+03_r8,     0.2818372E+03_r8,&  
         0.2813168E+03_r8,     0.2807967E+03_r8,     0.2802765E+03_r8,     0.2797560E+03_r8,     0.2792357E+03_r8,&  
         0.2787154E+03_r8,     0.2781944E+03_r8,     0.2776737E+03_r8,     0.2771523E+03_r8,     0.2766307E+03_r8,&  
         0.2761093E+03_r8,     0.2755871E+03_r8,     0.2750645E+03_r8,     0.2745415E+03_r8,     0.2740181E+03_r8,&  
         0.2734941E+03_r8,     0.2729694E+03_r8,     0.2724445E+03_r8,     0.2719185E+03_r8,     0.2713918E+03_r8,&  
         0.2708646E+03_r8,     0.2703363E+03_r8,     0.2698069E+03_r8,     0.2692770E+03_r8,     0.2687460E+03_r8,&  
         0.2682139E+03_r8,     0.2676806E+03_r8,     0.2671465E+03_r8,     0.2666106E+03_r8,     0.2660732E+03_r8,&  
         0.2655347E+03_r8,     0.2649950E+03_r8,     0.2644532E+03_r8,     0.2639097E+03_r8,     0.2633649E+03_r8,&  
         0.2628183E+03_r8,     0.2622693E+03_r8,     0.2617185E+03_r8,     0.2611659E+03_r8,     0.2606109E+03_r8,&  
         0.2600537E+03_r8,     0.2594940E+03_r8,     0.2589317E+03_r8,     0.2583670E+03_r8,     0.2577991E+03_r8,&  
         0.2572287E+03_r8,     0.2566551E+03_r8,     0.2560782E+03_r8,     0.2554980E+03_r8,     0.2549141E+03_r8,&  
         0.2543269E+03_r8,     0.2537360E+03_r8,     0.2531409E+03_r8,     0.2525417E+03_r8,     0.2519383E+03_r8,&  
         0.2513298E+03_r8,     0.2507173E+03_r8,     0.2500993E+03_r8,     0.2494767E+03_r8,     0.2488478E+03_r8,&  
         0.2482139E+03_r8,     0.2475738E+03_r8,     0.2469276E+03_r8,     0.2462746E+03_r8,     0.2456152E+03_r8,&  
         0.2449487E+03_r8,     0.2442744E+03_r8,     0.2435925E+03_r8,     0.2429027E+03_r8,     0.2422042E+03_r8/)  

    psaditmk(1:150,  146)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,&  
         0.3222098E+03_r8,     0.3216255E+03_r8,     0.3210467E+03_r8,     0.3204680E+03_r8,     0.3198916E+03_r8,&  
         0.3193160E+03_r8,     0.3187419E+03_r8,     0.3181684E+03_r8,     0.3175975E+03_r8,     0.3170271E+03_r8,&  
         0.3164586E+03_r8,     0.3158910E+03_r8,     0.3153254E+03_r8,     0.3147602E+03_r8,     0.3141969E+03_r8,&  
         0.3136345E+03_r8,     0.3130734E+03_r8,     0.3125137E+03_r8,     0.3119552E+03_r8,     0.3113979E+03_r8,&  
         0.3108418E+03_r8,     0.3102867E+03_r8,     0.3097330E+03_r8,     0.3091802E+03_r8,     0.3086289E+03_r8,&  
         0.3080784E+03_r8,     0.3075292E+03_r8,     0.3069812E+03_r8,     0.3064339E+03_r8,     0.3058876E+03_r8,&  
         0.3053432E+03_r8,     0.3047990E+03_r8,     0.3042560E+03_r8,     0.3037139E+03_r8,     0.3031725E+03_r8,&  
         0.3026324E+03_r8,     0.3020932E+03_r8,     0.3015551E+03_r8,     0.3010179E+03_r8,     0.3004817E+03_r8,&  
         0.2999460E+03_r8,     0.2994110E+03_r8,     0.2988772E+03_r8,     0.2983440E+03_r8,     0.2978121E+03_r8,&  
         0.2972803E+03_r8,     0.2967499E+03_r8,     0.2962194E+03_r8,     0.2956901E+03_r8,     0.2951619E+03_r8,&  
         0.2946341E+03_r8,     0.2941069E+03_r8,     0.2935804E+03_r8,     0.2930541E+03_r8,     0.2925289E+03_r8,&  
         0.2920042E+03_r8,     0.2914802E+03_r8,     0.2909564E+03_r8,     0.2904332E+03_r8,     0.2899105E+03_r8,&  
         0.2893884E+03_r8,     0.2888665E+03_r8,     0.2883455E+03_r8,     0.2878247E+03_r8,     0.2873041E+03_r8,&  
         0.2867840E+03_r8,     0.2862639E+03_r8,     0.2857441E+03_r8,     0.2852249E+03_r8,     0.2847059E+03_r8,&  
         0.2841876E+03_r8,     0.2836688E+03_r8,     0.2831507E+03_r8,     0.2826320E+03_r8,     0.2821139E+03_r8,&  
         0.2815959E+03_r8,     0.2810776E+03_r8,     0.2805594E+03_r8,     0.2800417E+03_r8,     0.2795236E+03_r8,&  
         0.2790051E+03_r8,     0.2784874E+03_r8,     0.2779686E+03_r8,     0.2774499E+03_r8,     0.2769313E+03_r8,&  
         0.2764122E+03_r8,     0.2758926E+03_r8,     0.2753730E+03_r8,     0.2748528E+03_r8,     0.2743322E+03_r8,&  
         0.2738112E+03_r8,     0.2732897E+03_r8,     0.2727676E+03_r8,     0.2722447E+03_r8,     0.2717218E+03_r8,&  
         0.2711973E+03_r8,     0.2706725E+03_r8,     0.2701471E+03_r8,     0.2696204E+03_r8,     0.2690927E+03_r8,&  
         0.2685646E+03_r8,     0.2680352E+03_r8,     0.2675049E+03_r8,     0.2669727E+03_r8,     0.2664397E+03_r8,&  
         0.2659059E+03_r8,     0.2653701E+03_r8,     0.2648329E+03_r8,     0.2642942E+03_r8,     0.2637541E+03_r8,&  
         0.2632123E+03_r8,     0.2626685E+03_r8,     0.2621228E+03_r8,     0.2615754E+03_r8,     0.2610259E+03_r8,&  
         0.2604742E+03_r8,     0.2599202E+03_r8,     0.2593639E+03_r8,     0.2588051E+03_r8,     0.2582438E+03_r8,&  
         0.2576797E+03_r8,     0.2571129E+03_r8,     0.2565429E+03_r8,     0.2559700E+03_r8,     0.2553937E+03_r8,&  
         0.2548140E+03_r8,     0.2542309E+03_r8,     0.2536441E+03_r8,     0.2530535E+03_r8,     0.2524585E+03_r8,&  
         0.2518596E+03_r8,     0.2512561E+03_r8,     0.2506479E+03_r8,     0.2500353E+03_r8,     0.2494167E+03_r8,&  
         0.2487936E+03_r8,     0.2481647E+03_r8,     0.2475300E+03_r8,     0.2468892E+03_r8,     0.2462420E+03_r8,&  
         0.2455885E+03_r8,     0.2449279E+03_r8,     0.2442599E+03_r8,     0.2435844E+03_r8,     0.2429011E+03_r8/)  

    psaditmk(1:150,  147)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,&  
         0.3223920E+03_r8,     0.3218184E+03_r8,     0.3212359E+03_r8,     0.3206581E+03_r8,     0.3200815E+03_r8,&  
         0.3195058E+03_r8,     0.3189323E+03_r8,     0.3183611E+03_r8,     0.3177897E+03_r8,     0.3172209E+03_r8,&  
         0.3166523E+03_r8,     0.3160858E+03_r8,     0.3155206E+03_r8,     0.3149567E+03_r8,     0.3143935E+03_r8,&  
         0.3138326E+03_r8,     0.3132721E+03_r8,     0.3127131E+03_r8,     0.3121553E+03_r8,     0.3115984E+03_r8,&  
         0.3110437E+03_r8,     0.3104893E+03_r8,     0.3099363E+03_r8,     0.3093847E+03_r8,     0.3088338E+03_r8,&  
         0.3082841E+03_r8,     0.3077359E+03_r8,     0.3071886E+03_r8,     0.3066425E+03_r8,     0.3060970E+03_r8,&  
         0.3055532E+03_r8,     0.3050099E+03_r8,     0.3044685E+03_r8,     0.3039274E+03_r8,     0.3033873E+03_r8,&  
         0.3028481E+03_r8,     0.3023096E+03_r8,     0.3017724E+03_r8,     0.3012364E+03_r8,     0.3007011E+03_r8,&  
         0.3001665E+03_r8,     0.2996330E+03_r8,     0.2991002E+03_r8,     0.2985684E+03_r8,     0.2980371E+03_r8,&  
         0.2975068E+03_r8,     0.2969777E+03_r8,     0.2964486E+03_r8,     0.2959206E+03_r8,     0.2953932E+03_r8,&  
         0.2948666E+03_r8,     0.2943412E+03_r8,     0.2938159E+03_r8,     0.2932915E+03_r8,     0.2927674E+03_r8,&  
         0.2922436E+03_r8,     0.2917210E+03_r8,     0.2911990E+03_r8,     0.2906777E+03_r8,     0.2901565E+03_r8,&  
         0.2896359E+03_r8,     0.2891154E+03_r8,     0.2885959E+03_r8,     0.2880764E+03_r8,     0.2875577E+03_r8,&  
         0.2870395E+03_r8,     0.2865211E+03_r8,     0.2860036E+03_r8,     0.2854859E+03_r8,     0.2849688E+03_r8,&  
         0.2844519E+03_r8,     0.2839352E+03_r8,     0.2834187E+03_r8,     0.2829025E+03_r8,     0.2823863E+03_r8,&  
         0.2818700E+03_r8,     0.2813543E+03_r8,     0.2808386E+03_r8,     0.2803225E+03_r8,     0.2798068E+03_r8,&  
         0.2792908E+03_r8,     0.2787749E+03_r8,     0.2782591E+03_r8,     0.2777428E+03_r8,     0.2772259E+03_r8,&  
         0.2767097E+03_r8,     0.2761928E+03_r8,     0.2756756E+03_r8,     0.2751584E+03_r8,     0.2746409E+03_r8,&  
         0.2741223E+03_r8,     0.2736036E+03_r8,     0.2730847E+03_r8,     0.2725649E+03_r8,     0.2720446E+03_r8,&  
         0.2715237E+03_r8,     0.2710022E+03_r8,     0.2704799E+03_r8,     0.2699566E+03_r8,     0.2694329E+03_r8,&  
         0.2689081E+03_r8,     0.2683822E+03_r8,     0.2678554E+03_r8,     0.2673278E+03_r8,     0.2667988E+03_r8,&  
         0.2662684E+03_r8,     0.2657370E+03_r8,     0.2652043E+03_r8,     0.2646701E+03_r8,     0.2641343E+03_r8,&  
         0.2635968E+03_r8,     0.2630582E+03_r8,     0.2625177E+03_r8,     0.2619752E+03_r8,     0.2614311E+03_r8,&  
         0.2608845E+03_r8,     0.2603360E+03_r8,     0.2597857E+03_r8,     0.2592325E+03_r8,     0.2586771E+03_r8,&  
         0.2581192E+03_r8,     0.2575585E+03_r8,     0.2569956E+03_r8,     0.2564294E+03_r8,     0.2558603E+03_r8,&  
         0.2552879E+03_r8,     0.2547122E+03_r8,     0.2541329E+03_r8,     0.2535500E+03_r8,     0.2529638E+03_r8,&  
         0.2523733E+03_r8,     0.2517785E+03_r8,     0.2511796E+03_r8,     0.2505762E+03_r8,     0.2499677E+03_r8,&  
         0.2493548E+03_r8,     0.2487364E+03_r8,     0.2481127E+03_r8,     0.2474834E+03_r8,     0.2468478E+03_r8,&  
         0.2462065E+03_r8,     0.2455585E+03_r8,     0.2449039E+03_r8,     0.2442420E+03_r8,     0.2435727E+03_r8/)  

    psaditmk(1:150,  148)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,&  
         0.8999999E+10_r8,     0.3219957E+03_r8,     0.3214211E+03_r8,     0.3208443E+03_r8,     0.3202687E+03_r8,&  
         0.3196947E+03_r8,     0.3191218E+03_r8,     0.3185506E+03_r8,     0.3179805E+03_r8,     0.3174117E+03_r8,&  
         0.3168443E+03_r8,     0.3162787E+03_r8,     0.3157137E+03_r8,     0.3151507E+03_r8,     0.3145882E+03_r8,&  
         0.3140279E+03_r8,     0.3134680E+03_r8,     0.3129100E+03_r8,     0.3123530E+03_r8,     0.3117970E+03_r8,&  
         0.3112428E+03_r8,     0.3106894E+03_r8,     0.3101374E+03_r8,     0.3095866E+03_r8,     0.3090367E+03_r8,&  
         0.3084883E+03_r8,     0.3079402E+03_r8,     0.3073940E+03_r8,     0.3068483E+03_r8,     0.3063043E+03_r8,&  
         0.3057611E+03_r8,     0.3052188E+03_r8,     0.3046780E+03_r8,     0.3041378E+03_r8,     0.3035990E+03_r8,&  
         0.3030611E+03_r8,     0.3025239E+03_r8,     0.3019873E+03_r8,     0.3014520E+03_r8,     0.3009177E+03_r8,&  
         0.3003843E+03_r8,     0.2998524E+03_r8,     0.2993205E+03_r8,     0.2987896E+03_r8,     0.2982596E+03_r8,&  
         0.2977302E+03_r8,     0.2972022E+03_r8,     0.2966748E+03_r8,     0.2961481E+03_r8,     0.2956219E+03_r8,&  
         0.2950964E+03_r8,     0.2945718E+03_r8,     0.2940479E+03_r8,     0.2935250E+03_r8,     0.2930023E+03_r8,&  
         0.2924807E+03_r8,     0.2919591E+03_r8,     0.2914380E+03_r8,     0.2909180E+03_r8,     0.2903985E+03_r8,&  
         0.2898795E+03_r8,     0.2893608E+03_r8,     0.2888427E+03_r8,     0.2883250E+03_r8,     0.2878075E+03_r8,&  
         0.2872910E+03_r8,     0.2867743E+03_r8,     0.2862586E+03_r8,     0.2857429E+03_r8,     0.2852274E+03_r8,&  
         0.2847122E+03_r8,     0.2841974E+03_r8,     0.2836828E+03_r8,     0.2831684E+03_r8,     0.2826542E+03_r8,&  
         0.2821405E+03_r8,     0.2816261E+03_r8,     0.2811125E+03_r8,     0.2805988E+03_r8,     0.2800851E+03_r8,&  
         0.2795711E+03_r8,     0.2790578E+03_r8,     0.2785440E+03_r8,     0.2780298E+03_r8,     0.2775163E+03_r8,&  
         0.2770021E+03_r8,     0.2764874E+03_r8,     0.2759734E+03_r8,     0.2754585E+03_r8,     0.2749431E+03_r8,&  
         0.2744278E+03_r8,     0.2739119E+03_r8,     0.2733956E+03_r8,     0.2728787E+03_r8,     0.2723618E+03_r8,&  
         0.2718438E+03_r8,     0.2713252E+03_r8,     0.2708062E+03_r8,     0.2702864E+03_r8,     0.2697661E+03_r8,&  
         0.2692446E+03_r8,     0.2687223E+03_r8,     0.2681992E+03_r8,     0.2676752E+03_r8,     0.2671496E+03_r8,&  
         0.2666234E+03_r8,     0.2660962E+03_r8,     0.2655674E+03_r8,     0.2650373E+03_r8,     0.2645060E+03_r8,&  
         0.2639732E+03_r8,     0.2634391E+03_r8,     0.2629030E+03_r8,     0.2623654E+03_r8,     0.2618262E+03_r8,&  
         0.2612849E+03_r8,     0.2607416E+03_r8,     0.2601963E+03_r8,     0.2596490E+03_r8,     0.2590991E+03_r8,&  
         0.2585475E+03_r8,     0.2579931E+03_r8,     0.2574361E+03_r8,     0.2568767E+03_r8,     0.2563139E+03_r8,&  
         0.2557484E+03_r8,     0.2551798E+03_r8,     0.2546082E+03_r8,     0.2540330E+03_r8,     0.2534542E+03_r8,&  
         0.2528719E+03_r8,     0.2522857E+03_r8,     0.2516955E+03_r8,     0.2511006E+03_r8,     0.2505020E+03_r8,&  
         0.2498985E+03_r8,     0.2492900E+03_r8,     0.2486766E+03_r8,     0.2480581E+03_r8,     0.2474336E+03_r8,&  
         0.2468038E+03_r8,     0.2461678E+03_r8,     0.2455257E+03_r8,     0.2448767E+03_r8,     0.2442207E+03_r8/)  

    psaditmk(1:150,  149)= (/ &
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,&  
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.3216098E+03_r8,     0.3210316E+03_r8,     0.3204551E+03_r8,&  
         0.3198815E+03_r8,     0.3193090E+03_r8,     0.3187384E+03_r8,     0.3181698E+03_r8,     0.3176014E+03_r8,&  
         0.3170342E+03_r8,     0.3164689E+03_r8,     0.3159052E+03_r8,     0.3153423E+03_r8,     0.3147809E+03_r8,&  
         0.3142211E+03_r8,     0.3136623E+03_r8,     0.3131050E+03_r8,     0.3125486E+03_r8,     0.3119938E+03_r8,&  
         0.3114399E+03_r8,     0.3108871E+03_r8,     0.3103359E+03_r8,     0.3097860E+03_r8,     0.3092370E+03_r8,&  
         0.3086890E+03_r8,     0.3081424E+03_r8,     0.3075968E+03_r8,     0.3070524E+03_r8,     0.3065092E+03_r8,&  
         0.3059666E+03_r8,     0.3054252E+03_r8,     0.3048851E+03_r8,     0.3043462E+03_r8,     0.3038076E+03_r8,&  
         0.3032712E+03_r8,     0.3027348E+03_r8,     0.3021996E+03_r8,     0.3016654E+03_r8,     0.3011318E+03_r8,&  
         0.3005995E+03_r8,     0.3000680E+03_r8,     0.2995380E+03_r8,     0.2990082E+03_r8,     0.2984794E+03_r8,&  
         0.2979513E+03_r8,     0.2974238E+03_r8,     0.2968976E+03_r8,     0.2963721E+03_r8,     0.2958473E+03_r8,&  
         0.2953235E+03_r8,     0.2947998E+03_r8,     0.2942771E+03_r8,     0.2937551E+03_r8,     0.2932341E+03_r8,&  
         0.2927137E+03_r8,     0.2921937E+03_r8,     0.2916741E+03_r8,     0.2911553E+03_r8,     0.2906373E+03_r8,&  
         0.2901192E+03_r8,     0.2896024E+03_r8,     0.2890860E+03_r8,     0.2885699E+03_r8,     0.2880543E+03_r8,&  
         0.2875389E+03_r8,     0.2870239E+03_r8,     0.2865094E+03_r8,     0.2859953E+03_r8,     0.2854819E+03_r8,&  
         0.2849686E+03_r8,     0.2844557E+03_r8,     0.2839431E+03_r8,     0.2834303E+03_r8,     0.2829178E+03_r8,&  
         0.2824058E+03_r8,     0.2818940E+03_r8,     0.2813818E+03_r8,     0.2808704E+03_r8,     0.2803591E+03_r8,&  
         0.2798475E+03_r8,     0.2793358E+03_r8,     0.2788243E+03_r8,     0.2783129E+03_r8,     0.2778011E+03_r8,&  
         0.2772891E+03_r8,     0.2767775E+03_r8,     0.2762652E+03_r8,     0.2757530E+03_r8,     0.2752407E+03_r8,&  
         0.2747277E+03_r8,     0.2742143E+03_r8,     0.2737008E+03_r8,     0.2731870E+03_r8,     0.2726726E+03_r8,&  
         0.2721575E+03_r8,     0.2716425E+03_r8,     0.2711263E+03_r8,     0.2706094E+03_r8,     0.2700920E+03_r8,&  
         0.2695742E+03_r8,     0.2690556E+03_r8,     0.2685357E+03_r8,     0.2680150E+03_r8,     0.2674936E+03_r8,&  
         0.2669713E+03_r8,     0.2664474E+03_r8,     0.2659225E+03_r8,     0.2653969E+03_r8,     0.2648696E+03_r8,&  
         0.2643409E+03_r8,     0.2638110E+03_r8,     0.2632797E+03_r8,     0.2627468E+03_r8,     0.2622119E+03_r8,&  
         0.2616756E+03_r8,     0.2611375E+03_r8,     0.2605975E+03_r8,     0.2600553E+03_r8,     0.2595112E+03_r8,&  
         0.2589649E+03_r8,     0.2584164E+03_r8,     0.2578654E+03_r8,     0.2573119E+03_r8,     0.2567557E+03_r8,&  
         0.2561970E+03_r8,     0.2556352E+03_r8,     0.2550704E+03_r8,     0.2545025E+03_r8,     0.2539311E+03_r8,&  
         0.2533563E+03_r8,     0.2527778E+03_r8,     0.2521961E+03_r8,     0.2516101E+03_r8,     0.2510200E+03_r8,&  
         0.2504253E+03_r8,     0.2498265E+03_r8,     0.2492232E+03_r8,     0.2486143E+03_r8,     0.2480006E+03_r8,&  
         0.2473815E+03_r8,     0.2467567E+03_r8,     0.2461262E+03_r8,     0.2454895E+03_r8,     0.2448461E+03_r8/)  

    psaditmk(1:150,  150)= (/ &
         0.9601257E+03_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,     0.8999999E+10_r8,&  
         0.8999999E+10_r8,     0.8999999E+10_r8,     0.3217826E+03_r8,     0.3212120E+03_r8,     0.3206398E+03_r8,&  
         0.3200671E+03_r8,     0.3194944E+03_r8,     0.3189249E+03_r8,     0.3183556E+03_r8,     0.3177885E+03_r8,&  
         0.3172225E+03_r8,     0.3166573E+03_r8,     0.3160941E+03_r8,     0.3155325E+03_r8,     0.3149719E+03_r8,&  
         0.3144120E+03_r8,     0.3138543E+03_r8,     0.3132976E+03_r8,     0.3127422E+03_r8,     0.3121874E+03_r8,&  
         0.3116346E+03_r8,     0.3110828E+03_r8,     0.3105323E+03_r8,     0.3099831E+03_r8,     0.3094348E+03_r8,&  
         0.3088878E+03_r8,     0.3083420E+03_r8,     0.3077973E+03_r8,     0.3072536E+03_r8,     0.3067110E+03_r8,&  
         0.3061696E+03_r8,     0.3056293E+03_r8,     0.3050899E+03_r8,     0.3045518E+03_r8,     0.3040143E+03_r8,&  
         0.3034784E+03_r8,     0.3029431E+03_r8,     0.3024090E+03_r8,     0.3018759E+03_r8,     0.3013437E+03_r8,&  
         0.3008123E+03_r8,     0.3002818E+03_r8,     0.2997523E+03_r8,     0.2992235E+03_r8,     0.2986957E+03_r8,&  
         0.2981691E+03_r8,     0.2976432E+03_r8,     0.2971177E+03_r8,     0.2965934E+03_r8,     0.2960695E+03_r8,&  
         0.2955470E+03_r8,     0.2950248E+03_r8,     0.2945035E+03_r8,     0.2939828E+03_r8,     0.2934625E+03_r8,&  
         0.2929433E+03_r8,     0.2924248E+03_r8,     0.2919068E+03_r8,     0.2913894E+03_r8,     0.2908727E+03_r8,&  
         0.2903564E+03_r8,     0.2898404E+03_r8,     0.2893255E+03_r8,     0.2888107E+03_r8,     0.2882969E+03_r8,&  
         0.2877831E+03_r8,     0.2872701E+03_r8,     0.2867572E+03_r8,     0.2862443E+03_r8,     0.2857328E+03_r8,&  
         0.2852211E+03_r8,     0.2847098E+03_r8,     0.2841987E+03_r8,     0.2836882E+03_r8,     0.2831778E+03_r8,&  
         0.2826672E+03_r8,     0.2821569E+03_r8,     0.2816476E+03_r8,     0.2811378E+03_r8,     0.2806280E+03_r8,&  
         0.2801187E+03_r8,     0.2796093E+03_r8,     0.2790996E+03_r8,     0.2785904E+03_r8,     0.2780812E+03_r8,&  
         0.2775717E+03_r8,     0.2770617E+03_r8,     0.2765526E+03_r8,     0.2760425E+03_r8,     0.2755322E+03_r8,&  
         0.2750219E+03_r8,     0.2745116E+03_r8,     0.2740007E+03_r8,     0.2734893E+03_r8,     0.2729778E+03_r8,&  
         0.2724658E+03_r8,     0.2719531E+03_r8,     0.2714400E+03_r8,     0.2709267E+03_r8,     0.2704124E+03_r8,&  
         0.2698973E+03_r8,     0.2693819E+03_r8,     0.2688656E+03_r8,     0.2683485E+03_r8,     0.2678302E+03_r8,&  
         0.2673114E+03_r8,     0.2667914E+03_r8,     0.2662707E+03_r8,     0.2657482E+03_r8,     0.2652249E+03_r8,&  
         0.2647006E+03_r8,     0.2641750E+03_r8,     0.2636478E+03_r8,     0.2631192E+03_r8,     0.2625889E+03_r8,&  
         0.2620574E+03_r8,     0.2615239E+03_r8,     0.2609890E+03_r8,     0.2604521E+03_r8,     0.2599131E+03_r8,&  
         0.2593724E+03_r8,     0.2588292E+03_r8,     0.2582837E+03_r8,     0.2577361E+03_r8,     0.2571862E+03_r8,&  
         0.2566334E+03_r8,     0.2560779E+03_r8,     0.2555198E+03_r8,     0.2549586E+03_r8,     0.2543944E+03_r8,&  
         0.2538272E+03_r8,     0.2532566E+03_r8,     0.2526821E+03_r8,     0.2521042E+03_r8,     0.2515225E+03_r8,&  
         0.2509366E+03_r8,     0.2503468E+03_r8,     0.2497522E+03_r8,     0.2491533E+03_r8,     0.2485495E+03_r8,&  
         0.2479407E+03_r8,     0.2473265E+03_r8,     0.2467070E+03_r8,     0.2460818E+03_r8,     0.2454503E+03_r8/)  

  END SUBROUTINE RD
!
!  END MODULE module_calc_cape
!
!
!  MODULE module_index_severe_storm
!
!-------------------------------------------------------------------------------
  FUNCTION K_index(T0,nCols,kMax)
    IMPLICIT NONE
    ! k-index   
    !(http://www.cimms.ou.edu/~schultz/papers/doswell_schultz_indices.pdf george (1960))
    INTEGER      , INTENT(in) :: nCols
    INTEGER      , INTENT(in) :: kMax
    REAL(KIND=r8), INTENT(in) :: T0(nCols,kMax)
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: K_index(SIZE(T0,dim=1))
    REAL(KIND=r8)             :: pvs700(SIZE(T0,dim=1))
    REAL(KIND=r8)             :: pvs850(SIZE(T0,dim=1))
    REAL(KIND=r8)             :: DewPoint850(SIZE(T0,dim=1))
    REAL(KIND=r8)             :: DewPoint700(SIZE(T0,dim=1))
    n=SIZE(T0)
    K_index=0.0_r8
    DO i=1, n
      pvs850(i)=fpvsq(T0(i,k850))            !saturation vapor pressure in Pascals
      pvs700(i)=fpvsq(T0(i,k700))            !saturation vapor pressure in Pascals
    ENDDO
    DewPoint850=DewPoint(pvs850,nCols)
    DewPoint700=DewPoint(pvs700,nCols)
    DO i=1, n
      K_index(i) = (T0(i,k850) -  T0(i,k500)) + DewPoint850(i) -  (T0(i,k700) -  DewPoint700(i))
    ENDDO
  END FUNCTION K_index
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
  FUNCTION SWEAT_index(T0,U0,V0,nCols,kMax)
    IMPLICIT NONE
    ! k-index   
    !(http://www.cimms.ou.edu/~schultz/papers/doswell_schultz_indices.pdf george (1960))
    INTEGER      , INTENT(in) :: nCols
    INTEGER      , INTENT(in) :: kMax
    REAL(KIND=r8), INTENT(in) :: T0(nCols,kMax)
    REAL(KIND=r8), INTENT(in) :: U0(nCols,kMax)
    REAL(KIND=r8), INTENT(in) :: V0(nCols,kMax)
    REAL(KIND=r8) :: Vknots(nCols,kMax)!1 kt = 0.51444 ms-1   -  
                                       !  y   = X 
                                       !  Y=  X * (1 kt)/ (0.51444 ms-1) 
    REAL(KIND=r8)             :: DirV(nCols,kMax)
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: SWEAT_index(nCols)
    REAL(KIND=r8)             :: pvs850(nCols)
    REAL(KIND=r8)             :: DewPoint850(nCols)
    REAL(KIND=r8)             :: TT(nCols)
    INTEGER :: k

    DO k=1,kMax
      DO i=1,nCols
         Vknots(i,k)  = sqrt( (U0(i,k)*U0(i,k))  +  (V0(i,k)*V0(i,k)))/0.51444_r8
         IF(V0(i,k) == 0.0_r8)THEN
            DirV(i,k)=atan2(U0(i,k),1.0e-12_r8)
         ELSE
            DirV(i,k)=atan2(U0(i,k),V0(i,k))
         END IF
         IF(DirV(i,k) < 0.0_r8)DirV(i,k)=DirV(i,k)+(2.0_r8*3.14_r8)
      END DO
    END DO
    n=nCols
    SWEAT_index=0.0_r8
    DO i=1, n
      pvs850(i)=fpvsq(T0(i,k850))            !saturation vapor pressure in Pascals
    ENDDO
    DewPoint850=DewPoint(pvs850,nCols)
    TT=TT_index(T0,nCols,kMax)
    DO i=1, n
       SWEAT_index(i) = 20.0_r8*MAX(((TT(i)) - 49.0_r8),0.0_r8) + 12.0_r8*(DewPoint850(i)-273.16_r8) + &
                         2.0_r8* Vknots(i,k850) + Vknots(i,k500) + 125.0_r8*(sin(DirV(i,k500) - DirV(i,k850) )+ 0.2_r8)
       SWEAT_index(i) = MAX(SWEAT_index(i),0.0_r8)  
    ENDDO
  END FUNCTION SWEAT_index
!-------------------------------------------------------------------------------

  FUNCTION CT_index(T0,nCols,kMax)
    IMPLICIT NONE
    ! Cross totals  
    !(http://www.cimms.ou.edu/~schultz/papers/doswell_schultz_indices.pdf Miller (1972))
    INTEGER      , INTENT(in) :: nCols
    INTEGER      , INTENT(in) :: kMax
    REAL(KIND=r8), INTENT(in) :: T0(nCols,kMax)
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: CT_index(nCols)
    REAL(KIND=r8)             :: pvs850(nCols)
    REAL(KIND=r8)             :: DewPoint850(nCols)
    n=nCols
    CT_index=0.0_r8
    DO i=1, n
      pvs850(i)=fpvsq(T0(i,k850))            !saturation vapor pressure in Pascals
    ENDDO
    DewPoint850=DewPoint(pvs850,nCols)
    DO i=1, n
      CT_index(i) = (DewPoint850(i) -  T0(i,k500))
    ENDDO
  END FUNCTION CT_index

!-------------------------------------------------------------------------------
  FUNCTION TT_index(T0,nCols,kMax)
    IMPLICIT NONE
    ! Total Totals  
    !(http://www.cimms.ou.edu/~schultz/papers/doswell_schultz_indices.pdf miller (1972))
    INTEGER      , INTENT(in) :: nCols
    INTEGER      , INTENT(in) :: kMax
    REAL(KIND=r8), INTENT(in) :: T0(nCols,kMax)
    REAL(KIND=r8)             :: TT_index(nCols)
    REAL(KIND=r8)             :: VT(nCols)
    REAL(KIND=r8)             :: CT(nCols)
    VT=VT_index(T0,nCols,kMax)
    CT=CT_index(T0,nCols,kMax)
    TT_index=0.0_r8
    TT_index = (VT + CT)
  END FUNCTION TT_index

!-------------------------------------------------------------------------------
  FUNCTION VT_index(T0,nCols,kMax)
    IMPLICIT NONE
    ! Vertical Totals   
    !(http://www.cimms.ou.edu/~schultz/papers/doswell_schultz_indices.pdf Miller (1972))
    INTEGER      , INTENT(in) :: nCols
    INTEGER      , INTENT(in) :: kMax
    REAL(KIND=r8), INTENT(in) :: T0(nCols,kMax)
    INTEGER                   :: i
    REAL(KIND=r8)             :: VT_index(nCols)
    VT_index=0.0_r8
    DO i=1, nCols
      VT_index(i) = (T0(i,k850) -  T0(i,k500))
    ENDDO
  END FUNCTION VT_index

  FUNCTION DewPoint(es,nCols)
    IMPLICIT NONE
    !  4 DEWPOINT AND FROSTPOINT FORMULAS
    !  Equations 2 and 3 are easily solved for vapor pressures at any given temperature, namely the
    !  dewpoint and frostpoint temperatures.  However, if vapor pressure is known with temperature
    !  as the unknown desired quantity, the solution immediately becomes complicated and must be
    !  solved by iteration.  For ease of computation, inverse equations have been developed to yield
    !  temperature at a given vapor pressure.
    !  4.1 Dewpoint Formula
    !  Equation 2 with ITS-90 coefficients was used to create a table of 201 data points from �100 to
    !  100�C, at 1 degree intervals.  The data was equally weighted and fit to equation 4.  Agreement
    !  between this dewpoint formula and equation 2 with ITS-90 coefficients is better than 0.3 mK
    !  over the range of -�100   to 100�C.
    !  over the range of 173.15 to 373.15�K.
    !  4.1 Dewpoint Formula
    ! (http://www.thunderscientific.com/tech_info/reflibrary/its90formulas.pdf)
    ! es is the saturation vapor pressure in Pa
    INTEGER      , INTENT(in) :: nCols
    REAL(KIND=r8), INTENT(in) :: es(nCols)!saturation vapor pressure in Pa
    INTEGER                   :: i
    REAL(KIND=r8)             :: DewPoint(nCols)! Td is dewpoint temperature in Kelvin
    REAL(KIND=r8), PARAMETER  :: c0=  2.0798233e2_r8
    REAL(KIND=r8), PARAMETER  :: c1= -2.0156028e1_r8
    REAL(KIND=r8), PARAMETER  :: c2=  4.6778925e-1_r8
    REAL(KIND=r8), PARAMETER  :: c3= -9.2288067e-6_r8
    REAL(KIND=r8), PARAMETER  :: d0=  1.0_r8
    REAL(KIND=r8), PARAMETER  :: d1= -1.3319669e-1_r8
    REAL(KIND=r8), PARAMETER  :: d2=  5.6577518e-3_r8
    REAL(KIND=r8), PARAMETER  :: d3= -7.5172865e-5_r8
    DewPoint=0.0_r8
    DO i=1, nCols
      DewPoint(i)= ((c0*(log(es(i)))**0) + (c1*(log(es(i)))**1) + (c2*(log(es(i)))**2) + (c3*(log(es(i)))**3))/ &
                   ((d0*(log(es(i)))**0) + (d1*(log(es(i)))**1) + (d2*(log(es(i)))**2) + (d3*(log(es(i)))**3))
    ENDDO
  END FUNCTION DewPoint

!
!  END MODULE module_index_severe_storm
!

END MODULE PhysicalFunctions



!
!  $Author: pkubota $
!  $Date: 2011/07/14 15:19:14 $
!  $Revision: 1.25 $
!
!
!MCGACPTEC : MEDIATION_LAYER:PHYSICS
!

MODULE MONIN_OBUKHOV
  IMPLICIT NONE
  ! Selecting Kinds
  INTEGER, PARAMETER :: r4 = SELECTED_REAL_KIND(6)  ! Kind for 32-bits Real Numbers
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)   ! Kind for 32-bits Integer Numbers
  INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)  ! Kind for 64-bits Integer Numbers
  INTEGER, PARAMETER :: r16 = SELECTED_REAL_KIND(15)! Kind for 128-bits Real Numbers

  REAL(KIND=r8), PARAMETER :: cts_vonKarman=0.4_r8
  REAL(KIND=r8), PARAMETER :: cts_gravity=9.81_r8
  REAL(KIND=r8), PARAMETER :: cts_gases=287.0_r8
  REAL(KIND=r8), PARAMETER :: cts_heatespecific=1004.6_r8
  REAL(KIND=r8), PARAMETER :: cts_heatevap    = 2.52e6_r8! heat of evaporation of water     (j/kg) 

  REAL(KIND=r8), PARAMETER :: alfa=0.74_r8 
  REAL(KIND=r8), PARAMETER :: beta=4.7_r8
  REAL(KIND=r8), PARAMETER :: gama1=15.0_r8
  REAL(KIND=r8), PARAMETER :: gama2=9.0_r8 
  REAL(KIND=r8), PARAMETER :: sens_limear=10.0_r8
  REAL(KIND=r8), PARAMETER :: ustar_limearUpper=0.8_r8
  REAL(KIND=r8), PARAMETER :: ustar_limearLower=0.09_r8
  REAL(KIND=r8), PARAMETER :: z0_limearUpper=10.0_r8
  REAL(KIND=r8), PARAMETER :: z0_limearLower=0.09_r8
CONTAINS
  FUNCTION ScaleLength2(P0,T0,H0,Ustar,nCols)
    IMPLICIT NONE
    ! length scale   
    !(equation 11.1 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(in) :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: ScaleLength2(SIZE(P0))
    REAL(KIND=r8)             :: rho

    n=SIZE(P0)
    ScaleLength2=0.0_r8
    DO i=1, n
       rho = P0(i)/(cts_gases*T0(i))
       IF(H0(i) == 0.0_r8)THEN
           ScaleLength2(i) = - (Ustar(i)**3)/((cts_vonKarman*(cts_gravity/T0(i)))*(H0(i)/(rho*cts_heatespecific)))
       ELSE
           ScaleLength2(i) = - (Ustar(i)**3)/((cts_vonKarman*(cts_gravity/T0(i)))*((H0(i)+0.001_r8)/(rho*cts_heatespecific)))
       END IF  
    ENDDO
  END FUNCTION ScaleLength2

  FUNCTION PSIh2(P0,T0,H0,Ustar,Zlev,nCols)
    IMPLICIT NONE
    ! Simularity function for  Temperature Profile   
    !(equation 11.7 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(in) :: nCols 
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: PSIh2(SIZE(P0))
    REAL(KIND=r8)             :: boyance_parameter(SIZE(P0))
    REAL(KIND=r8)             :: rho
    n=SIZE(P0)
    boyance_parameter = zLev/ScaleLength2(P0,T0,H0,Ustar,nCols)
    PSIh2=0.0_r8
    DO i=1, n
       rho = P0(i)/(cts_gases*T0(i))        
       IF(boyance_parameter(i)< 0.0_r8)THEN
          PSIh2(i) = alfa*((1.0_r8-gama2*boyance_parameter(i))**(-1.0_r8/2.0_r8)) ! unstable
       ELSE
          PSIh2(i) = alfa + beta*boyance_parameter(i)                             ! stable
       END IF
    ENDDO
  END FUNCTION PSIh2


  FUNCTION PSIm2(P0,T0,H0,Ustar,Zlev,nCols)
    IMPLICIT NONE
    ! Simularity function for wind and Temperature Profile   
    !(equation 11.6 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(in) :: nCols 
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: PSIm2(SIZE(P0))
    REAL(KIND=r8)             :: boyance_parameter(SIZE(P0))
    REAL(KIND=r8)             :: rho
    n=SIZE(P0)
    boyance_parameter = zLev/ScaleLength2(P0,T0,H0,Ustar,nCols)
    PSIm2=0.0_r8
    DO i=1, n
       rho = P0(i)/(cts_gases*T0(i))        
       IF(boyance_parameter(i)< 0.0_r8)THEN
          PSIm2(i) = alfa*((1.0_r8-gama1*boyance_parameter(i))**(-1.0_r8/4.0_r8)) ! unstable
       ELSE
          PSIm2(i) = alfa + beta*boyance_parameter(i)                             ! stable
       END IF
    ENDDO
  END FUNCTION PSIm2



  FUNCTION RichNumber2(P0,T0,H0,Ustar,Zlev,nCols)
    IMPLICIT NONE
    ! Richardson Number   
    !(equation 11.4 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(in) :: nCols 
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: RichNumber2(SIZE(P0))
    REAL(KIND=r8)             :: XX(SIZE(P0))
    REAL(KIND=r8)             :: YY(SIZE(P0))
    REAL(KIND=r8)             :: boyance_parameter(SIZE(P0))
    n=SIZE(P0)
    boyance_parameter = zLev/ScaleLength2(P0,T0,H0,Ustar,nCols)
    RichNumber2=0.0_r8
    XX = PSIh2(P0,T0,H0,Ustar,Zlev,nCols)
    YY = PSIm2(P0,T0,H0,Ustar,Zlev,nCols)
    DO i=1, n 
       RichNumber2(i) =  boyance_parameter(i) * (XX(i)/(YY(i)**2))

    END DO
  END FUNCTION RichNumber2





  FUNCTION QSIm2(P0,T0,H0,Ustar,Zlev,nCols)
    IMPLICIT NONE
    ! Simularity function for wind  Profile   
    !(equation 11.14 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(IN) :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: QSIm2(SIZE(P0))
    REAL(KIND=r8)             :: xx   (SIZE(P0))
    REAL(KIND=r8)             :: boyance_parameter(SIZE(P0))
    n=SIZE(P0)
    boyance_parameter = zLev/ScaleLength2(P0,T0,H0,Ustar,nCols)
    QSIm2=0.0_r8
    DO i=1, n
       IF(boyance_parameter(i)< 0.0_r8)THEN
          xx(i)=(1.0_r8 -15.0_r8*boyance_parameter(i))**(1.0_r8/4.0_r8)
          QSIm2(i) = LOG(((1.0_r8 + xx(i)**2)/2.0_r8) * (((1 + xx(i))/2.0_r8)**2) ) ! unstable
       ELSE
          QSIm2(i) = -5.0_r8*boyance_parameter(i)                             ! stable
       END IF
    ENDDO
  END FUNCTION QSIm2



  FUNCTION QSIh2(P0,T0,H0,Ustar,Zlev,nCols)
    IMPLICIT NONE
    ! Simularity function for Temperature  Profile   
    !(equation 11.14 Arya P.S. Introduction to Micrometeorology)
    INTEGER                   :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: QSIh2(SIZE(P0))
    REAL(KIND=r8)             :: xx   (SIZE(P0))
    REAL(KIND=r8)             :: boyance_parameter(SIZE(P0))
    n=SIZE(P0)
    boyance_parameter = zLev/ScaleLength2(P0,T0,H0,Ustar,nCols)
    QSIh2=0.0_r8
    DO i=1, n
       IF(boyance_parameter(i)< 0.0_r8)THEN
          xx(i)=(1.0_r8 -15.0_r8*boyance_parameter(i))**(1.0_r8/4.0_r8)
          QSIh2(i) = 2.0_r8*LOG((1.0_r8 + xx(i)**2)/2.0_r8 ) ! unstable
       ELSE
          QSIh2(i) = -5.0_r8*boyance_parameter(i)                             ! stable
       END IF
    ENDDO
  END FUNCTION QSIh2


  
  FUNCTION T2MT(P0,T0,H0,Ustar,Zlev,z0,nCols)
    IMPLICIT NONE
    ! Temperature Profile Method  
    !(equation 11.17 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(IN) :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(in) :: z0(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: T2MT(SIZE(P0))
    REAL(KIND=r8)             :: xx       (SIZE(P0))
    REAL(KIND=r8)             :: sensi    (SIZE(P0))
    REAL(KIND=r8)             :: ustar2   (SIZE(P0))
    REAL(KIND=r8)             :: zrough    (SIZE(P0))
    REAL(KIND=r8)             :: rho
    REAL(KIND=r8)             :: tstar

    n=SIZE(P0)
    DO i=1, n
       IF(ABS(H0(i)) > sens_limear) THEN
          sensi (i) = H0(i)
          zrough(i) = MAX(MIN(z0(i),z0_limearUpper),z0_limearLower)
          ustar2(i) = MAX(MIN(Ustar(i),ustar_limearUpper),ustar_limearLower)
       ELSE
          IF(H0(i) /= 0.0_r8)THEN
             ustar2(i)= ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = (H0(i)/H0(i))*sens_limear
          ELSE
             ustar2(i)= ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = sens_limear
          END IF     
       END IF 
    END DO 
    xx                =  QSIh2(P0,T0,sensi,ustar2,Zlev,nCols)
    T2MT=0.0_r8
    DO i=1, n
       rho    = P0(i)/(cts_gases*T0(i))
       tstar  = -sensi(i)/(rho*cts_heatespecific*ustar2(i))

       T2MT(i) = (tstar/cts_vonKarman) * ( LOG(Zlev) -  xx(i)) &
            + (T0(i) -  (tstar/cts_vonKarman)*LOG(z0(i)))
       T2MT(i) = MIN(MAX(T2MT(i),T0(i)*0.97_r8),T0(i)*1.03_r8)

    ENDDO

  END FUNCTION T2MT


  FUNCTION Q2MT(P0,T0,Q0,H0,E0,Ustar,Zlev,z0,nCols)
    IMPLICIT NONE
    ! Humidity Profile Method  
    !(equation 11.17 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(IN) :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: Q0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: E0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(in) :: z0(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: Q2MT (SIZE(P0))
    REAL(KIND=r8)             :: xx   (SIZE(P0))
    REAL(KIND=r8)             :: sensi    (SIZE(P0))
    REAL(KIND=r8)             :: ustar2   (SIZE(P0))
    REAL(KIND=r8)             :: zrough    (SIZE(P0))
    REAL(KIND=r8)             :: rho
    REAL(KIND=r8)             :: qstar
    n=SIZE(P0)
    DO i=1, n
       IF(ABS(H0(i)) > sens_limear) THEN
          sensi (i) = H0(i)
          zrough(i) = MAX(MIN(z0(i),z0_limearUpper),z0_limearLower)
          ustar2(i) = MAX(MIN(Ustar(i),ustar_limearUpper),ustar_limearLower)
       ELSE
          IF(H0(i) /= 0.0_r8)THEN
             ustar2(i)= ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = (H0(i)/H0(i))*sens_limear
          ELSE
             ustar2(i)= ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = sens_limear
          END IF     
       END IF 
    END DO 
    xx =  QSIh2(P0,T0,sensi,ustar2,Zlev,nCols)
    Q2MT=0.0_r8
    DO i=1, n
       rho    = P0(i)/(cts_gases*T0(i))
       qstar  = -E0(i)/(rho*cts_heatevap*ustar2(i))

       Q2MT(i) = (qstar/cts_vonKarman) * ( LOG(Zlev) -  xx(i)) &
            + (Q0(i) -  (qstar/cts_vonKarman)*LOG(zrough(i)))

       Q2MT(i) = MIN(MAX(Q2MT(i),Q0(i)*0.5_r8),Q0(i)*1.5_r8)
    ENDDO
  END FUNCTION Q2MT



  FUNCTION W2MT(P0,T0,Speedm,H0,Ustar,Zlev,z0,nCols,MaskLand)
    IMPLICIT NONE
    ! Wind Profile Method  
    !(equation 11.17 Arya P.S. Introduction to Micrometeorology)
    INTEGER      , INTENT(IN) :: nCols
    REAL(KIND=r8), INTENT(in) :: P0(nCols)
    REAL(KIND=r8), INTENT(in) :: T0(nCols)
    REAL(KIND=r8), INTENT(in) :: H0(nCols)
    REAL(KIND=r8), INTENT(in) :: Ustar(nCols)
    REAL(KIND=r8), INTENT(in) :: z0(nCols)
    REAL(KIND=r8), INTENT(in) :: Speedm(nCols)
    REAL(KIND=r8), INTENT(IN) :: ZLev
    INTEGER                   :: n, i
    REAL(KIND=r8)             :: W2MT (SIZE(P0))
    REAL(KIND=r8)             :: xx   (SIZE(P0))
    REAL(KIND=r8)             :: sensi    (SIZE(P0))
    REAL(KIND=r8)             :: ustar2   (SIZE(P0))
    REAL(KIND=r8)             :: zrough    (SIZE(P0))
    INTEGER                   :: MaskLand
    n=SIZE(P0)
    DO i=1, n
       IF(ABS(H0(i)) > sens_limear) THEN
          sensi (i) = H0(i)
          zrough(i) = MAX(MIN(z0(i),z0_limearUpper),z0_limearLower)
          ustar2(i) = MAX(MIN(Ustar(i),ustar_limearUpper),ustar_limearLower)
       ELSE
          IF(H0(i) /= 0.0_r8)THEN
             ustar2(i)= ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = (H0(i)/H0(i))*sens_limear
          ELSE
             ustar2(i) = ustar_limearLower
             zrough(i) = z0_limearLower
             sensi (i) = sens_limear
          END IF     
       END IF 
    END DO 
    xx =  QSIm2(P0,T0,sensi,ustar2,Zlev,nCols)
    W2MT=0.0_r8
    DO i=1, n
       IF(MaskLand == 1)THEN
          W2MT(i) = (Ustar(i)/cts_vonKarman) * ( LOG(Zlev) -  xx(i)) &
            +(sqrt(Speedm(i)) - ( (Ustar(i)/cts_vonKarman)*LOG(z0(i))))
       ELSE
          W2MT(i) = (Ustar(i)/cts_vonKarman) * ( LOG(Zlev) -  xx(i)) &
            +(sqrt(Speedm(i)) - ( (Ustar(i)/cts_vonKarman)*LOG(z0(i))))
       END IF
       W2MT(i) = MIN(MAX(W2MT(i),Speedm(i)*0.1_r8),Speedm(i)*1.5_r8)
    ENDDO
  END FUNCTION W2MT


END MODULE MONIN_OBUKHOV

MODULE wv_saturation
 IMPLICIT NONE
    INTEGER, PARAMETER :: r8 = SELECTED_REAL_KIND(15) ! Kind for 64-bits Real Numbers
    REAL(KIND=r8),PARAMETER :: SHR_CONST_MWWV   = 18.016_r8       ! molecular weight water vapor
    REAL(KIND=r8),PARAMETER :: SHR_CONST_MWDAIR = 28.966_r8       ! molecular weight dry air ~ kg/kmole
    REAL(KIND=r8),PARAMETER :: epsilo = shr_const_mwwv/shr_const_mwdair ! ratio of h2o to dry air molecular weights 
    REAL(KIND=r8),PARAMETER :: SHR_CONST_LATVAP = 2.501e6_r8      ! latent heat of evaporation ~ J/kg
    REAL(KIND=r8),PARAMETER :: SHR_CONST_TKFRZ  = 273.16_r8       ! freezing T of fresh water ~ K (intentionally made == to TKTRIP)
    REAL(KIND=r8),PARAMETER :: SHR_CONST_AVOGAD = 6.02214e26_r8   ! Avogadro's number ~ molecules/kmole
    REAL(KIND=r8),PARAMETER :: SHR_CONST_BOLTZ  = 1.38065e-23_r8  ! Boltzmann's constant ~ J/K/molecule
    REAL(KIND=r8),PARAMETER :: SHR_CONST_LATICE = 3.337e5_r8      ! latent heat of fusion ~ J/kg
    REAL(KIND=r8), PARAMETER :: latice = shr_const_latice ! Latent heat of fusion

    REAL(KIND=r8),PARAMETER :: trice  =  20.00_r8         ! Trans range from es over h2o to es over ice
    REAL(KIND=r8),PARAMETER :: ttrice=trice
    REAL(KIND=r8),PARAMETER :: SHR_CONST_CPDAIR = 1.00464e3_r8    ! specific heat of dry air ~ J/kg/K
    REAL(KIND=r8),PARAMETER :: cpair = shr_const_cpdair  ! specific heat of dry air (J/K/kg)
    REAL(KIND=r8),PARAMETER :: cp    =cpair
    REAL(KIND=r8),PARAMETER :: latvap = shr_const_latvap ! Latent heat of vaporization
    REAL(KIND=r8),PARAMETER :: SHR_CONST_RGAS   = SHR_CONST_AVOGAD*SHR_CONST_BOLTZ ! Universal gas constant ~ J/K/kmole
    REAL(KIND=r8),PARAMETER :: SHR_CONST_RWV    = SHR_CONST_RGAS/SHR_CONST_MWWV    ! Water vapor gas constant ~ J/K/kg
    INTEGER, PARAMETER :: plenest = 250! length of saturation vapor pressure table
    REAL(KIND=r8), PUBLIC, PARAMETER :: tmelt = shr_const_tkfrz   ! Freezing point of water

    REAL(KIND=r8)           :: estbl(plenest)      ! table values of saturation vapor pressure
    REAL(KIND=r8),PARAMETER :: tmn  = 173.16_r8          ! Minimum temperature entry in table
    REAL(KIND=r8),PARAMETER :: tmx  = 375.16_r8          ! Maximum temperature entry in table
    REAL(KIND=r8),PARAMETER :: tmin=tmn       ! min temperature (K) for table
    REAL(KIND=r8),PARAMETER :: tmax= tmx      ! max temperature (K) for table
    LOGICAL ,PARAMETER :: icephs=.TRUE.  ! false => saturation vapor press over water only
    INTEGER ,PARAMETER ::  iterp =2             ! #iterations for precipitation calculation
    INTEGER  ::  k1mb  =1  ! index of the eta level near 1 mb
    REAL(KIND=r8)           :: pcf(6)     ! polynomial coeffs -> es transition water to ice
    REAL(KIND=r8), PRIVATE :: hlatf  = latice
    REAL(KIND=r8), PRIVATE :: hlatv  = latvap
    REAL(KIND=r8), PRIVATE :: rgasv  = SHR_CONST_RWV    ! Gas constant for water vapor
    REAL(KIND=r8), PRIVATE :: t0 = tmelt                ! approximate freezing temp

CONTAINS

  SUBROUTINE findsp (nCols,pver, q, t, p, tsp, qsp)
    IMPLICIT NONE
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !     find the wet bulb temperature for a given t and q
    !     in a longitude height section
    !     wet bulb temp is the temperature and spec humidity that is 
    !     just saturated and has the same enthalpy
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
    !
    ! Method: 
    ! a Newton method is used
    ! first guess uses an algorithm provided by John Petch from the UKMO
    ! we exclude points where the physical situation is unrealistic
    ! e.g. where the temperature is outside the range of validity for the
    !      saturation vapor pressure, or where the water vapor pressure
    !      exceeds the ambient pressure, or the saturation specific humidity is 
    !      unrealistic
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    !
    !     input arguments
    !
    INTEGER, INTENT(in) :: nCols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels

    REAL(KIND=r8), INTENT(in) :: q(nCols,pver)        ! water vapor (kg/kg)
    REAL(KIND=r8), INTENT(in) :: t(nCols,pver)        ! temperature (K)
    REAL(KIND=r8), INTENT(in) :: p(nCols,pver)        ! pressure    (Pa)
    !
    ! output arguments
    !
    REAL(KIND=r8), INTENT(out) :: tsp(nCols,pver)      ! saturation temp (K)
    REAL(KIND=r8), INTENT(out) :: qsp(nCols,pver)      ! saturation mixing ratio (kg/kg)
    !
    ! local variables
    !
    INTEGER i                 ! work variable
    INTEGER k                 ! work variable
    LOGICAL lflg              ! work variable
    INTEGER iter              ! work variable
    INTEGER l                 ! work variable
    LOGICAL :: error_found

    REAL(KIND=r8) omeps                ! 1 minus epsilon
    REAL(KIND=r8) trinv                ! work variable
    REAL(KIND=r8) es                   ! sat. vapor pressure
    REAL(KIND=r8) desdt                ! change in sat vap pressure wrt temperature
    !     real(KIND=r8) desdp                ! change in sat vap pressure wrt pressure
    REAL(KIND=r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(KIND=r8) dgdt                 ! work variable
    REAL(KIND=r8) g                    ! work variable
    REAL(KIND=r8) weight(nCols)        ! work variable
    REAL(KIND=r8) hlatsb               ! (sublimation)
    REAL(KIND=r8) hlatvp               ! (vaporization)
    REAL(KIND=r8) hltalt(nCols,pver)   ! lat. heat. of vap.
    REAL(KIND=r8) tterm                ! work var.
    REAL(KIND=r8) qs                   ! spec. hum. of water vapor
    REAL(KIND=r8) tc                   ! crit temp of transition to ice

    ! work variables
    REAL(KIND=r8) t1, q1, dt, dq
    REAL(KIND=r8) dtm, dqm
    REAL(KIND=r8) qvd, a1, tmp
    REAL(KIND=r8) rair
    REAL(KIND=r8) r1b, c1, c2, c3
    REAL(KIND=r8) denom
    REAL(KIND=r8) dttol
    REAL(KIND=r8) dqtol
    INTEGER doit(nCols) 
    REAL(KIND=r8) enin(nCols), enout(nCols)
    REAL(KIND=r8) tlim(nCols)
    REAL(KIND=r8) epsqs
    epsqs = epsilo
    k1mb=1
    omeps = 1.0_r8 - epsqs
    trinv = 1.0_r8/ttrice
    a1 = 7.5_r8*LOG(10.0_r8)
    rair =  287.04_r8
    c3 = rair*a1/cp
    dtm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dqm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    !  tmin = 173.16_r8 ! the coldest temperature we can deal with
    !
    ! max number of times to iterate the calculation
    iter = 50
    !
    DO k = k1mb,pver

       !
       ! first guess on the wet bulb temperature
       !
       DO i = 1,ncols

          ! limit the temperature range to that relevant to the sat vap pres tables

          tlim(i) = MIN(MAX(t(i,k),173.0_r8),373.0_r8)

          es = estblf(tlim(i))

          denom = p(i,k) - omeps*es
          qs = epsqs*es/denom
          doit(i) = 0
          enout(i) = 1.0_r8
          ! make sure a meaningful calculation is possible
          IF (p(i,k) > 5.0_r8*es .AND. qs > 0.0_r8 .AND. qs < 0.5_r8) THEN
             !
             ! Saturation specific humidity
             !
             qs = MIN(epsqs*es/denom,1.0_r8)
             !
             ! "generalized" analytic expression for t derivative of es
             !  accurate to within 1 percent for 173.16 < t < 373.16
             !
             ! Weighting of hlat accounts for transition from water to ice
             ! polynomial expression approximates difference between es over
             ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
             ! -40): required for accurate estimate of es derivative in transition
             ! range from ice to water also accounting for change of hlatv with t
             ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
             !
             tc     = tlim(i) - t0
             lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
             weight(i) = MIN(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0_r8*tc
             IF (tlim(i) < t0) THEN
                hltalt(i,k) = hlatsb
             ELSE
                hltalt(i,k) = hlatvp
             END IF
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

             ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36.0_r8)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
             !#ifdef DEBUG
             !             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
             !                write (6,*) ' relative humidity ', q(i,k)/qs
             !                write (6,*) ' first guess ', tsp(i,k)
             !             endif
             !#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
          ELSE
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1.0_r8
          ENDIF
       END DO   ! end do i
       !
       ! now iterate on first guess
       !
       DO l = 1, iter
          dtm = 0.0_r8
          dqm = 0.0_r8
          DO i = 1,ncols
             IF (doit(i) == 0) THEN
                es = estblf(tsp(i,k))
                !
                ! Saturation specific humidity
                !
                qs = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                !
                ! "generalized" analytic expression for t derivative of es
                ! accurate to within 1 percent for 173.16 < t < 373.16
                !
                ! Weighting of hlat accounts for transition from water to ice
                ! polynomial expression approximates difference between es over
                ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                ! -40): required for accurate estimate of es derivative in transition
                ! range from ice to water also accounting for change of hlatv with t
                ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                !
                tc     = tsp(i,k) - t0
                lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
                weight(i) = MIN(-tc*trinv,1.0_r8)
                hlatsb = hlatv + weight(i)*hlatf
                hlatvp = hlatv - 2369.0_r8*tc
                IF (tsp(i,k) < t0) THEN
                   hltalt(i,k) = hlatsb
                ELSE
                   hltalt(i,k) = hlatvp
                END IF
                IF (lflg) THEN
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                ELSE
                   tterm = 0.0_r8
                END IF
                desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                !              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                dgdt = -(cp + hltalt(i,k)*dqsdt)
                t1 = tsp(i,k) - g/dgdt
                dt = ABS(t1 - tsp(i,k))/t1
                tsp(i,k) = MAX(t1,tmin)
                es = estblf(tsp(i,k))
                q1 = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                q1=MAX(q1,1.e-12_r8)
                dq = ABS(q1 - qsp(i,k))/MAX(q1,1.e-12_r8)
                qsp(i,k) = q1
                !#ifdef DEBUG
                !               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                !                  write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                !               endif
                !#endif
                dtm = MAX(dtm,dt)
                dqm = MAX(dqm,dq)
                ! if converged at this point, exclude it from more iterations
                IF (dt < dttol .AND. dq < dqtol) THEN
                   doit(i) = 2
                ENDIF
                enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                ! bail out if we are too near the end of temp range
                IF (tsp(i,k) < 174.16_r8) THEN
                   doit(i) = 4
                ENDIF
             ELSE
             ENDIF
          END DO              ! do i = 1,ncols

          IF (dtm < dttol .AND. dqm < dqtol) THEN
             go to 10
          ENDIF

       END DO                 ! do l = 1,iter
10     CONTINUE

       error_found = .FALSE.
       IF (dtm > dttol .OR. dqm > dqtol) THEN
          DO i = 1,ncols
             IF (doit(i) == 0) error_found = .TRUE.
          END DO
          IF (error_found) THEN
             DO i = 1,ncols
                IF (doit(i) == 0) THEN
                   WRITE (6,*) ' findsp not converging at point i, k ', i, k
                   WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                   WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                   STOP 'FINDSP'
                ENDIF
             END DO
          ENDIF
       ENDIF
       DO i = 1,ncols
          IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
             error_found = .TRUE.
          ENDIF
       END DO
       IF (error_found) THEN
          DO i = 1,ncols
             IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
                WRITE (6,*) ' the enthalpy is not conserved for point ', &
                     i, k, enin(i), enout(i)
                WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                STOP 'FINDSP'
             ENDIF
          END DO
       ENDIF

    END DO                    ! level loop (k=1,pver)

    RETURN
  END SUBROUTINE findsp

  SUBROUTINE findsp_mask (istart,iend,nCols,pver, q, t, p,ierr, tsp, qsp)
    IMPLICIT NONE
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    !     find the wet bulb temperature for a given t and q
    !     in a longitude height section
    !     wet bulb temp is the temperature and spec humidity that is 
    !     just saturated and has the same enthalpy
    !     if q > qs(t) then tsp > t and qsp = qs(tsp) < q
    !     if q < qs(t) then tsp < t and qsp = qs(tsp) > q
    !
    ! Method: 
    ! a Newton method is used
    ! first guess uses an algorithm provided by John Petch from the UKMO
    ! we exclude points where the physical situation is unrealistic
    ! e.g. where the temperature is outside the range of validity for the
    !      saturation vapor pressure, or where the water vapor pressure
    !      exceeds the ambient pressure, or the saturation specific humidity is 
    !      unrealistic
    ! 
    ! Author: P. Rasch
    ! 
    !-----------------------------------------------------------------------
    !
    !     input arguments
    !

    INTEGER, INTENT(in) :: nCols                 ! number of columns (max)
    INTEGER, INTENT(in) :: pver                  ! number of vertical levels
    INTEGER, INTENT(IN   )    :: istart
    INTEGER, INTENT(IN   )    :: iend
    REAL(KIND=r8), INTENT(in) :: q(nCols,pver)        ! water vapor (kg/kg)
    REAL(KIND=r8), INTENT(in) :: t(nCols,pver)        ! temperature (K)
    REAL(KIND=r8), INTENT(in) :: p(nCols,pver)        ! pressure    (Pa)
    INTEGER, INTENT(IN   )    :: ierr(nCols)
    !
    ! output arguments
    !
    REAL(KIND=r8), INTENT(out) :: tsp(nCols,pver)      ! saturation temp (K)
    REAL(KIND=r8), INTENT(out) :: qsp(nCols,pver)      ! saturation mixing ratio (kg/kg)
    !
    ! local variables
    !
    INTEGER i                 ! work variable
    INTEGER k                 ! work variable
    LOGICAL lflg              ! work variable
    INTEGER iter              ! work variable
    INTEGER l                 ! work variable
    LOGICAL :: error_found

    REAL(KIND=r8) omeps                ! 1 minus epsilon
    REAL(KIND=r8) trinv                ! work variable
    REAL(KIND=r8) es                   ! sat. vapor pressure
    REAL(KIND=r8) desdt                ! change in sat vap pressure wrt temperature
    !     real(KIND=r8) desdp                ! change in sat vap pressure wrt pressure
    REAL(KIND=r8) dqsdt                ! change in sat spec. hum. wrt temperature
    REAL(KIND=r8) dgdt                 ! work variable
    REAL(KIND=r8) g                    ! work variable
    REAL(KIND=r8) weight(nCols)        ! work variable
    REAL(KIND=r8) hlatsb               ! (sublimation)
    REAL(KIND=r8) hlatvp               ! (vaporization)
    REAL(KIND=r8) hltalt(nCols,pver)   ! lat. heat. of vap.
    REAL(KIND=r8) tterm                ! work var.
    REAL(KIND=r8) qs                   ! spec. hum. of water vapor
    REAL(KIND=r8) tc                   ! crit temp of transition to ice

    ! work variables
    REAL(KIND=r8) t1, q1, dt, dq
    REAL(KIND=r8) dtm, dqm
    REAL(KIND=r8) qvd, a1, tmp
    REAL(KIND=r8) rair
    REAL(KIND=r8) r1b, c1, c2, c3
    REAL(KIND=r8) denom
    REAL(KIND=r8) dttol
    REAL(KIND=r8) dqtol
    INTEGER doit(nCols) 
    REAL(KIND=r8) enin(nCols), enout(nCols)
    REAL(KIND=r8) tlim(nCols)
    REAL(KIND=r8) epsqs
    epsqs = epsilo
    k1mb=1
    omeps = 1.0_r8 - epsqs
    trinv = 1.0_r8/ttrice
    a1 = 7.5_r8*LOG(10.0_r8)
    rair =  287.04_r8
    c3 = rair*a1/cp
    dtm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dqm = 0.0_r8    ! needed for iter=0 blowup with f90 -ei
    dttol = 1.e-4_r8 ! the relative temp error tolerance required to quit the iteration
    dqtol = 1.e-4_r8 ! the relative moisture error tolerance required to quit the iteration
    !  tmin = 173.16_r8 ! the coldest temperature we can deal with
    !
    ! max number of times to iterate the calculation
    iter = 50
    !
    DO k = k1mb,pver

       !
       ! first guess on the wet bulb temperature
       !
       DO i =istart,iend 
       IF(ierr(i) == 0)THEN    
          ! limit the temperature range to that relevant to the sat vap pres tables

          tlim(i) = MIN(MAX(t(i,k),173.0_r8),373.0_r8)

          es = estblf(tlim(i))

          denom = p(i,k) - omeps*es
          qs = epsqs*es/denom
          doit(i) = 0
          enout(i) = 1.0_r8
          ! make sure a meaningful calculation is possible
          IF (p(i,k) > 5.0_r8*es .AND. qs > 0.0_r8 .AND. qs < 0.5_r8) THEN
             !
             ! Saturation specific humidity
             !
             qs = MIN(epsqs*es/denom,1.0_r8)
             !
             ! "generalized" analytic expression for t derivative of es
             !  accurate to within 1 percent for 173.16 < t < 373.16
             !
             ! Weighting of hlat accounts for transition from water to ice
             ! polynomial expression approximates difference between es over
             ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
             ! -40): required for accurate estimate of es derivative in transition
             ! range from ice to water also accounting for change of hlatv with t
             ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
             !
             tc     = tlim(i) - t0
             lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
             weight(i) = MIN(-tc*trinv,1.0_r8)
             hlatsb = hlatv + weight(i)*hlatf
             hlatvp = hlatv - 2369.0_r8*tc
             IF (tlim(i) < t0) THEN
                hltalt(i,k) = hlatsb
             ELSE
                hltalt(i,k) = hlatvp
             END IF
             enin(i) = cp*tlim(i) + hltalt(i,k)*q(i,k)

             ! make a guess at the wet bulb temp using a UKMO algorithm (from J. Petch)
             tmp =  q(i,k) - qs
             c1 = hltalt(i,k)*c3
             c2 = (tlim(i) + 36.0_r8)**2
             r1b    = c2/(c2 + c1*qs)
             qvd   = r1b*tmp
             tsp(i,k) = tlim(i) + ((hltalt(i,k)/cp)*qvd)
             !#ifdef DEBUG
             !             if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
             !                write (6,*) ' relative humidity ', q(i,k)/qs
             !                write (6,*) ' first guess ', tsp(i,k)
             !             endif
             !#endif
             es = estblf(tsp(i,k))
             qsp(i,k) = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
          ELSE
             doit(i) = 1
             tsp(i,k) = tlim(i)
             qsp(i,k) = q(i,k)
             enin(i) = 1.0_r8
          ENDIF
       END IF
       END DO   ! end do i
       !
       ! now iterate on first guess
       !
       DO l = 1, iter
          dtm = 0.0_r8
          dqm = 0.0_r8
          DO i =istart,iend
          IF(ierr(i) == 0)THEN 
             IF (doit(i) == 0) THEN
                es = estblf(tsp(i,k))
                !
                ! Saturation specific humidity
                !
                qs = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                !
                ! "generalized" analytic expression for t derivative of es
                ! accurate to within 1 percent for 173.16 < t < 373.16
                !
                ! Weighting of hlat accounts for transition from water to ice
                ! polynomial expression approximates difference between es over
                ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
                ! -40): required for accurate estimate of es derivative in transition
                ! range from ice to water also accounting for change of hlatv with t
                ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
                !
                tc     = tsp(i,k) - t0
                lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
                weight(i) = MIN(-tc*trinv,1.0_r8)
                hlatsb = hlatv + weight(i)*hlatf
                hlatvp = hlatv - 2369.0_r8*tc
                IF (tsp(i,k) < t0) THEN
                   hltalt(i,k) = hlatsb
                ELSE
                   hltalt(i,k) = hlatvp
                END IF
                IF (lflg) THEN
                   tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3)+tc*(pcf(4) + tc*pcf(5))))
                ELSE
                   tterm = 0.0_r8
                END IF
                desdt = hltalt(i,k)*es/(rgasv*tsp(i,k)*tsp(i,k)) + tterm*trinv
                dqsdt = (epsqs + omeps*qs)/(p(i,k) - omeps*es)*desdt
                !              g = cp*(tlim(i)-tsp(i,k)) + hltalt(i,k)*q(i,k)- hltalt(i,k)*qsp(i,k)
                g = enin(i) - (cp*tsp(i,k) + hltalt(i,k)*qsp(i,k))
                dgdt = -(cp + hltalt(i,k)*dqsdt)
                t1 = tsp(i,k) - g/dgdt
                dt = ABS(t1 - tsp(i,k))/t1
                tsp(i,k) = MAX(t1,tmin)
                es = estblf(tsp(i,k))
                q1 = MIN(epsqs*es/(p(i,k) - omeps*es),1.0_r8)
                q1=MAX(q1,1.e-12_r8)
                dq = ABS(q1 - qsp(i,k))/MAX(q1,1.e-12_r8)
                qsp(i,k) = q1
                !#ifdef DEBUG
                !               if ( (lchnk == lchnklook(nlook) ) .and. (i == icollook(nlook) ) ) then
                !                  write (6,*) ' rel chg lev, iter, t, q ', k, l, dt, dq, g
                !               endif
                !#endif
                dtm = MAX(dtm,dt)
                dqm = MAX(dqm,dq)
                ! if converged at this point, exclude it from more iterations
                IF (dt < dttol .AND. dq < dqtol) THEN
                   doit(i) = 2
                ENDIF
                enout(i) = cp*tsp(i,k) + hltalt(i,k)*qsp(i,k)
                ! bail out if we are too near the end of temp range
                IF (tsp(i,k) < 174.16_r8) THEN
                   doit(i) = 4
                ENDIF
             ELSE
             ENDIF
          END IF
          END DO              ! do i = 1,ncols

          IF (dtm < dttol .AND. dqm < dqtol) THEN
             go to 10
          ENDIF

       END DO                 ! do l = 1,iter
10     CONTINUE

       error_found = .FALSE.
       IF (dtm > dttol .OR. dqm > dqtol) THEN
          DO i =istart,iend
            IF( ierr(i) == 0)THEN
             IF (doit(i) == 0) error_found = .TRUE.
            END IF
          END DO
          IF (error_found) THEN
             DO i =istart,iend
               IF(ierr(i) == 0)THEN
                IF (doit(i) == 0) THEN
                   WRITE (6,*) ' findsp not converging at point i, k ', i, k
                   WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                   WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                   STOP 'FINDSP'
                ENDIF
               END IF
             END DO
          ENDIF
       ENDIF
       DO i =istart,iend
        IF(ierr(i) == 0)THEN
          IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
             error_found = .TRUE.
          ENDIF
        END IF
       END DO
       IF (error_found) THEN
          DO i =istart,iend
            IF(ierr(i) == 0)THEN  
             IF (doit(i) == 2 .AND. ABS((enin(i)-enout(i))/(enin(i)+enout(i))) > 1.e-4_r8) THEN
                WRITE (6,*) ' the enthalpy is not conserved for point ', &
                     i, k, enin(i), enout(i)
                WRITE (6,*) ' t, q, p, enin ', t(i,k), q(i,k), p(i,k), enin(i)
                WRITE (6,*) ' tsp, qsp, enout ', tsp(i,k), qsp(i,k), enout(i)
                STOP 'FINDSP'
             ENDIF
             END IF
          END DO
       ENDIF

    END DO                    ! level loop (k=1,pver)

    RETURN
  END SUBROUTINE findsp_mask


  SUBROUTINE gestbl(       )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Builds saturation vapor pressure table for later lookup procedure.
    ! 
    ! Method: 
    ! Uses Goff & Gratch (1946) relationships to generate the table
    ! according to a set of free parameters defined below.  Auxiliary
    ! routines are also included for making rapid estimates (well with 1%)
    ! of both es and d(es)/dt for the particular table configuration.
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use pmgrid, only: masterproc
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(KIND=r8)  epsqs
    REAL(KIND=r8) t             ! Temperature
    INTEGER n          ! Increment counter
    INTEGER lentbl     ! Calculated length of lookup table
    INTEGER itype      ! Ice phase: 0 -> no ice phase
    !            1 -> ice phase, no transitiong
    !           -x -> ice phase, x degree transition
    !
    !-----------------------------------------------------------------------
    !
    ! Set es table parameters
    !
    !   tmin   = tmn       ! Minimum temperature entry in table
    !   tmax   = tmx       ! Maximum temperature entry in table
    !   ttrice = trice     ! Trans. range from es over h2o to es over ice
    !   icephs = ip        ! Ice phase (true or false)
    !
    ! Set physical constants required for es calculation
    !
    epsqs  = epsilo
    !
    lentbl = INT(tmax-tmin+2.000001_r8)
    IF (lentbl .GT. plenest) THEN
       WRITE(6,9000) tmax, tmin, plenest
       STOP 'GESTBL'    ! Abnormal termination
    END IF
    !
    ! Begin building es table.
    ! Check whether ice phase requested.
    ! If so, set appropriate transition range for temperature
    !
    IF (icephs) THEN
       IF (ttrice /= 0.0_r8) THEN
          itype = -ttrice
       ELSE
          itype = 1
       END IF
    ELSE
       itype = 0
    END IF
    !
    !tmn  = 173.16_r8! Minimum temperature entry in table
    !tmx  = 375.16_r8! Maximum temperature entry in table

    t = tmin - 1.0_r8
    DO n=1,lentbl
       t = t + 1.0_r8
       CALL gffgch(t,estbl(n),itype)
    END DO
    !
    DO n=lentbl+1,plenest
       estbl(n) = -99999.0_r8
    END DO
    !
    ! Table complete -- Set coefficients for polynomial approximation of
    ! difference between saturation vapor press over water and saturation
    ! pressure over ice for -ttrice < t < 0 (degrees C). NOTE: polynomial
    ! is valid in the range -40 < t < 0 (degrees C).
    !
    !                  --- Degree 5 approximation ---
    !
    pcf(1) =  5.04469588506e-01_r8
    pcf(2) = -5.47288442819e+00_r8
    pcf(3) = -3.67471858735e-01_r8
    pcf(4) = -8.95963532403e-03_r8
    pcf(5) = -7.78053686625e-05_r8
    !
    !                  --- Degree 6 approximation ---
    !
    !-----pcf(1) =  7.63285250063e-02
    !-----pcf(2) = -5.86048427932e+00
    !-----pcf(3) = -4.38660831780e-01
    !-----pcf(4) = -1.37898276415e-02
    !-----pcf(5) = -2.14444472424e-04
    !-----pcf(6) = -1.36639103771e-06
    !
    !   if (masterproc) then
    !      write(6,*)' *** SATURATION VAPOR PRESSURE TABLE COMPLETED ***'
    !   end if

    RETURN
    !
9000 FORMAT('GESTBL: FATAL ERROR *********************************',/, &
         ' TMAX AND TMIN REQUIRE A LARGER DIMENSION ON THE LENGTH', &
         ' OF THE SATURATION VAPOR PRESSURE TABLE ESTBL(PLENEST)',/, &
         ' TMAX, TMIN, AND PLENEST => ', 2f7.2, i3)
    !
  END SUBROUTINE gestbl
  
  
  REAL(KIND=r8) FUNCTION estblf( td )
    !
    ! Saturation vapor pressure table lookup
    !
    REAL(KIND=r8), INTENT(in) :: td         ! Temperature for saturation lookup
    !
    REAL(KIND=r8) :: e          ! intermediate variable for es look-up
    REAL(KIND=r8) :: tmin       ! min temperature (K) for table
    REAL(KIND=r8) :: tmax       ! max temperature (K) for table

    REAL(KIND=r8) :: ai
    INTEGER  :: i

    tmin=tmn  
    tmax=tmx  

    !
    e = MAX(MIN(td,tmax),tmin)   ! partial pressure
    i = INT(e-tmin)+1
    ai = AINT(e-tmin)
    estblf = (tmin+ai-e+1.0_r8)* &
         estbl(i)-(tmin+ai-e)* &
         estbl(i+1)
  END FUNCTION estblf
  
  
  SUBROUTINE aqsat(t       ,p       ,es      ,qs        ,ii      , &
       ILEN    ,kk      ,kstart  ,kend      )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g),for input arrays of temperature and pressure (dimensioned ii,kk)
    ! This routine is useful for evaluating only a selected region in the
    ! vertical.
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER  , INTENT(in) :: ii             ! I dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: kk             ! K dimension of arrays t, p, es, qs
    INTEGER  , INTENT(in) :: ILEN           ! Length of vectors in I direction which
    INTEGER  , INTENT(in) :: kstart         ! Starting location in K direction
    INTEGER  , INTENT(in) :: kend           ! Ending location in K direction
    REAL(KIND=r8), INTENT(in) :: t(ii,kk)          ! Temperature [K]
    REAL(KIND=r8), INTENT(in) :: p(ii,kk)          ! Pressure [Pa]
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(out) :: es(ii,kk)         ! Saturation vapor pressure [Pa]
    REAL(KIND=r8), INTENT(out) :: qs(ii,kk)         ! Saturation specific humidity [kg/kg]
    !
    !---------------------------Local workspace-----------------------------
    !
    REAL(KIND=r8) epsqs      ! Ratio of h2o to dry air molecular weights 
    REAL(KIND=r8) omeps             ! 1 - 0.622
    INTEGER i, k           ! Indices
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo
    omeps = 1.0_r8 - epsqs
    DO k=kstart,kend
       DO i=1,ILEN
          es(i,k) = estblf(t(i,k))
          !
          ! Saturation specific humidity
          !
          qs(i,k) = epsqs*es(i,k)/(p(i,k) - omeps*es(i,k))
          !
          ! The following check is to avoid the generation of negative values
          ! that can occur in the upper stratosphere and mesosphere
          !
          qs(i,k) = MIN(1.0_r8,qs(i,k))
          !
          IF (qs(i,k) < 0.0_r8) THEN
             qs(i,k) = 1.0_r8
             es(i,k) = p(i,k)
          END IF
       END DO
    END DO
    !
    RETURN
  END SUBROUTINE aqsat
  SUBROUTINE vqsatd(t       ,p       ,es      ,qs      ,gam      , &
       len     )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Utility procedure to look up and return saturation vapor pressure from
    ! precomputed table, calculate and return saturation specific humidity
    ! (g/g), and calculate and return gamma (l/cp)*(d(qsat)/dT).  The same
    ! function as qsatd, but operates on vectors of temperature and pressure
    ! 
    ! Method: 
    ! 
    ! Author: J. Hack
    ! 
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    INTEGER, INTENT(in) :: len       ! vector length
    REAL(KIND=r8), INTENT(in) :: t(len)       ! temperature [K]
    REAL(KIND=r8), INTENT(in) :: p(len)       ! pressure [Pa]
    !
    ! Output arguments
    !
    REAL(KIND=r8), INTENT(out) :: es(len)   ! saturation vapor pressure
    REAL(KIND=r8), INTENT(out) :: qs(len)   ! saturation specific humidity
    REAL(KIND=r8), INTENT(out) :: gam(len)  ! (l/cp)*(d(qs)/dt)
    !
    !--------------------------Local Variables------------------------------
    !
    LOGICAL lflg        ! true if in temperature transition region
    !
    INTEGER i           ! index for vector calculations
    !
    REAL(KIND=r8) omeps     ! 1. - 0.622
    REAL(KIND=r8) trinv     ! reciprocal of ttrice (transition range)
    REAL(KIND=r8) tc        ! temperature (in degrees C)
    REAL(KIND=r8) weight    ! weight for es transition from water to ice
    REAL(KIND=r8) hltalt    ! appropriately modified hlat for T derivatives
    !
    REAL(KIND=r8) hlatsb    ! hlat weighted in transition region
    REAL(KIND=r8) hlatvp    ! hlat modified for t changes above freezing
    REAL(KIND=r8) tterm     ! account for d(es)/dT in transition region
    REAL(KIND=r8) desdt     ! d(es)/dT
    REAL(KIND=r8) epsqs
    !
    !-----------------------------------------------------------------------
    !
    epsqs = epsilo

    omeps = 1.0_r8 - epsqs
    DO i=1,len
       es(i) = estblf(t(i))
       !
       ! Saturation specific humidity
       !
       qs(i) = epsqs*es(i)/(p(i) - omeps*es(i))
       !
       ! The following check is to avoid the generation of negative
       ! values that can occur in the upper stratosphere and mesosphere
       !
       qs(i) = MIN(1.0_r8,qs(i))
       !
       IF (qs(i) < 0.0_r8) THEN
          qs(i) = 1.0_r8
          es(i) = p(i)
       END IF
    END DO
    !
    ! "generalized" analytic expression for t derivative of es
    ! accurate to within 1 percent for 173.16 < t < 373.16
    !
    trinv = 0.0_r8
    IF ((.NOT. icephs) .OR. (ttrice.EQ.0.0_r8)) go to 10
    trinv = 1.0_r8/ttrice
    DO i=1,len
       !
       ! Weighting of hlat accounts for transition from water to ice
       ! polynomial expression approximates difference between es over
       ! water and es over ice from 0 to -ttrice (C) (min of ttrice is
       ! -40): required for accurate estimate of es derivative in transition
       ! range from ice to water also accounting for change of hlatv with t
       ! above freezing where const slope is given by -2369 j/(kg c) = cpv - cw
       !
       tc     = t(i) - tmelt
       lflg   = (tc >= -ttrice .AND. tc < 0.0_r8)
       weight = MIN(-tc*trinv,1.0_r8)
       hlatsb = hlatv + weight*hlatf
       hlatvp = hlatv - 2369.0_r8*tc
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       IF (lflg) THEN
          tterm = pcf(1) + tc*(pcf(2) + tc*(pcf(3) + tc*(pcf(4) + tc*pcf(5))))
       ELSE
          tterm = 0.0_r8
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i)) + tterm*trinv
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    RETURN
    !
    ! No icephs or water to ice transition
    !
10  DO i=1,len
       !
       ! Account for change of hlatv with t above freezing where
       ! constant slope is given by -2369 j/(kg c) = cpv - cw
       !
       hlatvp = hlatv - 2369.0_r8*(t(i)-tmelt)
       IF (icephs) THEN
          hlatsb = hlatv + hlatf
       ELSE
          hlatsb = hlatv
       END IF
       IF (t(i) < tmelt) THEN
          hltalt = hlatsb
       ELSE
          hltalt = hlatvp
       END IF
       desdt  = hltalt*es(i)/(rgasv*t(i)*t(i))
       gam(i) = hltalt*qs(i)*p(i)*desdt/(cp*es(i)*(p(i) - omeps*es(i)))
       IF (qs(i) == 1.0_r8) gam(i) = 0.0_r8
    END DO
    !
    RETURN
    !
  END SUBROUTINE vqsatd

 SUBROUTINE gffgch(t       ,es      ,itype   )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Computes saturation vapor pressure over water and/or over ice using
    ! Goff & Gratch (1946) relationships. 
    ! <Say what the routine does> 
    ! 
    ! Method: 
    ! T (temperature), and itype are input parameters, while es (saturation
    ! vapor pressure) is an output parameter.  The input parameter itype
    ! serves two purposes: a value of zero indicates that saturation vapor
    ! pressures over water are to be returned (regardless of temperature),
    ! while a value of one indicates that saturation vapor pressures over
    ! ice should be returned when t is less than freezing degrees.  If itype
    ! is negative, its absolute value is interpreted to define a temperature
    ! transition region below freezing in which the returned
    ! saturation vapor pressure is a weighted average of the respective ice
    ! and water value.  That is, in the temperature range 0 => -itype
    ! degrees c, the saturation vapor pressures are assumed to be a weighted
    ! average of the vapor pressure over supercooled water and ice (all
    ! water at 0 c; all ice at -itype c).  Maximum transition range => 40 c
    ! 
    ! Author: J. Hack
    ! 
    !-----------------------------------------------------------------------
    !   use shr_kind_mod, only: r8 => shr_kind_r8
    !   use physconst, only: tmelt
    !   use abortutils, only: endrun

    IMPLICIT NONE
    !------------------------------Arguments--------------------------------
    !
    ! Input arguments
    !
    REAL(KIND=r8), INTENT(in) :: t          ! Temperature
    !
    ! Output arguments
    !
    INTEGER, INTENT(inout) :: itype   ! Flag for ice phase and associated transition

    REAL(KIND=r8), INTENT(out) :: es         ! Saturation vapor pressure
    !
    !---------------------------Local variables-----------------------------
    !
    REAL(KIND=r8) e1         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) e2         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) eswtr      ! Saturation vapor pressure over water
    REAL(KIND=r8) f          ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f1         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f2         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f3         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f4         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) f5         ! Intermediate scratch variable for es over water
    REAL(KIND=r8) ps         ! Reference pressure (mb)
    REAL(KIND=r8) t0         ! Reference temperature (freezing point of water)
    REAL(KIND=r8) term1      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) term2      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) term3      ! Intermediate scratch variable for es over ice
    REAL(KIND=r8) tr         ! Transition range for es over water to es over ice
    REAL(KIND=r8) ts         ! Reference temperature (boiling point of water)
    REAL(KIND=r8) weight     ! Intermediate scratch variable for es transition
    INTEGER itypo   ! Intermediate scratch variable for holding itype
    !
    !-----------------------------------------------------------------------
    !
    ! Check on whether there is to be a transition region for es
    !
    IF (itype < 0) THEN
       tr    = ABS(real(itype,kind=r8))
       itypo = itype
       itype = 1
    ELSE
       tr    = 0.0_r8
       itypo = itype
    END IF
    IF (tr > 40.0_r8) THEN
       WRITE(6,900) tr
       STOP 'GFFGCH'                ! Abnormal termination
    END IF
    !
    IF(t < (tmelt - tr) .AND. itype == 1) go to 10
    !
    ! Water
    !
    ps = 1013.246_r8
    ts = 373.16_r8
    e1 = 11.344_r8*(1.0_r8 - t/ts)
    e2 = -3.49149_r8*(ts/t - 1.0_r8)
    f1 = -7.90298_r8*(ts/t - 1.0_r8)
    f2 = 5.02808_r8*LOG10(ts/t)
    f3 = -1.3816_r8*(10.0_r8**e1 - 1.0_r8)/10000000.0_r8
    f4 = 8.1328_r8*(10.0_r8**e2 - 1.0_r8)/1000.0_r8
    f5 = LOG10(ps)
    f  = f1 + f2 + f3 + f4 + f5
    es = (10.0_r8**f)*100.0_r8
    eswtr = es
    !
    IF(t >= tmelt .OR. itype == 0) go to 20
    !
    ! Ice
    !
10  CONTINUE
    t0    = tmelt
    term1 = 2.01889049_r8/(t0/t)
    term2 = 3.56654_r8*LOG(t0/t)
    term3 = 20.947031_r8*(t0/t)
    es    = 575.185606e10_r8*EXP(-(term1 + term2 + term3))
    !
    IF (t < (tmelt - tr)) go to 20
    !
    ! Weighted transition between water and ice
    !
    weight = MIN((tmelt - t)/tr,1.0_r8)
    es = weight*es + (1.0_r8 - weight)*eswtr
    !
20  CONTINUE
    itype = itypo
    RETURN
    !
900 FORMAT('GFFGCH: FATAL ERROR ******************************',/, &
         'TRANSITION RANGE FOR WATER TO ICE SATURATION VAPOR', &
         ' PRESSURE, TR, EXCEEDS MAXIMUM ALLOWABLE VALUE OF', &
         ' 40.0 DEGREES C',/, ' TR = ',f7.2)
    !
  END SUBROUTINE gffgch
END MODULE wv_saturation
