!program test_rna
module test_rna
    use rna_class
    use rna_arch

    implicit none

contains
 
!subroutine pblh_rna(plond,    &
!                    plev,     &
!                    pcnst,    &
!                    th, q, z, u, v, t, pmid, kvf, shflx, taux, tauy, tsk, qsfc, psomc, pblh)
!
!    IMPLICIT NONE
!    INTEGER, INTENT(IN   ) :: plond      ! slt extended domain longitude
!    INTEGER, INTENT(IN   ) :: plev       ! number of vertical levels
!    INTEGER, INTENT(IN   ) :: pcnst      ! number of constituents (including water vapor)
!
!    REAL(KIND=r8), INTENT(IN   ) ::  th(plond,plev)          ! potential temperature [K]
!    REAL(KIND=r8), INTENT(IN   ) ::  q(plond,plev,pcnst)     ! specific humidity [kg/kg]
!    REAL(KIND=r8), INTENT(IN   ) ::  z(plond,plev)           ! height above surface [m]
!    REAL(KIND=r8), INTENT(IN   ) ::  u(plond,plev)           ! windspeed x-direction [m/s]
!    REAL(KIND=r8), INTENT(IN   ) ::  v(plond,plev)           ! windspeed y-direction [m/s]
!    REAL(KIND=r8), INTENT(IN   ) ::  t(plond,plev)           ! temperature (used for density)
!    REAL(KIND=r8), INTENT(IN   ) ::  pmid(plond,plev)        ! midpoint pressures
!    REAL(KIND=r8), INTENT(IN   ) ::  kvf(plond,plev + 1)        ! free atmospheric eddy diffsvty [m2/s]
!    REAL(KIND=r8), INTENT(IN   ) ::  shflx(plond)            ! surface heat flux (W/m2)
!    REAL(KIND=r8), INTENT(IN   ) ::  taux(plond)             ! surface u stress (N)
!    REAL(KIND=r8), INTENT(IN   ) ::  tauy(plond)             ! surface v stress (N)
!    REAL(KIND=r8), INTENT(IN   ) ::  tsk   (plond)
!    REAL(KIND=r8), INTENT(IN   ) ::  qsfc  (plond)
!    REAL(KIND=r8), INTENT(IN   ) ::  psomc(plond,plev)      ! (psm1/pmidm1)**cappa
!    !
!    ! Output arguments
!    !
!    REAL(KIND=r8), INTENT(INOUT ) ::  pblh   (plond)             ! boundary-layer height [m]
!
!
!    type(rna_AI) :: AI
!    integer      :: activationFunction
!    integer      :: nxl, nxc, nwh1l, nwh1c, nwsl, nwsc, nbh1l, nbh1c, nbsl, nbsc 
!
!    print*, 'inside pblh_rna' 
!
!    !read entries of the RNA
!    !AI%ch_id = 'pblh      '
!    AI%ch_id = '      '
!    call AI%get
!    allocate(valores_normalizados(max(1,AI%nlines),max(1,AI%ncolumns)))
!    allocate(valores_denormalizados(max(1,AI%nlines),max(1,AI%ncolumns)))
!    !save some data
!    nxl=AI%nlines
!    nxc=AI%ncolumns
!    !Onxl=AI%nlines
!    !Onxc=AI%ncolumns
!
!    call normaliza(AI%nlines, AI%ncolumns, AI%sp_entries, valores_normalizados)
!
!    !read the weights (wh1)
!    !AI%ch_id = 'pblh/wh1'
!    AI%ch_id = 'wh1'
!    call AI%weights
!    allocate(wh1(max(1,AI%nlines),max(1,AI%ncolumns)))
!    wh1=AI%sp_weights
!    !save some data
!    nwh1l=AI%nlines
!    nwh1c=AI%ncolumns
!
!    !read the weights (ws)
!    !AI%ch_id = 'pblh/ws'
!    AI%ch_id = 'ws'
!    call AI%weights
!    allocate(ws(max(1,AI%nlines),max(1,AI%ncolumns)))
!    ws=AI%sp_weights
!    !save some data
!    nwsl=AI%nlines
!    nwsc=AI%ncolumns
!
!    !read the weights (bh1)
!    !AI%ch_id = 'pblh/bh1'
!    AI%ch_id = 'bh1'
!    call AI%weights
!    allocate(bh1(max(1,AI%nlines),max(1,AI%ncolumns)))
!    bh1=AI%sp_weights
!    nbh1l=AI%nlines
!    nbh1c=AI%ncolumns
!
!    !read the weights (bs)
!    !AI%ch_id = 'pblh/bs'
!    AI%ch_id = 'bs'
!    call AI%weights
!    allocate(bs(max(1,AI%nlines),max(1,AI%ncolumns)))
!    bs=AI%sp_weights
!    nbsl=AI%nlines
!    nbsc=AI%ncolumns
!
!
!    !integer      :: OactivationFunction
!    !integer      :: Onxl, Onxc, Onwh1l, Onwh1c, Onwsl, Onwsc, Onbh1l, Onbh1c, Onbsl, Onbsc 
!
!   !with saved data allocate
!   allocate(vh1(nwh1l,nxc))
!   allocate(yh1(nwh1l,nxc))
!   allocate(vs(nxl, nwh1c))
!   allocate(ys(nxl, nwsc))
!
!   call feedfoward(activationFunction, NXL, NXC, NWH1L, NWH1C,&
!             NBH1L, NBH1C, NWSL, NWSC, NBSL,&
!             NBSC, VALORES_NORMALIZADOS, WH1, BH1, WS, BS, VH1, YH1, VS, YS, VALORES_DENORMALIZADOS)
!
!   
!   deallocate(vh1)
!   deallocate(yh1)
!   deallocate(vs)
!   deallocate(ys)
!end subroutine pblh_rna

subroutine pblh_rna_old
    type(rna_AI) :: AI
    integer      :: activationFunction
    integer      :: nxl, nxc, nwh1l, nwh1c, nwsl, nwsc, nbh1l, nbh1c, nbsl, nbsc 
    !integer      :: OactivationFunction
    !integer      :: Onxl, Onxc, Onwh1l, Onwh1c, Onwsl, Onwsc, Onbh1l, Onbh1c, Onbsl, Onbsc 

    print*, 'inside pblh_rna' 

    !read entries of the RNA
    !AI%ch_id = 'pblh      '
    AI%ch_id = '      '
    call AI%get
    allocate(valores_normalizados(max(1,AI%nlines),max(1,AI%ncolumns)))
    allocate(valores_denormalizados(max(1,AI%nlines),max(1,AI%ncolumns)))
    !save some data
    nxl=AI%nlines
    nxc=AI%ncolumns
    !Onxl=AI%nlines
    !Onxc=AI%ncolumns

    call normaliza(AI%nlines, AI%ncolumns, AI%sp_entries, valores_normalizados)

    !read the weights (wh1)
    !AI%ch_id = 'pblh/wh1'
    AI%ch_id = 'wh1'
    call AI%weights
    allocate(wh1(max(1,AI%nlines),max(1,AI%ncolumns)))
    wh1=AI%sp_weights
    !save some data
    nwh1l=AI%nlines
    nwh1c=AI%ncolumns

    !read the weights (ws)
    !AI%ch_id = 'pblh/ws'
    AI%ch_id = 'ws'
    call AI%weights
    allocate(ws(max(1,AI%nlines),max(1,AI%ncolumns)))
    ws=AI%sp_weights
    !save some data
    nwsl=AI%nlines
    nwsc=AI%ncolumns

    !read the weights (bh1)
    !AI%ch_id = 'pblh/bh1'
    AI%ch_id = 'bh1'
    call AI%weights
    allocate(bh1(max(1,AI%nlines),max(1,AI%ncolumns)))
    bh1=AI%sp_weights
    nbh1l=AI%nlines
    nbh1c=AI%ncolumns

    !read the weights (bs)
    !AI%ch_id = 'pblh/bs'
    AI%ch_id = 'bs'
    call AI%weights
    allocate(bs(max(1,AI%nlines),max(1,AI%ncolumns)))
    bs=AI%sp_weights
    nbsl=AI%nlines
    nbsc=AI%ncolumns


    !integer      :: OactivationFunction
    !integer      :: Onxl, Onxc, Onwh1l, Onwh1c, Onwsl, Onwsc, Onbh1l, Onbh1c, Onbsl, Onbsc 

   !with saved data allocate
   allocate(vh1(nwh1l,nxc))
   allocate(yh1(nwh1l,nxc))
   allocate(vs(nxl, nwh1c))
   allocate(ys(nxl, nwsc))

   call feedfoward(activationFunction, NXL, NXC, NWH1L, NWH1C,&
             NBH1L, NBH1C, NWSL, NWSC, NBSL,&
             NBSC, VALORES_NORMALIZADOS, WH1, BH1, WS, BS, VH1, YH1, VS, YS, VALORES_DENORMALIZADOS)

   
   deallocate(vh1)
   deallocate(yh1)
   deallocate(vs)
   deallocate(ys)
end subroutine pblh_rna_old

!contains

       SUBROUTINE FEEDFOWARD(activationFunction, NXL, NXC, NWH1L, NWH1C, NBH1L, NBH1C, NWSL, NWSC,&
                             NBSL, NBSC, VN, P1, B1, PES, BIS, V1, Y1, VSS, YSS, DENORMALIZADOS)
       INTEGER, INTENT(IN) :: activationFunction                       
       INTEGER, INTENT(IN) :: NXL, NXC 
       INTEGER, INTENT(IN) :: NWH1L, NWH1C 
       INTEGER, INTENT(IN) :: NBH1L, NBH1C 
       INTEGER, INTENT(IN) :: NWSL, NWSC
       INTEGER, INTENT(IN) :: NBSL, NBSC
       REAL, DIMENSION(NXL,NXC), INTENT(IN) :: VN
       REAL, DIMENSION(NWH1L,NWH1C), INTENT(IN) :: P1 
       REAL, DIMENSION(NBH1L,NBH1C), INTENT(IN) :: B1
       REAL, DIMENSION(NWSL, NWSC), INTENT(IN) :: PES
       REAL, DIMENSION(NBSL,NBSC), INTENT(IN) :: BIS
       REAL(8), DIMENSION(NWH1L,NXC), INTENT(OUT) :: V1
       REAL(8), DIMENSION(NWH1L,NXC), INTENT(OUT) :: Y1
       REAL(8), DIMENSION(NWSL,NXC), INTENT(OUT) :: VSS
       REAL(8), DIMENSION(NBSL,NXC), INTENT(OUT) :: YSS
       REAL(8), DIMENSION(NBSL,NXC), INTENT(OUT) :: DENORMALIZADOS
       
              
       INTEGER :: I,J
             
      !----------------------------------------------------------------------!
      ! INICIO DA REDE: FEEDFORWARD
      !----------------------------------------------------------------------!
      	!ACTIVATING HIDDEN LAYER 
	!Enver   open(6, file = 'result.out', STATUS = 'UNKNOWN', access = 'append')
	   open(666, file = 'result.out', STATUS = 'UNKNOWN', access = 'append')
	  
	   V1 = matmul(P1,VN) - B1
           !PRINT*,'V1=' 
	   !WRITE(6, *) V1
           !PRINT*,size(V1,dim=1),size(V1,dim=2)
           Y1 = 1.d0/(1.d0 + exp(-V1))	   
           !WRITE(6, *) Y1
           !PRINT*,'Y1='
           !PRINT*,'Y1=' 
           !PRINT*,size(Y1,dim=1),size(Y1,dim=2)
         
                                  
           VSS = matmul(PES,Y1) - BIS
           !PRINT*,'VSS='
	   !WRITE(6, *) VSS
           !PRINT*,size(VSS,dim=1),size(VSS,dim=2)
           YSS = 1.d0/(1.d0 + exp(-VSS))
           DENORMALIZADOS = YSS*2698.059178
           !PRINT*,'YSS='
           !Enver WRITE(6, *) DENORMALIZADOS           
           WRITE(666, *) DENORMALIZADOS           
	   !PRINT*,size(YSS,dim=1),size(YSS,dim=2)
           !Enver close(6)
           close(666)
                  	            
       END SUBROUTINE FEEDFOWARD 

       SUBROUTINE NORMALIZA (M, N, ET, NORMALIZADOS)
       INTEGER, INTENT(IN) :: M, N 
       REAL, DIMENSION(M,N), INTENT(IN) :: ET
       REAL, DIMENSION(M,N), INTENT(OUT) :: NORMALIZADOS
       INTEGER :: I
       REAL, PARAMETER :: th_max = 308.5164024	
       REAL, PARAMETER :: q_max = 43013845.97
       REAL, PARAMETER :: z_max = 1925.200431  
       REAL, PARAMETER :: u_max = -0.000283314
       REAL, PARAMETER :: v_max = 0.228466744
       REAL, PARAMETER :: t_max = 303.0733005
       REAL, PARAMETER :: pmid_max = 98501.42629
       REAL, PARAMETER :: kvf_max = 7.164266287
       REAL, PARAMETER :: shflx_max = 495.5330865
       REAL, PARAMETER :: taux_max = 0.30693858
       REAL, PARAMETER :: tauy_max = 0.53825646
       REAL, PARAMETER :: tsk_max = 301.922878
       REAL, PARAMETER :: qsfc_max = 0.032072326
       REAL, PARAMETER :: psomc_max = 1.06529817
       REAL, PARAMETER :: pblh_max = 2698.059178	
       REAL, PARAMETER :: th_min = 291.8784202
       REAL, PARAMETER :: q_min = 0
       REAL, PARAMETER :: z_min = 18.24950519
       REAL, PARAMETER :: u_min = -8.746260007 
       REAL, PARAMETER :: v_min = -4.937401516
       REAL, PARAMETER :: t_min = 285.837311
       REAL, PARAMETER :: pmid_min = 79210.10226
       REAL, PARAMETER :: kvf_min = 0.01
       REAL, PARAMETER :: shflx_min = -189.4587537
       REAL, PARAMETER :: taux_min = 3.80286E-05
       REAL, PARAMETER :: tauy_min = -0.002669614
       REAL, PARAMETER :: tsk_min =  294.446586
       REAL, PARAMETER :: qsfc_min = 0.016325435
       REAL, PARAMETER :: psomc_min = 1.001434148
       REAL, PARAMETER :: pblh_min = 18.61083813

       
       
           NORMALIZADOS(1,1) = (ET(1,1)-th_min)/(th_max-th_min)
           NORMALIZADOS(2,1) = (ET(2,1)-q_min)/(q_max-q_min)
           NORMALIZADOS(3,1) = (ET(3,1)-z_min)/(z_max-z_min)
           NORMALIZADOS(4,1) = (ET(4,1)-u_min)/(u_max-u_min)
           NORMALIZADOS(5,1) = (ET(5,1)-v_min)/(v_max-v_min)
           NORMALIZADOS(6,1) = (ET(6,1)-t_min)/(t_max-t_min)
           NORMALIZADOS(7,1) = (ET(7,1)-pmid_min)/(pmid_max-pmid_min)
           NORMALIZADOS(8,1) = (ET(8,1)-kvf_min)/(kvf_max-kvf_min)
           NORMALIZADOS(9,1) = (ET(9,1)-shflx_min)/(shflx_max-shflx_min)
           NORMALIZADOS(10,1) = (ET(10,1)-taux_min)/(taux_max-taux_min)
           NORMALIZADOS(11,1) = (ET(11,1)-tauy_min)/(tauy_max-tauy_min)
           NORMALIZADOS(12,1) = (ET(12,1)-tsk_min)/(tsk_max-tsk_min)
           NORMALIZADOS(13,1) = (ET(13,1)-qsfc_min)/(qsfc_max-qsfc_min)
           NORMALIZADOS(14,1) = (ET(14,1)-psomc_min)/(psomc_max-psomc_min)
           NORMALIZADOS(15,1) = (ET(15,1)-pblh_min)/(pblh_max-pblh_min)
          
       PRINT*, 'VALORES_NORMALIZADOS=' 
       PRINT*, NORMALIZADOS         
       END SUBROUTINE NORMALIZA


!end program test_rna
end module test_rna
