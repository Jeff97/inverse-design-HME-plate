      MODULE GLOBAL
C     ------------------------------------------------------------------
C     Global module to transfer SDVs from UEL to UVARM for visualization
C     on a dummy mesh.
C     ------------------------------------------------------------------
      IMPLICIT NONE
      SAVE

C     Parameters
      INTEGER, PARAMETER :: NELEMENT = 5000
      INTEGER, PARAMETER :: ELEMOFFSET = 100000
      INTEGER, PARAMETER :: nUvrm = 5

C     Global Array
C     GloUVar(ElementID, GaussPointID, VarID)
      REAL*8, ALLOCATABLE :: GloUVar(:,:,:)
      INTEGER :: err

      END MODULE GLOBAL

C     ==================================================================
C     Abaqus User Subroutine UEL for Magneto-Elastic Material
C     
C     Description:
C       Coupled implementation of deformation, pressure, and magnetic
C       potential.
C       Elements: C3D20R or similar 20-node hexahedra
C       Formulation: Mixed formulation (u-p-phi)
C
C     Developer: Zhanfeng Li
C     Last Update: 2026-01-11
C     ==================================================================
! ---------------------------------------------------------------------
!
!  8-node     8-----------7
!  brick     /|          /|       zeta
!           / |         / |       
!          5-----------6  |       |     eta
!          |  |        |  |       |   /
!          |  |        |  |       |  /
!          |  4--------|--3       | /
!          | /         | /        |/
!          |/          |/         O--------- xi
!          1-----------2        origin at cube center
!
!  20-node     8-----15----7       
!  brick      /|          /|       
!           16 |        14 |       
!           /  20       /  19      
!          5-----13----6   |       
!          |   |       |   |       
!          |   4----11-|---3       
!         17  /       18  /        
!          | 12        | 10        
!          |/          |/          
!          1-----9-----2           
! ---------------------------------------------------------------------
! User coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
C
      USE GLOBAL
      INCLUDE 'ABA_PARAM.INC'
C
C     ------------------------------------------------------------------
C     1. Argument Declarations
C     ------------------------------------------------------------------
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

C     ------------------------------------------------------------------
C     2. Local Variable Declarations
C     ------------------------------------------------------------------
C     
C     -- Geometry and Integration Points --
      INTEGER :: gp, ii, jj, kk, ll, pp
      INTEGER :: degree, nGP, EleType, nDim, nDofDisp, nDofp, nDofphi
      INTEGER :: nGP_p, nNode_p, degree_p
      INTEGER :: nGP_phi, nNode_phi, degree_phi
      
      DOUBLE PRECISION :: param(3), Jac, Jac_p, Jac_phi
      DOUBLE PRECISION :: GaussXc, GaussYc, GaussZc
      DOUBLE PRECISION :: dvol0, dvol, b1, b2, b3
      
C     -- Nodal Arrays and Shape Functions --
      DOUBLE PRECISION :: xNode(NNODE), yNode(NNODE), zNode(NNODE)
      DOUBLE PRECISION :: xNode_p(8), yNode_p(8), zNode_p(8)
      DOUBLE PRECISION :: xNode_phi(NNODE), yNode_phi(NNODE), 
     &                    zNode_phi(NNODE)
      
      DOUBLE PRECISION :: N(NNODE), dN_dX(NNODE), dN_dY(NNODE), 
     &                    dN_dZ(NNODE)
      DOUBLE PRECISION :: Np(8), dNp_dX(8), dNp_dY(8), dNp_dZ(8)
      DOUBLE PRECISION :: Nphi(NNODE), dNphi_dX(NNODE), dNphi_dY(NNODE),
     &                    dNphi_dZ(NNODE)
      DOUBLE PRECISION :: dN_dxx(NNODE), dN_dyy(NNODE), dN_dzz(NNODE)
      DOUBLE PRECISION :: dNphi_dxx(NNODE), dNphi_dyy(NNODE), 
     &                    dNphi_dzz(NNODE)
      
C     -- Material Properties --
      DOUBLE PRECISION :: C0, Nu, MValue, bValue, Theta_M, Theta_b, h
      DOUBLE PRECISION :: Gm, Kappa, Mu0, Epsilon, Pi
      
C     -- Field Variables --
      DOUBLE PRECISION :: Pressure, Phi, DetF, Ib_F
      DOUBLE PRECISION :: Cur_time, TotalT
      DOUBLE PRECISION :: Energy_Mech, Energy_Mag, Energy_Sum
      DOUBLE PRECISION :: vonMisesSS
      
C     -- Vectors and Tensors --
      DOUBLE PRECISION :: MRef(3), MCur(3), bExt(3), dphi(3),
     &                    G_phi_d(3)
      DOUBLE PRECISION :: F(3,3), FInv(3,3), b_F(3,3), Ae(3,3), G(3,3)
      DOUBLE PRECISION :: MyStrain(9), CauchySS(9)
      DOUBLE PRECISION :: Elasticity(9,9), Gc(3,9)
      DOUBLE PRECISION :: pCouple1(9,3), pCouple2(9,3)
      DOUBLE PRECISION :: Magneticity(3,3)
      DOUBLE PRECISION :: Gu_p(3,3), Gu_p_Tmp(3,3)
      
      DOUBLE PRECISION :: pElem(8), phiElem(NNODE)
      
C     -- Allocatable Arrays --
      DOUBLE PRECISION, ALLOCATABLE :: gpts1(:), gpts2(:), gpts3(:), 
     &                                 gwts(:)
      DOUBLE PRECISION, ALLOCATABLE :: gpts1_phi(:), gpts2_phi(:), 
     &                                 gpts3_phi(:), gwts_phi(:)
      DOUBLE PRECISION, ALLOCATABLE :: gpts1_p(:), gpts2_p(:), 
     &                                 gpts3_p(:), gwts_p(:)
      DOUBLE PRECISION, ALLOCATABLE :: dispElem(:), ResU(:), Res_phi(:)
      DOUBLE PRECISION, ALLOCATABLE :: ResP(:)
      
C     -- Stiffness Matrices --
      DOUBLE PRECISION, ALLOCATABLE :: Kuu(:,:), Kup(:,:), Kpu(:,:), 
     &                                 Kpp(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: K_u_phi(:,:), K_u_phi_Tmp(:,:), 
     &                                 K_phi_u(:,:), K_phi_phi(:,:)
      
C     ------------------------------------------------------------------
C     3. Initialization and Parameters
C     ------------------------------------------------------------------
      
      degree     = 2
      nGP        = 27
      EleType    = 8
      nDim       = 3
      nDofDisp   = NNODE * nDim
      nDofphi    = NNODE
      
      nGP_p      = 8
      nNode_p    = 8
      degree_p   = 1
      nDofp      = 8
      
      nGP_phi    = 27
      nNode_phi  = 20
      degree_phi = 2
      
      Pi = 3.14159265359

C     Allocate Local Arrays
      ALLOCATE(gpts1(nGP), gpts2(nGP), gpts3(nGP), gwts(nGP))
      ALLOCATE(gpts1_phi(nGP_phi), gpts2_phi(nGP_phi), 
     &         gpts3_phi(nGP_phi), gwts_phi(nGP_phi))
      ALLOCATE(gpts1_p(nGP_p), gpts2_p(nGP_p), 
     &         gpts3_p(nGP_p), gwts_p(nGP_p))
      
      ALLOCATE(dispElem(nDofDisp))
      ALLOCATE(ResU(nDofDisp))
      ALLOCATE(Res_phi(nDofphi))
      ALLOCATE(ResP(nDofp))
      
      ALLOCATE(Kuu(nDofDisp,nDofDisp))
      ALLOCATE(Kup(nDofDisp,nDofp))
      ALLOCATE(Kpu(nDofp,nDofDisp))
      ALLOCATE(Kpp(nDofp,nDofp))
      
      ALLOCATE(K_u_phi(nDofDisp,nDofphi))
      ALLOCATE(K_u_phi_Tmp(nDofDisp,nDofphi))
      ALLOCATE(K_phi_u(nDofphi,nDofDisp))
      ALLOCATE(K_phi_phi(nDofphi,nDofphi))
      
C     Allocate Global Variable for Visualization
      ALLOCATE(GloUVar(NELEMENT,nGP,nUvrm),stat=err)

C     ------------------------------------------------------------------
C     4. Material Properties Input (PROPS)
C     ------------------------------------------------------------------
C     PROPS(1) - C0: Material Parameter (Pa)
C     PROPS(2) - Nu: Poisson's ratio
C     PROPS(3) - MValue: Magnetization magnitude (kA/m)
C     PROPS(4) - bValue: Magnetic induction magnitude (mT)
C     PROPS(5) - Theta_M: Magnetization direction angle (Rad)
C     PROPS(6) - Theta_b: Magnetic induction direction angle (Rad)
C     PROPS(7) - h: Half thickness (m)
C     PROPS(8) - TotalT: Total loading time (s)

      C0      = PROPS(1)
      Nu      = PROPS(2)
      MValue  = PROPS(3) 
      bValue  = PROPS(4) 
      Theta_M = PROPS(5) 
      Theta_b = PROPS(6) 
      h       = PROPS(7) 
      TotalT  = PROPS(8)
      
      Cur_time = TIME(2)/TotalT
      
C     Derived Parameters
      Epsilon = bValue * MValue / (C0)
      Gm      = 2.0 * C0
      Kappa   = 2.0 * Gm * (1.0+Nu) / (3.0 * (1.0 - 2.0*Nu))
      Mu0     = 1.2566371 ! N (kA)^{-2}

C     Initial Print (Debug)
      IF (KSTEP .EQ. 1 .AND. KINC .EQ. 1) THEN
        IF (JELEM .EQ. 1) THEN
           WRITE(*,*) "Info: Coupled Magneto-Elastic UEL Initialized."
           WRITE(*,*) "C0:", C0, " MValue:", MValue, " bValue:", bValue
           WRITE(*,*) "Epsilon:", Epsilon
        END IF  
      END IF

C     External Fields (If coordinate is needed, can use GaussXc, GaussYc, GaussZc)
      bExt(1) = bValue * COS(Theta_b) * Cur_time
      bExt(2) = 0.0
      bExt(3) = bValue * SIN(Theta_b) * Cur_time

      MRef(1) = MValue * COS(Theta_M)
      MRef(2) = 0.0
      MRef(3) = MValue * SIN(Theta_M)

C     ------------------------------------------------------------------
C     5. Prepare Nodal Coordinates and DOF Variables
C     ------------------------------------------------------------------
      
      DO jj=1,NNODE
        xNode(jj) = COORDS(1,jj)
        yNode(jj) = COORDS(2,jj)
        zNode(jj) = COORDS(3,jj)
        
        xNode_phi(jj) = COORDS(1,jj)
        yNode_phi(jj) = COORDS(2,jj)
        zNode_phi(jj) = COORDS(3,jj)
      END DO

      DO jj=1,nNode_p
        xNode_p(jj) = COORDS(1,jj)
        yNode_p(jj) = COORDS(2,jj)
        zNode_p(jj) = COORDS(3,jj)
      END DO

C     Extract displacements (u), pressure (p), scalar potential (phi)
C     Structure of U:
C       Nodes 1-8:   u_x, u_y, u_z, phi, p (5 DOFs)
C       Nodes 9-20:  u_x, u_y, u_z, phi    (4 DOFs)

      DO ii=1, nNode_p ! 1 to 8
        dispElem(3*ii-2) = U(5*ii-4) ! x
        dispElem(3*ii-1) = U(5*ii-3) ! y
        dispElem(3*ii)   = U(5*ii-2) ! z
        phiElem(ii)      = U(5*ii-1) ! phi 
        pElem(ii)        = U(5*ii)   ! p
      END DO

      DO ii=1+nNode_p, NNODE ! 9 to 20
        dispElem(3*ii-2) = U(nDofp + 4*ii-3)
        dispElem(3*ii-1) = U(nDofp + 4*ii-2)
        dispElem(3*ii)   = U(nDofp + 4*ii-1)
        phiElem(ii)      = U(nDofp + 4*ii) 
      END DO

C     Initialize Global Matrices/Vectors
      Kuu = 0.0
      Kup = 0.0
      Kpu = 0.0
      K_u_phi = 0.0
      K_u_phi_Tmp = 0.0
      K_phi_u = 0.0
      Kpp = 0.0
      K_phi_phi = 0.0

      ResU = 0.0
      ResP = 0.0
      Res_phi = 0.0

C     Get Gauss Points
      call getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts)
      call getGaussPointsHexa(nGP_phi, gpts1_phi, gpts2_phi, 
     &                        gpts3_phi, gwts_phi)
      call getGaussPointsHexa(nGP_p, gpts1_p, gpts2_p, gpts3_p, gwts_p)

C     ------------------------------------------------------------------
C     6. Integration Loop
C     ------------------------------------------------------------------
      DO gp=1, nGP
        
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

C       -- Shape Functions --
        ! Displacement (u)
        call computeBasisFunctions3D(NNODE, EleType, degree, param,
     &    xNode, yNode, zNode, N, dN_dX, dN_dY, dN_dZ, Jac)
        
        IF(Jac < 0.0) THEN
          WRITE(*,*) "Warning: Negative Jacobian at Elem", JELEM
        END IF

        ! Potential (phi) - uses same shape functions as u usually
        call computeBasisFunctions3D(NNODE, EleType, degree, 
     &    param, xNode, yNode, zNode, 
     &    Nphi, dNphi_dX, dNphi_dY, dNphi_dZ, Jac_phi)

        ! Pressure (p)
        call computeBasisFunctions3D(nNode_p, EleType, degree_p, param,
     &    xNode_p, yNode_p, zNode_p, Np, dNp_dX, dNp_dY, dNp_dZ, Jac_p)

C       -- Interpolate Fields at Gauss Point --
        Pressure = 0.0
        DO ii=1, nNode_p
          Pressure = Pressure + pElem(ii) * Np(ii)
        END DO

        Phi = 0.0
        DO ii=1, NNODE
          Phi = Phi + phiElem(ii) * Nphi(ii)
        END DO

C       -- Deformation Gradient F --
        F = 0.0
        F(1,1) = 1.0; F(2,2) = 1.0; F(3,3) = 1.0

        DO ii=1, NNODE
          F(1,1) = F(1,1) + dispElem(3*ii-2) * dN_dX(ii)
          F(1,2) = F(1,2) + dispElem(3*ii-2) * dN_dY(ii)
          F(1,3) = F(1,3) + dispElem(3*ii-2) * dN_dZ(ii)
          F(2,1) = F(2,1) + dispElem(3*ii-1) * dN_dX(ii)
          F(2,2) = F(2,2) + dispElem(3*ii-1) * dN_dY(ii)
          F(2,3) = F(2,3) + dispElem(3*ii-1) * dN_dZ(ii)
          F(3,1) = F(3,1) + dispElem(3*ii)   * dN_dX(ii)
          F(3,2) = F(3,2) + dispElem(3*ii)   * dN_dY(ii)
          F(3,3) = F(3,3) + dispElem(3*ii)   * dN_dZ(ii)
        END DO

        call DetMatrix(F, DetF)
        call InvMatrix(FInv, F, DetF)

C       -- Current Magnetization MCur = J^(-1) * F . MRef --
        MCur(1) = (F(1,1)*MRef(1)+F(1,2)*MRef(2)+F(1,3)*MRef(3))/DetF
        MCur(2) = (F(2,1)*MRef(1)+F(2,2)*MRef(2)+F(2,3)*MRef(3))/DetF
        MCur(3) = (F(3,1)*MRef(1)+F(3,2)*MRef(2)+F(3,3)*MRef(3))/DetF

C       -- Left Cauchy-Green b = F.F^T --
        b_F = 0.0
        DO ii=1, 3
          DO jj=1, 3
            DO KK=1, 3
              b_F(ii,jj) = b_F(ii,jj) + F(ii,KK)*F(jj,KK)
            END DO
          END DO
        END DO
        Ib_F = b_F(1,1) + b_F(2,2) + b_F(3,3)

C       -- Gradients in Current Configuration --

        DO ii=1, NNODE
          ! For displacement N
          dN_dxx(ii) = dN_dX(ii)*FInv(1,1) + dN_dY(ii)*FInv(2,1) + 
     &                 dN_dZ(ii)*FInv(3,1)
          dN_dyy(ii) = dN_dX(ii)*FInv(1,2) + dN_dY(ii)*FInv(2,2) + 
     &                 dN_dZ(ii)*FInv(3,2)
          dN_dzz(ii) = dN_dX(ii)*FInv(1,3) + dN_dY(ii)*FInv(2,3) + 
     &                 dN_dZ(ii)*FInv(3,3)
          
          ! For potential Nphi
          dNphi_dxx(ii) = dNphi_dX(ii)*FInv(1,1) + 
     &                    dNphi_dY(ii)*FInv(2,1) + 
     &                    dNphi_dZ(ii)*FInv(3,1)
          dNphi_dyy(ii) = dNphi_dX(ii)*FInv(1,2) + 
     &                    dNphi_dY(ii)*FInv(2,2) + 
     &                    dNphi_dZ(ii)*FInv(3,2)
          dNphi_dzz(ii) = dNphi_dX(ii)*FInv(1,3) + 
     &                    dNphi_dY(ii)*FInv(2,3) + 
     &                    dNphi_dZ(ii)*FInv(3,3)
        END DO

        dphi = 0.0
        DO ii=1, NNODE
          dphi(1) = dphi(1) + dNphi_dxx(ii)*phiElem(ii)
          dphi(2) = dphi(2) + dNphi_dyy(ii)*phiElem(ii)
          dphi(3) = dphi(3) + dNphi_dzz(ii)*phiElem(ii)
        END DO

C       -- Constitutive Tensors --
        Elasticity = 0.0 
        call ElasticityTensor(Elasticity, Pressure, Gm, Kappa,
     &     DetF, Ib_F, b_F, Mu0, dphi, MCur)

        pCouple1 = 0.0 
        call CouplingTensor1(pCouple1, Mu0, dphi, DetF, F, MCur)

        pCouple2 = 0.0 
        call CouplingTensor2(pCouple2, Mu0, dphi, DetF, F, MCur)

        Magneticity = 0.0 
        call MagneticityTensor(Magneticity, Mu0)

C       -- Cauchy stress calculation --
        CauchySS = 0.0
        call cal_CauchyStress(CauchySS, Pressure, Gm, Kappa,
     &     DetF, Ib_F, b_F, dphi, F, MCur, bExt)

        vonMisesSS = Sqrt(0.5*((CauchySS(1)-CauchySS(5))**2
     &    +(CauchySS(5)-CauchySS(9))**2+(CauchySS(9)-CauchySS(1))**2)
     &    +3.0*(CauchySS(2)**2 + CauchySS(6)**2 + CauchySS(7)**2))

C       -- Volume Elements --
        dvol0 = gwts(gp) * Jac 
        dvol  = gwts(gp) * Jac * DetF

C       -- Energy Calculation --
        Energy_Mech = ( C0*(DetF**(-2.0/3.0)*Ib_F-3.0)
     &         +Pressure*(DetF-1.0) -Pressure*Pressure/(2.0*Kappa) )
        Energy_Mag = -MCur(1)*bExt(1) -MCur(2)*bExt(2) -MCur(3)*bExt(3)
        Energy_Sum = Energy_Mech + Energy_Mag

C       -- TODO: SDVs for Visualization --
        GloUVar(JELEM,gp,1)  = vonMisesSS
        GloUVar(JELEM,gp,2)  = Pressure
        GloUVar(JELEM,gp,3)  = Phi
      !   GloUVar(JELEM,gp,4)  = bExt(1)
      !   GloUVar(JELEM,gp,5)  = bExt(2)
      !   GloUVar(JELEM,gp,6)  = bExt(3)
      !   GloUVar(JELEM,gp,7)  = MRef(1)
      !   GloUVar(JELEM,gp,8)  = MRef(2)
      !   GloUVar(JELEM,gp,9)  = MRef(3)
        GloUVar(JELEM,gp,4) = Energy_Mech
        GloUVar(JELEM,gp,5) = Energy_Mag

C       -- Assembly Loop (Over Nodes) --
        Gc = 0.0
        DO ii=1, NNODE
          b1 = dN_dxx(ii)*dvol
          b2 = dN_dyy(ii)*dvol
          b3 = dN_dzz(ii)*dvol
          
          ! Elasticity Contribution Gc
          DO pp=1, 9
            Gc(1,pp) = b1*Elasticity(1,pp) + b2*Elasticity(4,pp) + 
     &                 b3*Elasticity(7,pp)
            Gc(2,pp) = b1*Elasticity(2,pp) + b2*Elasticity(5,pp) + 
     &                 b3*Elasticity(8,pp)
            Gc(3,pp) = b1*Elasticity(3,pp) + b2*Elasticity(6,pp) + 
     &                 b3*Elasticity(9,pp)
          END DO

          ! Coupling Contribution Gu_p
          DO pp=1, 3
            Gu_p(1,pp) = b1*pCouple1(1,pp) + b2*pCouple1(4,pp) + 
     &                   b3*pCouple1(7,pp)
            Gu_p(2,pp) = b1*pCouple1(2,pp) + b2*pCouple1(5,pp) + 
     &                   b3*pCouple1(8,pp)
            Gu_p(3,pp) = b1*pCouple1(3,pp) + b2*pCouple1(6,pp) + 
     &                   b3*pCouple1(9,pp)
          END DO

          ! Coupling Contribution Gu_p_Tmp
          DO pp=1, 3
            Gu_p_Tmp(1,pp) = b1*pCouple2(1,pp) + b2*pCouple2(4,pp) + 
     &                       b3*pCouple2(7,pp)
            Gu_p_Tmp(2,pp) = b1*pCouple2(2,pp) + b2*pCouple2(5,pp) + 
     &                       b3*pCouple2(8,pp)
            Gu_p_Tmp(3,pp) = b1*pCouple2(3,pp) + b2*pCouple2(6,pp) + 
     &                       b3*pCouple2(9,pp)
          END DO
          
          ! Assembly: Kuu (Stiffness)
          DO jj=1, NNODE
             Kuu(3*ii-2, 3*jj-2) = Kuu(3*ii-2, 3*jj-2) 
     &         + Gc(1,1)*dN_dxx(jj) + Gc(1,4)*dN_dyy(jj) 
     &         + Gc(1,7)*dN_dzz(jj)
             Kuu(3*ii-2, 3*jj-1) = Kuu(3*ii-2, 3*jj-1)
     &         + Gc(1,2)*dN_dxx(jj) + Gc(1,5)*dN_dyy(jj) 
     &         + Gc(1,8)*dN_dzz(jj)
             Kuu(3*ii-2, 3*jj)   = Kuu(3*ii-2, 3*jj)  
     &         + Gc(1,3)*dN_dxx(jj) + Gc(1,6)*dN_dyy(jj) 
     &         + Gc(1,9)*dN_dzz(jj)

             Kuu(3*ii-1, 3*jj-2) = Kuu(3*ii-1, 3*jj-2)
     &         + Gc(2,1)*dN_dxx(jj) + Gc(2,4)*dN_dyy(jj) 
     &         + Gc(2,7)*dN_dzz(jj)
             Kuu(3*ii-1, 3*jj-1) = Kuu(3*ii-1, 3*jj-1)
     &         + Gc(2,2)*dN_dxx(jj) + Gc(2,5)*dN_dyy(jj) 
     &         + Gc(2,8)*dN_dzz(jj)
             Kuu(3*ii-1, 3*jj)   = Kuu(3*ii-1, 3*jj)  
     &         + Gc(2,3)*dN_dxx(jj) + Gc(2,6)*dN_dyy(jj) 
     &         + Gc(2,9)*dN_dzz(jj)

             Kuu(3*ii, 3*jj-2) = Kuu(3*ii, 3*jj-2)
     &         + Gc(3,1)*dN_dxx(jj) + Gc(3,4)*dN_dyy(jj) 
     &         + Gc(3,7)*dN_dzz(jj)
             Kuu(3*ii, 3*jj-1) = Kuu(3*ii, 3*jj-1)
     &         + Gc(3,2)*dN_dxx(jj) + Gc(3,5)*dN_dyy(jj) 
     &         + Gc(3,8)*dN_dzz(jj)
             Kuu(3*ii, 3*jj)   = Kuu(3*ii, 3*jj)  
     &         + Gc(3,3)*dN_dxx(jj) + Gc(3,6)*dN_dyy(jj) 
     &         + Gc(3,9)*dN_dzz(jj)
          END DO
          
          ! Assembly: K_u_phi
          DO jj=1, NNODE
            K_u_phi(3*ii-2,jj) = K_u_phi(3*ii-2,jj)      
     &        + Gu_p(1,1)*dNphi_dxx(jj) + Gu_p(1,2)*dNphi_dyy(jj) 
     &        + Gu_p(1,3)*dNphi_dzz(jj)
            K_u_phi(3*ii-1,jj) = K_u_phi(3*ii-1,jj)       
     &        + Gu_p(2,1)*dNphi_dxx(jj) + Gu_p(2,2)*dNphi_dyy(jj) 
     &        + Gu_p(2,3)*dNphi_dzz(jj)
            K_u_phi(3*ii,jj)   = K_u_phi(3*ii,jj)       
     &        + Gu_p(3,1)*dNphi_dxx(jj) + Gu_p(3,2)*dNphi_dyy(jj) 
     &        + Gu_p(3,3)*dNphi_dzz(jj)
          END DO

          ! Assembly: K_u_phi_Tmp (Used for K_phi_u later)
          DO jj=1, NNODE
            K_u_phi_Tmp(3*ii-2,jj) = K_u_phi_Tmp(3*ii-2,jj)      
     &        + Gu_p_Tmp(1,1)*dNphi_dxx(jj) + 
     &          Gu_p_Tmp(1,2)*dNphi_dyy(jj) + 
     &          Gu_p_Tmp(1,3)*dNphi_dzz(jj)
            K_u_phi_Tmp(3*ii-1,jj) = K_u_phi_Tmp(3*ii-1,jj)       
     &        + Gu_p_Tmp(2,1)*dNphi_dxx(jj) + 
     &          Gu_p_Tmp(2,2)*dNphi_dyy(jj) + 
     &          Gu_p_Tmp(2,3)*dNphi_dzz(jj)
            K_u_phi_Tmp(3*ii,jj)   = K_u_phi_Tmp(3*ii,jj)       
     &        + Gu_p_Tmp(3,1)*dNphi_dxx(jj) + 
     &          Gu_p_Tmp(3,2)*dNphi_dyy(jj) + 
     &          Gu_p_Tmp(3,3)*dNphi_dzz(jj)
          END DO

          ! Assembly: Kup
          DO jj=1, nNode_p
            Kup(3*ii-2,jj) = Kup(3*ii-2,jj) + b1*Np(jj)
            Kup(3*ii-1,jj) = Kup(3*ii-1,jj) + b2*Np(jj)
            Kup(3*ii,jj)   = Kup(3*ii,jj)   + b3*Np(jj)
          END DO

          ! Assembly: K_phi_phi
          DO pp=1, 3
            G_phi_d(pp) = dNphi_dxx(ii)*dvol*Magneticity(1,pp) + 
     &                    dNphi_dyy(ii)*dvol*Magneticity(2,pp) + 
     &                    dNphi_dzz(ii)*dvol*Magneticity(3,pp) 
          END DO

          DO jj=1, NNODE
             K_phi_phi(ii,jj) = K_phi_phi(ii,jj) + 
     &         G_phi_d(1)*dNphi_dxx(jj) + 
     &         G_phi_d(2)*dNphi_dyy(jj) + 
     &         G_phi_d(3)*dNphi_dzz(jj)
          END DO

          ! Residuals: ResU
          ResU(3*ii-2) = ResU(3*ii-2) - b1*CauchySS(1) - 
     &                   b2*CauchySS(2) - b3*CauchySS(3)
          ResU(3*ii-1) = ResU(3*ii-1) - b1*CauchySS(4) - 
     &                   b2*CauchySS(5) - b3*CauchySS(6)
          ResU(3*ii)   = ResU(3*ii)   - b1*CauchySS(7) - 
     &                   b2*CauchySS(8) - b3*CauchySS(9)

          ! Residuals: Res_phi
          Res_phi(ii) = Res_phi(ii) 
     &       - dNphi_dxx(ii)*dvol*Mu0*(MCur(1)-dphi(1))
     &       - dNphi_dyy(ii)*dvol*Mu0*(MCur(2)-dphi(2))
     &       - dNphi_dzz(ii)*dvol*Mu0*(MCur(3)-dphi(3))

        END DO ! End Node Loop

        ! Residuals: ResP
        DO ii=1, nNode_p
          ResP(ii) = ResP(ii) - Np(ii)*(DetF -1.0 -Pressure/kappa)
     &               *dvol0
        END DO

        ! Assembly: Kpp
        DO ii=1, nNode_p
          DO jj=1, nNode_p
            Kpp(ii,jj) = Kpp(ii,jj) - 
     &                   Np(ii)*Np(jj)/Kappa*dvol0
          END DO
        END DO

      END DO ! End Gauss Point Loop

C     ------------------------------------------------------------------
C     7. Matrix Transpositions and Final Assembly
C     ------------------------------------------------------------------
      
      ! K_phi_u = Transpose(K_u_phi_Tmp)
      DO ii=1, NNODE
        DO jj=1, NNODE
          K_phi_u(jj,3*ii-2) = K_u_phi_Tmp(3*ii-2,jj)      
          K_phi_u(jj,3*ii-1) = K_u_phi_Tmp(3*ii-1,jj)      
          K_phi_u(jj,3*ii)   = K_u_phi_Tmp(3*ii,jj)  
        END DO
      END DO

      ! Kpu = Transpose(Kup)
      DO ii=1, NNODE
        DO jj=1, nNode_p
          Kpu(jj,3*ii-2) = Kup(3*ii-2,jj) 
          Kpu(jj,3*ii-1) = Kup(3*ii-1,jj) 
          Kpu(jj,3*ii)   = Kup(3*ii,jj) 
        END DO
      END DO

C     Assemble AMATRX
C     Structure:
C     u-nodes 1-8: [Kuu  K_u_phi Kup]
C     phi-nodes 1-8:[K_phi_u K_phi_phi 0]
C     p-nodes 1-8: [Kpu  0         Kpp]
C     ... similar for 9-20 ...

      ! Part A: Nodes 1-8 (Heavy coupling u, phi, p)
      DO ii=1, nDofp
        DO jj=1, nDofp
          ! Disp-Disp
          AMATRX(5*ii-4,5*jj-4) = Kuu(3*ii-2,3*jj-2)
          AMATRX(5*ii-4,5*jj-3) = Kuu(3*ii-2,3*jj-1)
          AMATRX(5*ii-4,5*jj-2) = Kuu(3*ii-2,3*jj)
          ! Disp-Phi
          AMATRX(5*ii-4,5*jj-1) = K_u_phi(3*ii-2,jj)
          ! Disp-P
          AMATRX(5*ii-4,5*jj)   = Kup(3*ii-2,jj)

          AMATRX(5*ii-3,5*jj-4) = Kuu(3*ii-1,3*jj-2)
          AMATRX(5*ii-3,5*jj-3) = Kuu(3*ii-1,3*jj-1)
          AMATRX(5*ii-3,5*jj-2) = Kuu(3*ii-1,3*jj)
          AMATRX(5*ii-3,5*jj-1) = K_u_phi(3*ii-1,jj)
          AMATRX(5*ii-3,5*jj)   = Kup(3*ii-1,jj)

          AMATRX(5*ii-2,5*jj-4) = Kuu(3*ii,3*jj-2)
          AMATRX(5*ii-2,5*jj-3) = Kuu(3*ii,3*jj-1)
          AMATRX(5*ii-2,5*jj-2) = Kuu(3*ii,3*jj)
          AMATRX(5*ii-2,5*jj-1) = K_u_phi(3*ii,jj)
          AMATRX(5*ii-2,5*jj)   = Kup(3*ii,jj)

          ! Phi-Disp
          AMATRX(5*ii-1,5*jj-4) = K_phi_u(ii,3*jj-2)
          AMATRX(5*ii-1,5*jj-3) = K_phi_u(ii,3*jj-1)
          AMATRX(5*ii-1,5*jj-2) = K_phi_u(ii,3*jj)
          ! Phi-Phi
          AMATRX(5*ii-1,5*jj-1) = K_phi_phi(ii,jj)
          ! Phi-P
          AMATRX(5*ii-1,5*jj)   = 0.0

          ! P-Disp
          AMATRX(5*ii,5*jj-4) = Kpu(ii,3*jj-2)
          AMATRX(5*ii,5*jj-3) = Kpu(ii,3*jj-1)
          AMATRX(5*ii,5*jj-2) = Kpu(ii,3*jj)
          ! P-Phi
          AMATRX(5*ii,5*jj-1) = 0.0
          ! P-P
          AMATRX(5*ii,5*jj)   = Kpp(ii,jj)
        END DO 
      END DO

      ! Part B: Nodes 1-8 Interactions with Nodes 9-20
      DO ii=1, nDofp
        DO jj=1+nDofp, NNODE
          ! Columns for nodes 9-20 have 4 DOFs (u_x, u_y, u_z, phi)
          AMATRX(5*ii-4,nDofp+4*jj-3) = Kuu(3*ii-2,3*jj-2)
          AMATRX(5*ii-4,nDofp+4*jj-2) = Kuu(3*ii-2,3*jj-1)
          AMATRX(5*ii-4,nDofp+4*jj-1) = Kuu(3*ii-2,3*jj)
          AMATRX(5*ii-4,nDofp+4*jj)   = K_u_phi(3*ii-2,jj)

          AMATRX(5*ii-3,nDofp+4*jj-3) = Kuu(3*ii-1,3*jj-2)
          AMATRX(5*ii-3,nDofp+4*jj-2) = Kuu(3*ii-1,3*jj-1)
          AMATRX(5*ii-3,nDofp+4*jj-1) = Kuu(3*ii-1,3*jj)
          AMATRX(5*ii-3,nDofp+4*jj)   = K_u_phi(3*ii-1,jj)

          AMATRX(5*ii-2,nDofp+4*jj-3) = Kuu(3*ii,3*jj-2)
          AMATRX(5*ii-2,nDofp+4*jj-2) = Kuu(3*ii,3*jj-1)
          AMATRX(5*ii-2,nDofp+4*jj-1) = Kuu(3*ii,3*jj)
          AMATRX(5*ii-2,nDofp+4*jj)   = K_u_phi(3*ii,jj)

          AMATRX(5*ii-1,nDofp+4*jj-3) = K_phi_u(ii,3*jj-2)
          AMATRX(5*ii-1,nDofp+4*jj-2) = K_phi_u(ii,3*jj-1)
          AMATRX(5*ii-1,nDofp+4*jj-1) = K_phi_u(ii,3*jj)
          AMATRX(5*ii-1,nDofp+4*jj)   = K_phi_phi(ii,jj)

          AMATRX(5*ii,nDofp+4*jj-3)   = Kpu(ii,3*jj-2)
          AMATRX(5*ii,nDofp+4*jj-2)   = Kpu(ii,3*jj-1)
          AMATRX(5*ii,nDofp+4*jj-1)   = Kpu(ii,3*jj)
          AMATRX(5*ii,nDofp+4*jj)     = 0.0
        END DO 
      END DO

      ! Part C: Nodes 9-20 Interactions with Nodes 1-8
      DO ii=1+nDofp, NNODE
        DO jj=1, nDofp
          AMATRX(nDofp+4*ii-3,5*jj-4) = Kuu(3*ii-2,3*jj-2)
          AMATRX(nDofp+4*ii-3,5*jj-3) = Kuu(3*ii-2,3*jj-1)
          AMATRX(nDofp+4*ii-3,5*jj-2) = Kuu(3*ii-2,3*jj)
          AMATRX(nDofp+4*ii-3,5*jj-1) = K_u_phi(3*ii-2,jj)
          AMATRX(nDofp+4*ii-3,5*jj)   = Kup(3*ii-2,jj)

          AMATRX(nDofp+4*ii-2,5*jj-4) = Kuu(3*ii-1,3*jj-2)
          AMATRX(nDofp+4*ii-2,5*jj-3) = Kuu(3*ii-1,3*jj-1)
          AMATRX(nDofp+4*ii-2,5*jj-2) = Kuu(3*ii-1,3*jj)
          AMATRX(nDofp+4*ii-2,5*jj-1) = K_u_phi(3*ii-1,jj)
          AMATRX(nDofp+4*ii-2,5*jj)   = Kup(3*ii-1,jj)

          AMATRX(nDofp+4*ii-1,5*jj-4) = Kuu(3*ii,3*jj-2)
          AMATRX(nDofp+4*ii-1,5*jj-3) = Kuu(3*ii,3*jj-1)
          AMATRX(nDofp+4*ii-1,5*jj-2) = Kuu(3*ii,3*jj)
          AMATRX(nDofp+4*ii-1,5*jj-1) = K_u_phi(3*ii,jj)
          AMATRX(nDofp+4*ii-1,5*jj)   = Kup(3*ii,jj)

          AMATRX(nDofp+4*ii,5*jj-4)   = K_phi_u(ii,3*jj-2)
          AMATRX(nDofp+4*ii,5*jj-3)   = K_phi_u(ii,3*jj-1)
          AMATRX(nDofp+4*ii,5*jj-2)   = K_phi_u(ii,3*jj)
          AMATRX(nDofp+4*ii,5*jj-1)   = K_phi_phi(ii,jj)
          AMATRX(nDofp+4*ii,5*jj)     = 0.0
        END DO 
      END DO

      ! Part D: Nodes 9-20 Interactions with Nodes 9-20
      DO ii=1+nDofp, NNODE
        DO jj=1+nDofp, NNODE
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-3) = Kuu(3*ii-2,3*jj-2)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-2) = Kuu(3*ii-2,3*jj-1)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-1) = Kuu(3*ii-2,3*jj)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj)   = K_u_phi(3*ii-2,jj)

          AMATRX(nDofp+4*ii-2,nDofp+4*jj-3) = Kuu(3*ii-1,3*jj-2)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj-2) = Kuu(3*ii-1,3*jj-1)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj-1) = Kuu(3*ii-1,3*jj)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj)   = K_u_phi(3*ii-1,jj)

          AMATRX(nDofp+4*ii-1,nDofp+4*jj-3) = Kuu(3*ii,3*jj-2)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj-2) = Kuu(3*ii,3*jj-1)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj-1) = Kuu(3*ii,3*jj)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj)   = K_u_phi(3*ii,jj)

          AMATRX(nDofp+4*ii,nDofp+4*jj-3)   = K_phi_u(ii,3*jj-2)
          AMATRX(nDofp+4*ii,nDofp+4*jj-2)   = K_phi_u(ii,3*jj-1)
          AMATRX(nDofp+4*ii,nDofp+4*jj-1)   = K_phi_u(ii,3*jj)
          AMATRX(nDofp+4*ii,nDofp+4*jj)     = K_phi_phi(ii,jj)
        END DO 
      END DO

C     Assemble RHS (Residuals)
C     Nodes 1-8
      DO ii=1, nDofp
        RHS(5*ii-4,1) =  ResU(3*ii-2) 
        RHS(5*ii-3,1) =  ResU(3*ii-1) 
        RHS(5*ii-2,1) =  ResU(3*ii) 
        RHS(5*ii-1,1) =  Res_phi(ii)
        RHS(5*ii  ,1) =  ResP(ii) 
      END DO

C     Nodes 9-20
      DO ii=1+nDofp, NNODE
        RHS(nDofp+4*ii-3,1) = ResU(3*ii-2)
        RHS(nDofp+4*ii-2,1) = ResU(3*ii-1)
        RHS(nDofp+4*ii-1,1) = ResU(3*ii) 
        RHS(nDofp+4*ii  ,1) = Res_phi(ii)
      END DO

      DEALLOCATE(gpts1, gpts2, gpts3, gwts)
      DEALLOCATE(gpts1_phi, gpts2_phi, gpts3_phi, gwts_phi)
      DEALLOCATE(gpts1_p, gpts2_p, gpts3_p, gwts_p)
      DEALLOCATE(dispElem, ResU, Res_phi, ResP)
      DEALLOCATE(Kuu, Kup, Kpu, Kpp)
      DEALLOCATE(K_u_phi, K_u_phi_Tmp, K_phi_u, K_phi_phi)

      RETURN
      END SUBROUTINE UEL

C     ==================================================================
C     Helper Subroutines
C     ==================================================================

      SUBROUTINE MagneticityTensor(ElecTensor, Mu0)
      IMPLICIT NONE
      DOUBLE PRECISION :: ElecTensor(3,3), Mu0
      INTEGER :: ii, jj

      DO ii=1, 3
        DO jj=1, 3
          IF (ii==jj) THEN
           ElecTensor(ii,jj) = ElecTensor(ii,jj) - Mu0
          END IF
        END DO
      END DO

      END SUBROUTINE MagneticityTensor
!=======================================================================
      SUBROUTINE CouplingTensor1(CoupTensor, Mu0, dphi, DetF, F, MCur)
      IMPLICIT NONE
      DOUBLE PRECISION :: CoupTensor(9,3), Mu0, dphi(3), DetF, F(3,3),
     &                    MCur(3)
      INTEGER :: ii, jj, kk

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (ii==jj) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          - Mu0*dphi(kk)
            END IF
          END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (ii==kk) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          + Mu0*dphi(jj)
            END IF
          END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (jj==kk) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          + Mu0*dphi(ii)
            END IF
          END DO
        END DO
      END DO

      END SUBROUTINE CouplingTensor1
!=======================================================================
      SUBROUTINE CouplingTensor2(CoupTensor, Mu0, dphi, DetF, F, MCur)
      IMPLICIT NONE
      DOUBLE PRECISION :: CoupTensor(9,3), Mu0, dphi(3), DetF, F(3,3),
     &                    MCur(3)
      INTEGER :: ii, jj, kk

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (ii==jj) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          - Mu0*dphi(kk)
            END IF
          END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (ii==kk) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          + Mu0*dphi(jj)
            END IF
          END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
          DO kk=1, 3
            IF (jj==kk) THEN
             CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
     &          + Mu0*dphi(ii)
            END IF
          END DO
        END DO
      END DO

      END SUBROUTINE CouplingTensor2
!=======================================================================
      SUBROUTINE ElasticityTensor(ElsTensor, intPressure, 
     1     Gm, Kappa, DetF, Ib, b, Mu0, dphi, MCur)
      IMPLICIT NONE
      DOUBLE PRECISION :: ElsTensor(9,9), b(3,3), dphi(3), MCur(3)
      INTEGER :: ii, jj, kk, ll
      DOUBLE PRECISION :: Gm, Kappa, Mu0
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi

      ! mechanical terms
      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0D0/3.0D0*Gm*DetF**(-5.0D0/3.0D0)*b(kk,ll)
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj .AND. ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0D0/3.0D0*Gm*DetF**(-5.0D0/3.0D0)*(-1.0D0/3.0D0)*Ib
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (ii==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Gm*DetF**(-5.0D0/3.0D0)*b(jj,ll)
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Gm*DetF**(-5.0D0/3.0D0)*(-2.0D0/3.0D0)*b(ii,jj)
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (jj==kk .AND. ll==ii) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Gm*DetF**(-5.0D0/3.0D0)*(1.0D0/3.0D0)*Ib
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (jj==ii .AND. ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + intPressure
               END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
               IF (ll==ii .AND. jj==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            - intPressure
               END IF
              END DO
           END DO
        END DO
      END DO

      ! coupling terms
      Sq_dphi = dphi(1)**2 + dphi(2)**2 + dphi(3)**2

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (jj==kk .AND. ii==ll) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               + 0.5D0*Mu0*Sq_dphi
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (ii==jj .AND. ll==kk) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               - 0.5D0*Mu0*Sq_dphi
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (jj==ii) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu0*dphi(kk)*dphi(ll)
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (ii==ll) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            - Mu0*dphi(jj)*dphi(kk)
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu0*dphi(ii)*dphi(jj)
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (jj==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            - Mu0*dphi(ii)*dphi(ll)
                END IF
              END DO
           END DO
        END DO
      END DO

      DO ii=1, 3
        DO jj=1, 3
           DO kk=1, 3
              DO ll=1, 3
                IF (jj==ll) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            - Mu0*dphi(ii)*dphi(kk)
                END IF
              END DO
           END DO
        END DO
      END DO

      END SUBROUTINE ElasticityTensor
!=======================================================================
      SUBROUTINE cal_CauchyStress(CauchySS, intPressure, 
     1     Gm, Kappa, DetF, Ib, b, dphi, F, MCur, b_apply)
      IMPLICIT NONE
      DOUBLE PRECISION :: b(3,3), CauchySS(9), dphi(3), F(3,3), MCur(3)
      DOUBLE PRECISION :: Gm, Kappa, b_apply(3)
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi, Coef1

      Sq_dphi = dphi(1)**2 + dphi(2)**2 + dphi(3)**2
      Coef1 = Gm*DetF**(-5.0D0/3.0D0)

      CauchySS(1) = Coef1*(b(1,1) - Ib/3.0D0) + intPressure
     &            - MCur(1)*b_apply(1)
      CauchySS(2) = Coef1*b(2,1) - MCur(2)*b_apply(1)
      CauchySS(3) = Coef1*b(3,1)- MCur(3)*b_apply(1)

      CauchySS(4) = Coef1*b(1,2)- MCur(1)*b_apply(2)
      CauchySS(5) = Coef1*(b(2,2) - Ib/3.0D0) + intPressure
     &            - MCur(2)*b_apply(2)
      CauchySS(6) = Coef1*b(3,2)- MCur(3)*b_apply(2)

      CauchySS(7) = Coef1*b(1,3)- MCur(1)*b_apply(3)
      CauchySS(8) = Coef1*b(2,3)- MCur(2)*b_apply(3)
      CauchySS(9) = Coef1*(b(3,3) - Ib/3.0D0) + intPressure
     &            - MCur(3)*b_apply(3)

      END SUBROUTINE cal_CauchyStress
!=======================================================================
      SUBROUTINE UVARM(uvar,direct,t,time,dtime,cmname,orname,
     1 nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     2 jmac,jmatyp,matlayo,laccfla)

      USE GLOBAL
      INCLUDE 'ABA_PARAM.INC'

      INTEGER :: nuvarm, noel, npt, layer, kspt, kstep, kinc, ndi, nshr
      INTEGER :: jmac(*), jmatyp(*), matlayo, laccfla
      DOUBLE PRECISION :: uvar(nuvarm), direct(3,3), t(3,3), time(2)
      DOUBLE PRECISION :: dtime, coord(*)
      CHARACTER*80 cmname,orname
      INTEGER :: NELEMAN, i

      NELEMAN = NOEL - ELEMOFFSET

      IF (NELEMAN > 0 .AND. NELEMAN <= NELEMENT) THEN
         ! Check bounds to prevent crash if mesh is not set up correctly
         DO i=1,nUvrm
             uvar(i) = GloUVar(NELEMAN,NPT,i)
         END DO
      ELSE
         uvar = 0.0D0
      END IF

      RETURN
      END SUBROUTINE UVARM
!=======================================================================
      SUBROUTINE GetGaussPoints1D(nGP, GaussPoints, GaussWeights)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nGP
      DOUBLE PRECISION :: GaussPoints(nGP), GaussWeights(nGP)

      SELECT CASE(nGP)
        CASE(1)  
          GaussPoints(1) = 0.0
          GaussWeights(1) = 2.0
        CASE(2)  
          GaussPoints(1) = -0.577350269189626
          GaussWeights(1) = 1.0
          GaussPoints(2) =  0.577350269189626
          GaussWeights(2) = 1.0
        CASE(3)  
          GaussPoints(1) = -0.774596669241483
          GaussWeights(1) = 0.555555555555556
          GaussPoints(2) =  0.0
          GaussWeights(2) = 0.888888888888889
          GaussPoints(3) =  0.774596669241483
          GaussWeights(3) = 0.555555555555556
        CASE(4) 
          GaussPoints(1) = -0.861136311594953
          GaussWeights(1) = 0.347854845137454
          GaussPoints(2) = -0.339981043584856
          GaussWeights(2) = 0.652145154862546
          GaussPoints(3) =  0.339981043584856
          GaussWeights(3) = 0.652145154862546
          GaussPoints(4) =  0.861136311594953
          GaussWeights(4) = 0.347854845137454
        CASE(5) 
          GaussPoints(1) = -0.906179845938664
          GaussWeights(1) = 0.236926885056189
          GaussPoints(2) = -0.538469310105683
          GaussWeights(2) = 0.478628670499366
          GaussPoints(3) =  0.0
          GaussWeights(3) = 0.568888888888889
          GaussPoints(4) =  0.538469310105683
          GaussWeights(4) = 0.478628670499366
          GaussPoints(5) =  0.906179845938664
          GaussWeights(5) = 0.236926885056189
        CASE(6) 
          GaussPoints(1) = -0.932469514203152
          GaussWeights(1) = 0.171324492379170
          GaussPoints(2) = -0.661209386466265
          GaussWeights(2) = 0.360761573048139
          GaussPoints(3) = -0.238619186083197
          GaussWeights(3) = 0.467913934572691
          GaussPoints(4) =  0.238619186083197
          GaussWeights(4) = 0.467913934572691
          GaussPoints(5) =  0.661209386466265
          GaussWeights(5) = 0.360761573048139
          GaussPoints(6) =  0.932469514203152
          GaussWeights(6) = 0.171324492379170
        CASE(7) 
          GaussPoints(1) = -0.949107912342759
          GaussWeights(1) = 0.129484966168870
          GaussPoints(2) = -0.741531185599394
          GaussWeights(2) = 0.279705391489277
          GaussPoints(3) = -0.405845151377397
          GaussWeights(3) = 0.381830050505119
          GaussPoints(4) =  0.0			        
          GaussWeights(4) = 0.417959183673469
          GaussPoints(5) =  0.405845151377397
          GaussWeights(5) = 0.381830050505119
          GaussPoints(6) =  0.741531185599394
          GaussWeights(6) = 0.279705391489277
          GaussPoints(7) =  0.949107912342759
          GaussWeights(7) = 0.129484966168870
        CASE(8) 
          GaussPoints(1) = -0.96028986
          GaussWeights(1) = 0.10122854
          GaussPoints(2) = -0.79666648
          GaussWeights(2) = 0.22238103
          GaussPoints(3) = -0.52553241
          GaussWeights(3) = 0.31370665
          GaussPoints(4) = -0.18343464
          GaussWeights(4) = 0.36268378
          GaussPoints(5) =  0.18343464
          GaussWeights(5) = 0.36268378
          GaussPoints(6) =  0.52553241
          GaussWeights(6) = 0.31370665
          GaussPoints(7) =  0.79666648
          GaussWeights(7) = 0.22238103
          GaussPoints(8) =  0.96028986
          GaussWeights(8) = 0.10122854
        CASE DEFAULT
          WRITE(*,*) "Error: Invalid nGP in GetGaussPoints1D"
      END SELECT
      END SUBROUTINE GetGaussPoints1D
!=======================================================================
      SUBROUTINE getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts)
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: nGP
      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP), 
     &                    gwts(nGP)
      INTEGER ::  nn, ii, jj, kk, ll
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gpoints1, 
     &                                               gweights1

      IF(nGP == 1) THEN
        nn = 1
      ELSE IF(nGP == 8) THEN
        nn = 2
      ELSE IF(nGP == 27) THEN
        nn = 3
      ELSE IF(nGP == 64) THEN
        nn = 4
      ELSE
        nn = 5
      ENDIF

      ALLOCATE(gpoints1(nn))
      ALLOCATE(gweights1(nn))

      call GetGaussPoints1D(nn, gpoints1, gweights1)

      ll=1
      DO kk=1, nn
        DO jj=1, nn
          DO ii=1, nn
            gpts1(ll) = gpoints1(ii)
            gpts2(ll) = gpoints1(jj)
            gpts3(ll) = gpoints1(kk)
            gwts(ll)  = gweights1(kk)*gweights1(jj)*gweights1(ii)
            ll = ll+1
          END DO
        END DO
      END DO

      DEALLOCATE(gpoints1)
      DEALLOCATE(gweights1)

      END SUBROUTINE getGaussPointsHexa
!=======================================================================
      SUBROUTINE LagrangeBasisFunsTet(degree, xi1, xi2, xi3,
     &  N, dN_dxi1, dN_dxi2, dN_dxi3)
      IMPLICIT NONE
      INTEGER :: degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi3

       SELECT CASE (degree)
        CASE (0)
          N(1) = 1.0
          dN_dxi1(1) = 0.0
          dN_dxi2(1) = 0.0
          dN_dxi3(1) = 0.0
        CASE (1)
          N(1) = xi1
          N(2) = xi2
          N(3) = 1.0 - xi1 - xi2 - xi3
          N(4) = xi3

          dN_dxi1(1) =  1.0; dN_dxi1(2) =  0.0
          dN_dxi1(3) = -1.0; dN_dxi1(4) =  0.0

          dN_dxi2(1) =  0.0; dN_dxi2(2) =  1.0
          dN_dxi2(3) = -1.0; dN_dxi2(4) =  0.0

          dN_dxi3(1) =  0.0; dN_dxi3(2) =  0.0
          dN_dxi3(3) = -1.0; dN_dxi3(4) =  1.0
      CASE DEFAULT
          WRITE(*,*) "Error: Invalid degree in LagrangeBasisFunsTet"
      END SELECT
      END SUBROUTINE LagrangeBasisFunsTet
!=======================================================================
      SUBROUTINE LagrangeBasisFunsPrism(degree, xi1, xi2, xi4, 
     &  N, dN_dxi1, dN_dxi2, dN_dxi4)
      IMPLICIT NONE
      INTEGER :: degree
      DOUBLE PRECISION, INTENT(IN) :: xi1, xi2, xi4
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi4
      DOUBLE PRECISION :: xi3

      xi3 = 1.0 - xi1 - xi2

      SELECT CASE (degree)
        case (1)
          N(1) = xi3*(0.5*(1.0-xi4))
          N(2) = xi1*(0.5*(1.0-xi4))
          N(3) = xi2*(0.5*(1.0-xi4))
          N(4) = xi3*(0.5*(1.0+xi4))
          N(5) = xi1*(0.5*(1.0+xi4))
          N(6) = xi2*(0.5*(1.0+xi4))

          dN_dxi1(1)  = (-1.0)*(0.5*(1.0-xi4))
          dN_dxi1(2)  = ( 1.0)*(0.5*(1.0-xi4))
          dN_dxi1(3)  = ( 0.0)*(0.5*(1.0-xi4))
          dN_dxi1(4)  = (-1.0)*(0.5*(1.0+xi4))
          dN_dxi1(5)  = ( 1.0)*(0.5*(1.0+xi4))
          dN_dxi1(6)  = ( 0.0)*(0.5*(1.0+xi4))

          dN_dxi2(1)  = (-1.0)*(0.5*(1.0-xi4))
          dN_dxi2(2)  = ( 0.0)*(0.5*(1.0-xi4))
          dN_dxi2(3)  = ( 1.0)*(0.5*(1.0-xi4))
          dN_dxi2(4)  = (-1.0)*(0.5*(1.0+xi4))
          dN_dxi2(5)  = ( 0.0)*(0.5*(1.0+xi4))
          dN_dxi2(6)  = ( 1.0)*(0.5*(1.0+xi4))

          dN_dxi4(1)  = xi3*(-0.5)
          dN_dxi4(2)  = xi1*(-0.5)
          dN_dxi4(3)  = xi2*(-0.5)
          dN_dxi4(4)  = xi3*( 0.5)
          dN_dxi4(5)  = xi1*( 0.5)
          dN_dxi4(6)  = xi2*( 0.5)
        CASE DEFAULT
          WRITE(*,*) "Error: Invalid degree in LagrangeBasisFunsPrism"
      END SELECT
      END SUBROUTINE LagrangeBasisFunsPrism
!=======================================================================
      SUBROUTINE LagrangeBasisFunsHex(nNode, degree, xi1, xi2, xi3, 
     &  N, dN_dxi1, dN_dxi2, dN_dxi3)
      IMPLICIT NONE
      INTEGER :: nNode, degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION :: N(nNode), dN_dxi1(nNode), dN_dxi2(nNode),
     &  dN_dxi3(nNode)
      DOUBLE PRECISION :: v11, v12, v21, v22, v31, v32

      SELECT CASE (degree)
        case (0)
          N(1) = 1.0
          dN_dxi1(1) = 0.0; dN_dxi2(1) = 0.0; dN_dxi3(1) = 0.0

        case (1) ! linear hexahedron element
          v11 = 1.0 - xi1; v12 = 1.0 + xi1
          v21 = 1.0 - xi2; v22 = 1.0 + xi2
          v31 = 1.0 - xi3; v32 = 1.0 + xi3
          
          N(1) = 0.125*v11*v21*v31
          N(2) = 0.125*v12*v21*v31
          N(4) = 0.125*v11*v22*v31
          N(3) = 0.125*v12*v22*v31
          N(5) = 0.125*v11*v21*v32
          N(6) = 0.125*v12*v21*v32
          N(8) = 0.125*v11*v22*v32
          N(7) = 0.125*v12*v22*v32

          dN_dxi1(1) = -0.125*v21*v31; dN_dxi1(2) =  0.125*v21*v31
          dN_dxi1(4) = -0.125*v22*v31; dN_dxi1(3) =  0.125*v22*v31
          dN_dxi1(5) = -0.125*v21*v32; dN_dxi1(6) =  0.125*v21*v32
          dN_dxi1(8) = -0.125*v22*v32; dN_dxi1(7) =  0.125*v22*v32

          dN_dxi2(1) = -0.125*v11*v31; dN_dxi2(2) = -0.125*v12*v31
          dN_dxi2(4) =  0.125*v11*v31; dN_dxi2(3) =  0.125*v12*v31
          dN_dxi2(5) = -0.125*v11*v32; dN_dxi2(6) = -0.125*v12*v32
          dN_dxi2(8) =  0.125*v11*v32; dN_dxi2(7) =  0.125*v12*v32

          dN_dxi3(1) = -0.125*v11*v21; dN_dxi3(2) = -0.125*v12*v21
          dN_dxi3(4) = -0.125*v11*v22; dN_dxi3(3) = -0.125*v12*v22
          dN_dxi3(5) =  0.125*v11*v21; dN_dxi3(6) =  0.125*v12*v21
          dN_dxi3(8) =  0.125*v11*v22; dN_dxi3(7) =  0.125*v12*v22

        case (2) ! serendipity hexahedron element
          N(1) =  0.125*(-1+xi1)*(-1+xi2)*(-1+xi3)*(2+xi1+xi2+xi3)
          N(2) = -0.125*(1+xi1)*(-1+xi2)*(-1+xi3)*(2-xi1+xi2+xi3)
          N(3) = -0.125*(1+xi1)*(1+xi2)*(-1+xi3)*(-2+xi1+xi2-xi3)
          N(4) = -0.125*(-1+xi1)*(1+xi2)*(-1+xi3)*(2+xi1-xi2+xi3)
          N(5) = -0.125*(-1+xi1)*(-1+xi2)*(2+xi1+xi2-xi3)*(1+xi3)
          N(6) = -0.125*(1+xi1)*(-1+xi2)*(1+xi3)*(-2+xi1-xi2+xi3)
          N(7) =  0.125*(1+xi1)*(1+xi2)*(1+xi3)*(-2+xi1+xi2+xi3)
          N(8) =  0.125*(-1+xi1)*(1+xi2)*(2+xi1-xi2-xi3)*(1+xi3)

          N(9)  = -0.25*(-1+xi1**2)*(-1+xi2)*(-1+xi3)
          N(10) =  0.25*(1+xi1)*(-1+xi2**2)*(-1+xi3)
          N(11) =  0.25*(-1+xi1**2)*(1+xi2)*(-1+xi3)
          N(12) = -0.25*(-1+xi1)*(-1+xi2**2)*(-1+xi3)
          N(13) =  0.25*(-1+xi1**2)*(-1+xi2)*(1+xi3)
          N(14) = -0.25*(1+xi1)*(-1+xi2**2)*(1+xi3)
          N(15) = -0.25*(-1+xi1**2)*(1+xi2)*(1+xi3)
          N(16) =  0.25*(-1+xi1)*(-1+xi2**2)*(1+xi3)
          N(17) = -0.25*(-1+xi1)*(-1+xi2)*(-1+xi3**2)
          N(18) =  0.25*(1+xi1)*(-1+xi2)*(-1+xi3**2)
          N(19) = -0.25*(1+xi1)*(1+xi2)*(-1+xi3**2)
          N(20) =  0.25*(-1+xi1)*(1+xi2)*(-1+xi3**2)

          ! dN_dxi1
          dN_dxi1(1) =  0.125*(-1+xi2)*(-1+xi3)*(1+2*xi1+xi2+xi3)
          dN_dxi1(2) = -0.125*(-1+xi2)*(-1+xi3)*(1-2*xi1+xi2+xi3)
          dN_dxi1(3) =  0.125*(1+xi2)*(-1+xi3)*(1-2*xi1-xi2+xi3)
          dN_dxi1(4) =  0.125*(1+xi2)*(-1-2*xi1+xi2-xi3)*(-1+xi3)
          dN_dxi1(5) = -0.125*(-1+xi2)*(1+2*xi1+xi2-xi3)*(1+xi3)
          dN_dxi1(6) =  0.125*(-1+xi2)*(1-2*xi1+xi2-xi3)*(1+xi3)
          dN_dxi1(7) =  0.125*(1+xi2)*(1+xi3)*(-1+2*xi1+xi2+xi3)
          dN_dxi1(8) = -0.125*(1+xi2)*(1+xi3)*(-1-2*xi1+xi2+xi3)
          dN_dxi1(9)  = -0.5*xi1*(-1+xi2)*(-1+xi3)
          dN_dxi1(10) =  0.25*(-1+xi2**2)*(-1+xi3)
          dN_dxi1(11) =  0.5*xi1*(1+xi2)*(-1+xi3)
          dN_dxi1(12) = -0.25*(-1+xi2**2)*(-1+xi3)
          dN_dxi1(13) =  0.5*xi1*(-1+xi2)*(1+xi3)
          dN_dxi1(14) = -0.25*(-1+xi2**2)*(1+xi3)
          dN_dxi1(15) = -0.5*xi1*(1+xi2)*(1+xi3)
          dN_dxi1(16) =  0.25*(-1+xi2**2)*(1+xi3)
          dN_dxi1(17) = -0.25*(-1+xi2)*(-1+xi3**2)
          dN_dxi1(18) =  0.25*(-1+xi2)*(-1+xi3**2)
          dN_dxi1(19) = -0.25*(1+xi2)*(-1+xi3**2)
          dN_dxi1(20) =  0.25*(1+xi2)*(-1+xi3**2)

          ! dN_dxi2
          dN_dxi2(1) =  0.125*(-1+xi1)*(-1+xi3)*(1+xi1+2*xi2+xi3)
          dN_dxi2(2) =  0.125*(1+xi1)*(-1+xi1-2*xi2-xi3)*(-1+xi3)
          dN_dxi2(3) =  0.125*(1+xi1)*(-1+xi3)*(1-xi1-2*xi2+xi3)
          dN_dxi2(4) = -0.125*(-1+xi1)*(-1+xi3)*(1+xi1-2*xi2+xi3)
          dN_dxi2(5) = -0.125*(-1+xi1)*(1+xi1+2*xi2-xi3)*(1+xi3)
          dN_dxi2(6) = -0.125*(1+xi1)*(1+xi3)*(-1+xi1-2*xi2+xi3)
          dN_dxi2(7) =  0.125*(1+xi1)*(1+xi3)*(-1+xi1+2*xi2+xi3)
          dN_dxi2(8) =  0.125*(-1+xi1)*(1+xi1-2*xi2-xi3)*(1+xi3)
          dN_dxi2(9)  = -0.25*(-1+xi1**2)*(-1+xi3)
          dN_dxi2(10) =  0.5*(1+xi1)*xi2*(-1+xi3)
          dN_dxi2(11) =  0.25*(-1+xi1**2)*(-1+xi3)
          dN_dxi2(12) = -0.5*(-1+xi1)*xi2*(-1+xi3)
          dN_dxi2(13) =  0.25*(-1+xi1**2)*(1+xi3)
          dN_dxi2(14) = -0.5*(1+xi1)*xi2*(1+xi3)
          dN_dxi2(15) = -0.25*(-1+xi1**2)*(1+xi3)
          dN_dxi2(16) =  0.5*(-1+xi1)*xi2*(1+xi3)
          dN_dxi2(17) = -0.25*(-1+xi1)*(-1+xi3**2)
          dN_dxi2(18) =  0.25*(1+xi1)*(-1+xi3**2)
          dN_dxi2(19) = -0.25*(1+xi1)*(-1+xi3**2)
          dN_dxi2(20) =  0.25*(-1+xi1)*(-1+xi3**2)

          ! dN_dxi3
          dN_dxi3(1) =  0.125*(-1+xi1)*(-1+xi2)*(1+xi1+xi2+2*xi3)
          dN_dxi3(2) =  0.125*(1+xi1)*(-1+xi2)*(-1+xi1-xi2-2*xi3)
          dN_dxi3(3) = -0.125*(1+xi1)*(1+xi2)*(-1+xi1+xi2-2*xi3)
          dN_dxi3(4) = -0.125*(-1+xi1)*(1+xi2)*(1+xi1-xi2+2*xi3)
          dN_dxi3(5) = -0.125*(-1+xi1)*(-1+xi2)*(1+xi1+xi2-2*xi3)
          dN_dxi3(6) =  0.125*(1+xi1)*(-1+xi2)*(1-xi1+xi2-2*xi3)
          dN_dxi3(7) =  0.125*(1+xi1)*(1+xi2)*(-1+xi1+xi2+2*xi3)
          dN_dxi3(8) =  0.125*(-1+xi1)*(1+xi2)*(1+xi1-xi2-2*xi3)
          dN_dxi3(9)  = -0.25*(-1+xi1**2)*(-1+xi2)
          dN_dxi3(10) =  0.25*(1+xi1)*(-1+xi2**2)
          dN_dxi3(11) =  0.25*(-1+xi1**2)*(1+xi2)
          dN_dxi3(12) = -0.25*(-1+xi1)*(-1+xi2**2)
          dN_dxi3(13) =  0.25*(-1+xi1**2)*(-1+xi2)
          dN_dxi3(14) = -0.25*(1+xi1)*(-1+xi2**2)
          dN_dxi3(15) = -0.25*(-1+xi1**2)*(1+xi2)
          dN_dxi3(16) =  0.25*(-1+xi1)*(-1+xi2**2)
          dN_dxi3(17) = -0.5*(-1+xi1)*(-1+xi2)*xi3
          dN_dxi3(18) =  0.5*(1+xi1)*(-1+xi2)*xi3
          dN_dxi3(19) = -0.5*(1+xi1)*(1+xi2)*xi3
          dN_dxi3(20) =  0.5*(-1+xi1)*(1+xi2)*xi3

        CASE DEFAULT
          WRITE(*,*) "Error: Invalid degree in LagrangeBasisFunsHex"
      END SELECT
      END SUBROUTINE LagrangeBasisFunsHex
!=======================================================================
      SUBROUTINE computeBasisFunctions3D(nNode, ETYPE, degree, 
     &  param, xNode, yNode, zNode, N, dN_dx, dN_dy, dN_dz, Jac)
      IMPLICIT NONE
      INTEGER :: nNode, ETYPE, degree
      DOUBLE PRECISION :: param(3), Jac
      DOUBLE PRECISION :: xNode(nNode), yNode(nNode), zNode(nNode)
      DOUBLE PRECISION :: N(nNode), dN_dx(nNode), dN_dy(nNode),
     &                    dN_dz(nNode)
      INTEGER :: ii
      DOUBLE PRECISION :: dN_du1(nNode), dN_du2(nNode), dN_du3(nNode)
      DOUBLE PRECISION :: xx, yy, zz, detinv
      DOUBLE PRECISION :: B(3,3), Binv(3,3)

      IF(ETYPE == 4) THEN ! tetrahedral
        call LagrangeBasisFunsTet(degree, param(1), param(2), param(3),
     &    N, dN_du1, dN_du2, dN_du3)
      ELSE IF(ETYPE == 6) THEN ! prism
        call LagrangeBasisFunsPrism(degree, param(1), param(2), param(3),
     &    N, dN_du1, dN_du2, dN_du3)
      ELSE ! hex
        call LagrangeBasisFunsHex(nNode, degree, 
     &    param(1), param(2), param(3), N, dN_du1, dN_du2, dN_du3)
      END IF

      ! Calculation of Jacobian Matrix B
      B = 0.0
      DO ii=1, nNode
        xx = xNode(ii)
        yy = yNode(ii)
        zz = zNode(ii)

        B(1,1) = B(1,1) + (xx * dN_du1(ii))
        B(2,1) = B(2,1) + (xx * dN_du2(ii))
        B(3,1) = B(3,1) + (xx * dN_du3(ii))

        B(1,2) = B(1,2) + (yy * dN_du1(ii))
        B(2,2) = B(2,2) + (yy * dN_du2(ii))
        B(3,2) = B(3,2) + (yy * dN_du3(ii))

        B(1,3) = B(1,3) + (zz * dN_du1(ii))
        B(2,3) = B(2,3) + (zz * dN_du2(ii))
        B(3,3) = B(3,3) + (zz * dN_du3(ii))
      END DO

      Jac = -B(1,3)*B(2,2)*B(3,1) + B(1,2)*B(2,3)*B(3,1)
     &      +B(1,3)*B(2,1)*B(3,2) - B(1,1)*B(2,3)*B(3,2)
     &      -B(1,2)*B(2,1)*B(3,3) + B(1,1)*B(2,2)*B(3,3)

      detinv = 1.0/Jac

      Binv(1,1) = +detinv * (B(2,2)*B(3,3) - B(2,3)*B(3,2))
      Binv(2,1) = -detinv * (B(2,1)*B(3,3) - B(2,3)*B(3,1))
      Binv(3,1) = +detinv * (B(2,1)*B(3,2) - B(2,2)*B(3,1))
      Binv(1,2) = -detinv * (B(1,2)*B(3,3) - B(1,3)*B(3,2))
      Binv(2,2) = +detinv * (B(1,1)*B(3,3) - B(1,3)*B(3,1))
      Binv(3,2) = -detinv * (B(1,1)*B(3,2) - B(1,2)*B(3,1))
      Binv(1,3) = +detinv * (B(1,2)*B(2,3) - B(1,3)*B(2,2))
      Binv(2,3) = -detinv * (B(1,1)*B(2,3) - B(1,3)*B(2,1))
      Binv(3,3) = +detinv * (B(1,1)*B(2,2) - B(1,2)*B(2,1))

      DO ii=1, nNode
        dN_dx(ii) = dN_du1(ii) * Binv(1,1) + dN_du2(ii) * Binv(1,2)
     &            + dN_du3(ii) * Binv(1,3)
        dN_dy(ii) = dN_du1(ii) * Binv(2,1) + dN_du2(ii) * Binv(2,2)
     &            + dN_du3(ii) * Binv(2,3)
        dN_dz(ii) = dN_du1(ii) * Binv(3,1) + dN_du2(ii) * Binv(3,2)
     &            + dN_du3(ii) * Binv(3,3)
      END DO

      END SUBROUTINE computeBasisFunctions3D
!=======================================================================
      SUBROUTINE DetMatrix(M, DetM)
      IMPLICIT NONE
      DOUBLE PRECISION :: M(3,3), DetM

      DetM = -M(1,3)*M(2,2)*M(3,1) + M(1,2)*M(2,3)*M(3,1)
     &       +M(1,3)*M(2,1)*M(3,2) - M(1,1)*M(2,3)*M(3,2)
     &       -M(1,2)*M(2,1)*M(3,3) + M(1,1)*M(2,2)*M(3,3)

      END SUBROUTINE DetMatrix
!=======================================================================
      SUBROUTINE InvMatrix(MInv, M, DetM)
      IMPLICIT NONE
      DOUBLE PRECISION :: MInv(3,3), M(3,3), DetM

      MInv(1,1) = (M(2,2)*M(3,3)-M(2,3)*M(3,2))/DetM
      MInv(1,2) = (M(1,3)*M(3,2)-M(1,2)*M(3,3))/DetM
      MInv(1,3) = (M(1,2)*M(2,3)-M(1,3)*M(2,2))/DetM
      MInv(2,1) = (M(2,3)*M(3,1)-M(2,1)*M(3,3))/DetM
      MInv(2,2) = (M(1,1)*M(3,3)-M(1,3)*M(3,1))/DetM
      MInv(2,3) = (M(1,3)*M(2,1)-M(1,1)*M(2,3))/DetM
      MInv(3,1) = (M(2,1)*M(3,2)-M(2,2)*M(3,1))/DetM
      MInv(3,2) = (M(1,2)*M(3,1)-M(1,1)*M(3,2))/DetM
      MInv(3,3) = (M(1,1)*M(2,2)-M(1,2)*M(2,1))/DetM

      END SUBROUTINE InvMatrix
!=======================================================================