! ---------------------------------------------------------------------
! ABAQUS user Subroutine UEL for incompressibel electro growth material
! in 3D case. Regarding to the formulation of this UEL, please refer to
! these papers: 
! A robust and computationally efficient finite element framework for
! coupled electromechanics. 10.1016/j.cma.2020.113443
! On the advantages of mixed formulation and higher-order elements for 
! computational morphoelasticity. 10.1016/j.jmps.2020.104289
! and a recent paper TODO: need to update link
! UEL by Zhanfeng Li 20241016 in Guangzhou
! ---------------------------------------------------------------------
! Anyone use this UEL or have inquires, please email Zhanfeng Li 
! through ctzf.li@mail.scut.edu.cn to get him noticed. Thank you! 
! A FORTRAN 77 statement may be up to 40 lines long! 
! Columns 7-72 contain the FORTRAN instructions.
!
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
!
! ---------------------------------------------------------------------
      module global

      ! This module is used to transfer SDV's from the UEL
      !  to the UVARM so that SDV's can be visualized on a
      !  dummy mesh
      !
      !  GloUVar(X,Y,Z)
      !   X - element pointer
      !   Y - integration point pointer
      !   Z - SDV pointer
      !
      !  NELEMENT
      !   Total number of elements in the real mesh, the dummy
      !   mesh needs to have the same number of elements, and
      !   the dummy mesh needs to have the same number of integ
      !   points.  You must set that parameter value here.
      !
      !  ELEMOFFSET
      !   Offset between element numbers on the real mesh and
      !    dummy mesh.  That is set in the input file, and
      !    that value must be set here the same.
      !
      !  nUvrm
      !   Total number of Uvarm variables

      integer NELEMENT, ELEMOFFSET, nUvrm, err

      ! Set the number of UEL elements used here
      parameter(NELEMENT=77280)

      ! Set the offset here for UVARM plotting, must match input file!
      parameter(ELEMOFFSET=100000)

      parameter(nUvrm=20)

      real*8, allocatable :: GloUVar(:,:,:)

      end module global
! ------------------------- UEL subroutine ----------------------------
! User coding to define RHS, AMATRX, SVARS, ENERGY, and PNEWDT
      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
!
      USE global
      INCLUDE 'ABA_PARAM.INC'
!
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

      DOUBLE PRECISION :: xNode(nNode), yNode(nNode), zNode(nNode), 
     & param(3), 
     & Jac, N(nNode), dN_dX(nNode), dN_dY(nNode), dN_dZ(nNode), 
     & GaussXc, GaussYc, GaussZc, 
     & Em, Nu, Mu, kappa, Lambda, Mu0, Im, MatFlag, gFactor, Pi, 
     & MRef(3), MCur(3), Ha(3), BExt(3), SizeFilm, 
     & MyStrain(9), CauchySS(9), vonMisesSS, Elasticity(9,9), 
     & F(3,3), FInv(3,3), b_A(3,3), Ib_A, DetF, 
     & Gc(3,9), dN_dxx(nNode), dN_dyy(nNode), dN_dzz(nNode), 
     & dvol0, dvol, ElVol0, b1, b2, b3, b4, 
     & gpts1_p(8), gpts2_p(8), gpts3_p(8), gwts_p(8), Np(8), 
     & xNode_p(8), yNode_p(8), zNode_p(8), pElem(8), ResP(8), 
     & dNp_dX(8), dNp_dY(8), dNp_dZ(8), Jac_p, Pressure, Phi, 
     & G(3,3), DetG, Ae(3,3), DetA, Cur_time, TotalT, 
     & phiElem(nNode), Nphi(nNode), Jac_phi,
     & xNode_phi(nNode), yNode_phi(nNode), zNode_phi(nNode), 
     & dNphi_dX(nNode), dNphi_dY(nNode), dNphi_dZ(nNode), 
     & dNphi_dxx(nNode), dNphi_dyy(nNode), dNphi_dzz(nNode), dphi(3),
     & pCouple1(9,3), pCouple2(9,3), Magneticity(3,3), 
     & Gu_p(3,3), Gu_p_Tmp(3,3), G_phi_d(3), 
     & Potential_Mech, Potential_Mag

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gpts1, gpts2, 
     &  gpts3, gwts, gpts1_phi, gpts2_phi, gpts3_phi, gwts_phi, 
     &  dispElem, ResU, Res_phi
      DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Kuu, Kup, 
     & Kpu, Kpp, K_u_phi, K_u_phi_Tmp, K_phi_u, K_phi_phi

      INTEGER :: gp,  ii, jj, kk, ll, PP, QQ, 
     &  degree=2, nGP=27, EleType=8, nDim=3, nDofDisp, nDofp, nDofphi,
     &  nGP_p=8, nNode_p=8, degree_p=1, 
     &  nGP_phi=27, nNode_phi=20, degree_phi=2

      ! UVAR(NELEMENT,1,NINPT) to UVAR(NELEMENT,9,NINPT) are 9 stress components 
      ! UVAR(NELEMENT,10,NINPT) is von-Mises stress
      ! UVAR(NELEMENT,11,NINPT) is internal pressure

      ! UVAR(NELEMENT,12,NINPT) is 
      ! UVAR(NELEMENT,13,NINPT) is 
      ! UVAR(NELEMENT,14,NINPT) is 

      ! UVAR(NELEMENT,15,NINPT) is 
      ! UVAR(NELEMENT,16,NINPT) is 
      ! UVAR(NELEMENT,17,NINPT) is 
      nDofDisp = nNode*nDim
      nDofphi = nNode
      nDofp = 8

      ALLOCATE(gpts1(nGP))
      ALLOCATE(gpts2(nGP))
      ALLOCATE(gpts3(nGP))
      ALLOCATE(gwts(nGP))

      ALLOCATE(gpts1_phi(nGP_phi))
      ALLOCATE(gpts2_phi(nGP_phi))
      ALLOCATE(gpts3_phi(nGP_phi))
      ALLOCATE(gwts_phi(nGP_phi))

      ALLOCATE(dispElem(nDofDisp))
      ALLOCATE(ResU(nDofDisp))
      ALLOCATE(Res_phi(nDofphi))

      ALLOCATE(Kuu(nDofDisp,nDofDisp))
      ALLOCATE(Kup(nDofDisp,nDofp))
      ALLOCATE(Kpu(nDofp,nDofDisp))
      ALLOCATE(Kpp(nDofp,nDofp))

      ALLOCATE(K_u_phi(nDofDisp,nDofphi))
      ALLOCATE(K_u_phi_Tmp(nDofDisp,nDofphi))
      ALLOCATE(K_phi_u(nDofphi,nDofDisp))
      ALLOCATE(K_phi_phi(nDofphi,nDofphi))

c Allocate memory for the global variables used for visualization
c results on the dummy mesh
      ALLOCATE(GloUVar(NELEMENT,nGP,nUvrm),stat=err)
! Parameters input ----------------------------------------------------
      ! WRITE(*,*) "Entering UEL, parameter input..."
      ! WRITE(*,*) "Rows of RHS are", MLVARX
      ! WRITE(*,*) "Rows/columns of AMATRX are", NDOFEL
      ! WRITE(*,*) "dofs of displacement", nDofDisp
      ! WRITE(*,*) "dofs of pressure", nDofp
      ! WRITE(*,*) "Number of nodes is", NNODE
      ! WRITE(*,*) "Number of Gauss points is", nGP
      ! WRITE(*,*) "Degree in UEL is", degree

      Cur_time = TIME(2)

      ! Material parameters 
      Em = PROPS(1)
      Nu = PROPS(2)
      Pi = 3.14159265359

      ! growth effect
      gFactor = PROPS(9)

      MatFlag = PROPS(10)

      ! Calculate Neo-Hookean parameters 
      Mu = Em/(2.0+2.0*Nu) ! use Pa
      kappa = Em/3.0/(1.0-2.0*Nu)
!       Lambda = Em*Nu/(1.0+Nu)/(1-2.0*Nu)

      ! Mu0 = 1.2566371/(10**6)
      ! customize the unit N (kA)^{-2}
      Mu0 = 1.2566371

      ! loading control
      TotalT = 1.0 

      IF (Cur_time .GT. TotalT) THEN
        ! assign magnetization vector Mr
        MRef(1) = PROPS(3)
        MRef(2) = PROPS(4)
        MRef(3) = PROPS(5)
        ! assign external magnetic induction vector Ba
        BExt(1) = PROPS(6)*(Cur_time - 1.0)
        BExt(2) = PROPS(7)*(Cur_time - 1.0)
        BExt(3) = PROPS(8)*(Cur_time - 1.0)

      ELSE
        MRef(1) = PROPS(3)*Cur_time
        MRef(2) = PROPS(4)*Cur_time
        MRef(3) = PROPS(5)*Cur_time
        
        BExt(1) = 0.0
        BExt(2) = 0.0
        BExt(3) = 0.0
      END IF

      ! Calculate external magnetic field Ha
      Ha(1) = BExt(1)/Mu0
      Ha(2) = BExt(2)/Mu0
      Ha(3) = BExt(3)/Mu0
      
      ! WRITE(*,*) "Mu, kappa, Mu0", Mu, kappa, Mu0
      ! WRITE(*,*) "Ha(3)", Ha(1),Ha(2),Ha(3)

! pass the original coordinates to xNode, yNode, zNode ----------------
      DO jj=1,nNode
        xNode(jj) = COORDS(1,jj)
        yNode(jj) = COORDS(2,jj)
        zNode(jj) = COORDS(3,jj)
        ! WRITE(*,*) "x,y,z of node",jj, xNode(jj), yNode(jj), zNode(jj)
      END DO

      DO jj=1,nNode
        xNode_phi(jj) = COORDS(1,jj)
        yNode_phi(jj) = COORDS(2,jj)
        zNode_phi(jj) = COORDS(3,jj)
      END DO

      ! obatin coordinates of nodes 1 to 8
      DO jj=1,nNode_p
        xNode_p(jj) = COORDS(1,jj)
        yNode_p(jj) = COORDS(2,jj)
        zNode_p(jj) = COORDS(3,jj)
      END DO

! obtain dispElem, phiElem and pElem from U(*) ------------------------
        ! for nodes 1 to 8
        DO ii=1, nNode_p
          dispElem(3*ii-2) = U(5*ii-4) ! x
          dispElem(3*ii-1) = U(5*ii-3) ! y
          dispElem(3*ii)   = U(5*ii-2) ! z
          phiElem(ii)      = U(5*ii-1) ! phi 
          pElem(ii)        = U(5*ii) ! p
        END DO

        ! for nodes 9 to 20
        DO ii=1+nNode_p,nNode
          dispElem(3*ii-2) = U(nDofp+4*ii-3)
          dispElem(3*ii-1) = U(nDofp+4*ii-2)
          dispElem(3*ii)   = U(nDofp+4*ii-1)
          phiElem(ii)      = U(nDofp+4*ii)
        END DO

!         DO ii=1,nNode_p
!           WRITE(*,*) "pElem(ii)", ii, pElem(ii)
!         END DO

!         DO ii=1,nNode
!           WRITE(*,*) "phiElem(ii)", ii, phiElem(ii)
!         END DO

! Gauss point coordinates and Weights ---------------------------------
      ! WRITE(*,*) "get Gauss points for displacement..."
      call getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts)

      ! WRITE(*,*) "get Gauss points for phi..." 
      call getGaussPointsHexa(nGP_phi, 
     &       gpts1_phi, gpts2_phi, gpts3_phi, gwts_phi)

      ! WRITE(*,*) "get Gauss points for pressure..."
      call getGaussPointsHexa(nGP_p, gpts1_p, gpts2_p, gpts3_p, gwts_p)

      Kuu = 0.0
      Kup = 0.0
      Kpu = 0.0
      K_u_phi = 0.0
      K_u_phi_Tmp = 0.0
      K_phi_u = 0.0
      Kpp = 0.0
      K_phi_phi = 0.0

      ResU = 0.0
      Resp = 0.0
      Res_phi = 0.0

! loop over Gauss points ----------------------------------------------
      DO gp=1, nGP
        ! WRITE(*,*) "Gauss point number = ", gp
        param(1) = gpts1(gp)
        param(2) = gpts2(gp)
        param(3) = gpts3(gp)

  ! get basis functions for displacement ------------------------------
        ! EleType == 8 for hexa elements
        ! WRITE(*,*) "Enter computeBasisFunctions3D..."
        call computeBasisFunctions3D(nNode, EleType, degree, param,
     &    xNode, yNode, zNode, N, dN_dX, dN_dY, dN_dZ, Jac)
        IF(Jac < 0.0) THEN
          WRITE(*,*) "Jac < 0.0 for the element, check orientation" 
        END IF

  ! get basis functions for phi ---------------------------------------
!         call computeBasisFunctions3D(nNode_phi, EleType, degree_phi, 
!      &    param, xNode_phi, yNode_phi, zNode_phi, 
!      &    Nphi, dNphi_dX, dNphi_dY, dNphi_dZ, Jac_phi)
        call computeBasisFunctions3D(nNode, EleType, degree, 
     &    param, xNode, yNode, zNode, 
     &    Nphi, dNphi_dX, dNphi_dY, dNphi_dZ, Jac_phi)

  ! get basis functions for pressure ----------------------------------
        call computeBasisFunctions3D(nNode_p, EleType, degree_p, param,
     &    xNode_p, yNode_p, zNode_p, Np, dNp_dX, dNp_dY, dNp_dZ, Jac_p)

  ! calculate the coordinate of this displacement Gauss point ---------
        GaussXc = 0
        GaussYc = 0
        GaussZc = 0

        DO ii=1, nNode
          GaussXc = GaussXc + xNode(ii) * N(ii)
          GaussYc = GaussYc + yNode(ii) * N(ii)
          GaussZc = GaussZc + zNode(ii) * N(ii)
        END DO

        ! calculate pressure in this displacement Gauss point
        Pressure = 0.0
        DO ii=1, nNode_p
          Pressure = Pressure + pElem(ii) * Np(ii)
        END DO

        ! calculate phi in this displacement Gauss point
        Phi = 0.0
        DO ii=1, nNode
          Phi = Phi + phiElem(ii) * Nphi(ii)
        END DO

        ! WRITE(*,*) "Np(ii), 1-4", Np(1),Np(2),Np(3),Np(4)
        ! WRITE(*,*) "Np(ii), 5-8", Np(5),Np(6),Np(7),Np(8)

  ! prescribe magnetic field TODO: ------------------------------------
      ! distribution of magnetization vector Mr, for bulge bendng 
      ! SizeFilm = 60.0
      ! IF (GaussXc .GT. 0) THEN
      !   MRef(1) = -PROPS(3)*Cur_time
      ! END IF

  ! Growth functions input --------------------------------------------
        G = 0
        G(1,1) = 1
        G(2,2) = 1
        G(3,3) = 1

        ! growth only happens in 1st stage 
        IF (Cur_time .GT. TotalT) THEN
          Cur_time = TotalT ! freeze growth
        END IF

        call GrowthFun(G, GaussXc,GaussYc,GaussZc, Cur_time, TotalT, 
     &                 gFactor)

  ! Calculate the deformation gradient tensor F on Gauss point --------
        F = 0.0
        F(1,1) = 1.0
        F(2,2) = 1.0
        F(3,3) = 1.0

        DO ii=1, nNode
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

  ! Inverse of F ------------------------------------------------------
        call DetMatrix(F, DetF)
        ! WRITE(*,*) "#gp, DetF = ", gp, DetF
        call InvMatrix(FInv, F, DetF)

  ! Calculate the Magnetization vector MCur on Gauss point ------------
        MCur(1) = (F(1,1)*MRef(1)+F(1,2)*MRef(2)+F(1,3)*MRef(3))/DetF
        MCur(2) = (F(2,1)*MRef(1)+F(2,2)*MRef(2)+F(2,3)*MRef(3))/DetF
        MCur(3) = (F(3,1)*MRef(1)+F(3,2)*MRef(2)+F(3,3)*MRef(3))/DetF

        ! WRITE(*,*) "MCur(3)", MCur(1),MCur(2),MCur(3)

  ! elastic defoemation gradient tensor A = F G^(-1) ------------------
        call CalAe(Ae, F, G)

  ! Det of Ae ---------------------------------------------------------
        call DetMatrix(Ae, DetA)

  ! left Cauchy-Green tensor b_A=A.A^T --------------------------------
        b_A = 0.0
        DO ii=1, 3
          DO jj=1, 3
            DO KK=1, 3
              b_A(ii,jj) = b_A(ii,jj) + Ae(ii,KK)*Ae(jj,KK)
            END DO
          END DO
        END DO
  ! the first invariant of tensor b_A
        Ib_A = b_A(1,1) + b_A(2,2) + b_A(3,3)

  ! calculate d(phi)d(i) 3×1
        dphi = 0
        DO ii=1, nNode
          dNphi_dxx(ii) = dNphi_dX(ii)*FInv(1,1) 
     &      + dNphi_dY(ii)*FInv(2,1) + dNphi_dZ(ii)*FInv(3,1)
          dNphi_dyy(ii) = dNphi_dX(ii)*FInv(1,2) 
     &      + dNphi_dY(ii)*FInv(2,2) + dNphi_dZ(ii)*FInv(3,2)
          dNphi_dzz(ii) = dNphi_dX(ii)*FInv(1,3) 
     &      + dNphi_dY(ii)*FInv(2,3) + dNphi_dZ(ii)*FInv(3,3)
        END DO

        DO ii=1, nNode
          dphi(1) = dphi(1) + dNphi_dxx(ii)*phiElem(ii)
          dphi(2) = dphi(2) + dNphi_dyy(ii)*phiElem(ii)
          dphi(3) = dphi(3) + dNphi_dzz(ii)*phiElem(ii)
        END DO
        ! WRITE(*,*) "#gp, dphi = ", gp, dphi(1), dphi(2), dphi(3)

  ! Elasticity tensor ee of the constitutive model --------------------
  ! Elasticity row = Sigma {11,21,31,12,22,32,13,23,33}
  ! Elasticity column = Mu0 {11,21,31,12,22,32,13,23,33}
  ! Notice Elasticity here is different from that in UMAT
        ! WRITE(*,*) "Calculating ElasticityTensor..."
        Elasticity = 0.0 
        call ElasticityTensor(Elasticity, Pressure, Mu, kappa,
     &     DetA, Ib_A, b_A, Mu0, dphi, MCur)

!         call ElasticityTensorGent(Elasticity, Pressure, Mu, kappa,
!      &     DetA, Ib_A, b_A, Mu0, dphi, Im)

!         DO ii=1, 9
!           WRITE(*,*) "ElasticityTensor", Elasticity(ii,1), 
!      &      Elasticity(ii,2), Elasticity(ii,3), Elasticity(ii,4), 
!      &      Elasticity(ii,5), Elasticity(ii,6), Elasticity(ii,7), 
!      &      Elasticity(ii,8), Elasticity(ii,9)
!         END DO

  ! Coupling tensor p_ijk ---------------------------------------------
  ! CouplingTensor row = ij {11,12,13,21,22,23,31,32,33}
  ! CouplingTensor column = k {1,2,3}
        ! WRITE(*,*) "Calculating CouplingTensor..."
        pCouple1 = 0.0 
        call CouplingTensor1(pCouple1, Mu0, dphi, DetA, Ae, MCur)

        pCouple2 = 0.0 
        call CouplingTensor2(pCouple2, Mu0, dphi, DetA, Ae, MCur)

  ! Magnetic permittivity tensor d_ij ---------------------------------
        ! WRITE(*,*) "Calculating MagneticityTensor..."
        Magneticity = 0.0 
        call MagneticityTensor(Magneticity, Mu0)

  ! Calculate my effective Cauchy stress {11,12,13,21,22,23,31,32,33}
        ! WRITE(*,*) "Calculating CauchyStress..."
        CauchySS = 0.0
        call cal_CauchyStress(CauchySS, Pressure, Mu, kappa,
     &     DetA, Ib_A, b_A, Mu0, dphi, Ae, MCur, Ha)

!         call cal_CauchyStressGent(CauchySS, Pressure, Mu, kappa,
!      &     DetA, Ib_A, b_A, Mu0, dphi, Im)

        vonMisesSS = Sqrt(0.5*((CauchySS(1)-CauchySS(5))**2
     &    +(CauchySS(5)-CauchySS(9))**2+(CauchySS(9)-CauchySS(1))**2)
     &    +3*(CauchySS(2)*CauchySS(2)+CauchySS(6)*CauchySS(6)
     &    +CauchySS(7)*CauchySS(7)))

  ! calculate derivative to the current coordinate denoted by dN_dxx, dN_dyy and dN_dzz
        dN_dxx = 0.0
        dN_dyy = 0.0
        dN_dzz = 0.0

        DO ii=1, nNode
          dN_dxx(ii) = dN_dX(ii)*FInv(1,1) 
     &               + dN_dY(ii)*FInv(2,1) + dN_dZ(ii)*FInv(3,1)
          dN_dyy(ii) = dN_dX(ii)*FInv(1,2) 
     &               + dN_dY(ii)*FInv(2,2) + dN_dZ(ii)*FInv(3,2)
          dN_dzz(ii) = dN_dX(ii)*FInv(1,3) 
     &               + dN_dY(ii)*FInv(2,3) + dN_dZ(ii)*FInv(3,3)
        END DO

  ! difference of volume
        DetG  = DetF/DetA
        dvol0 = gwts(gp) * Jac 
        dvol  = gwts(gp) * Jac * DetF

  ! calculate strain energy function
        Potential_Mech = DetG*( Gm/2.0*(DetA**(-2.0/3.0)*Ib_A-3.0)
     &         +Pressure*(DetA-1.0) -Pressure*Pressure/(2.0*kappa) )
        Potential_Mag = -Mu0*DetF*( (Ha(1)-dphi(1))*MCur(1) 
     &  + (Ha(2)-dphi(2))*MCur(2) + (Ha(3)-dphi(3))*MCur(3) ) 
     &  - 0.5*Mu0*DetF*( dphi(1)*dphi(1) + dphi(2)*dphi(2) 
     &  + dphi(3)*dphi(3) )

  ! TODO: pass stress components to global variable
        DO kk=1,9
          GloUVar(JELEM,gp,kk) = CauchySS(kk)
        END DO
        GloUVar(JELEM,gp,10) = vonMisesSS
        GloUVar(JELEM,gp,11) = Pressure
        GloUVar(JELEM,gp,12) = Phi
        GloUVar(JELEM,gp,13) = Sqrt((dphi(1)**2) 
     &                              + (dphi(2)**2) + (dphi(3)**2))
        GloUVar(JELEM,gp,14) = G(1,1)
        GloUVar(JELEM,gp,15) = MRef(3)
        GloUVar(JELEM,gp,16) = BExt(3)
        GloUVar(JELEM,gp,17) = Potential_Mech
        GloUVar(JELEM,gp,18) = Potential_Mag

        Gc = 0.0
        DO ii=1, nNode  ! summation over dN_d(ii) ---------------------
          b1 = dN_dxx(ii)*dvol
          b2 = dN_dyy(ii)*dvol
          b3 = dN_dzz(ii)*dvol
          ! b4 = N(ii)*dvol

          ! Gu^T c, c is the elasticity tensor 3×9
          DO pp=1, 9
            Gc(1,pp) = b1*Elasticity(1,pp) 
     &        + b2*Elasticity(4,pp) + b3*Elasticity(7,pp)
            Gc(2,pp) = b1*Elasticity(2,pp) 
     &        + b2*Elasticity(5,pp) + b3*Elasticity(8,pp)
            Gc(3,pp) = b1*Elasticity(3,pp) 
     &        + b2*Elasticity(6,pp) + b3*Elasticity(9,pp)
          END DO

          ! Gu^T p, p is the coupling tensor 3×3
          DO pp=1, 3
            Gu_p(1,pp) = b1*pCouple1(1,pp) 
     &        + b2*pCouple1(4,pp) + b3*pCouple1(7,pp)
            Gu_p(2,pp) = b1*pCouple1(2,pp) 
     &        + b2*pCouple1(5,pp) + b3*pCouple1(8,pp)
            Gu_p(3,pp) = b1*pCouple1(3,pp) 
     &        + b2*pCouple1(6,pp) + b3*pCouple1(9,pp)
          END DO

          ! Gu^T p, p is the coupling tensor 3×3
          DO pp=1, 3
            Gu_p_Tmp(1,pp) = b1*pCouple2(1,pp) 
     &        + b2*pCouple2(4,pp) + b3*pCouple2(7,pp)
            Gu_p_Tmp(2,pp) = b1*pCouple2(2,pp) 
     &        + b2*pCouple2(5,pp) + b3*pCouple2(8,pp)
            Gu_p_Tmp(3,pp) = b1*pCouple2(3,pp) 
     &        + b2*pCouple2(6,pp) + b3*pCouple2(9,pp)
          END DO

          ! stiffness matrix 1 Kuu  60×60
          DO jj=1, nNode
            Kuu(3*ii-2, 3*jj-2) = Kuu(3*ii-2, 3*jj-2) 
     &                + Gc(1,1)*dN_dxx(jj) + Gc(1,4)*dN_dyy(jj) 
     &                + Gc(1,7)*dN_dzz(jj)
            Kuu(3*ii-2, 3*jj-1) = Kuu(3*ii-2, 3*jj-1)
     &                + Gc(1,2)*dN_dxx(jj) + Gc(1,5)*dN_dyy(jj) 
     &                + Gc(1,8)*dN_dzz(jj)
            Kuu(3*ii-2, 3*jj)   = Kuu(3*ii-2, 3*jj)  
     &                + Gc(1,3)*dN_dxx(jj) + Gc(1,6)*dN_dyy(jj) 
     &                + Gc(1,9)*dN_dzz(jj)

            Kuu(3*ii-1, 3*jj-2) = Kuu(3*ii-1, 3*jj-2)
     &                + Gc(2,1)*dN_dxx(jj) + Gc(2,4)*dN_dyy(jj) 
     &                + Gc(2,7)*dN_dzz(jj)
            Kuu(3*ii-1, 3*jj-1) = Kuu(3*ii-1, 3*jj-1)
     &                + Gc(2,2)*dN_dxx(jj) + Gc(2,5)*dN_dyy(jj) 
     &                + Gc(2,8)*dN_dzz(jj)
            Kuu(3*ii-1, 3*jj)   = Kuu(3*ii-1, 3*jj)  
     &                + Gc(2,3)*dN_dxx(jj) + Gc(2,6)*dN_dyy(jj) 
     &                + Gc(2,9)*dN_dzz(jj)

            Kuu(3*ii, 3*jj-2) = Kuu(3*ii, 3*jj-2)
     &                + Gc(3,1)*dN_dxx(jj) + Gc(3,4)*dN_dyy(jj) 
     &                + Gc(3,7)*dN_dzz(jj)
            Kuu(3*ii, 3*jj-1) = Kuu(3*ii, 3*jj-1)
     &                + Gc(3,2)*dN_dxx(jj) + Gc(3,5)*dN_dyy(jj) 
     &                + Gc(3,8)*dN_dzz(jj)
            Kuu(3*ii, 3*jj)   = Kuu(3*ii, 3*jj)  
     &                + Gc(3,3)*dN_dxx(jj) + Gc(3,6)*dN_dyy(jj) 
     &                + Gc(3,9)*dN_dzz(jj)
          END DO

          ! stiffness matrix 2 K_u_phi 60×20
          DO jj=1, nNode
            K_u_phi(3*ii-2,jj) = K_u_phi(3*ii-2,jj)      
     &          + Gu_p(1,1)*dNphi_dxx(jj) + Gu_p(1,2)*dNphi_dyy(jj) 
     &          + Gu_p(1,3)*dNphi_dzz(jj)
            K_u_phi(3*ii-1,jj) = K_u_phi(3*ii-1,jj)       
     &          + Gu_p(2,1)*dNphi_dxx(jj) + Gu_p(2,2)*dNphi_dyy(jj) 
     &          + Gu_p(2,3)*dNphi_dzz(jj)
            K_u_phi(3*ii,jj)   = K_u_phi(3*ii,jj)       
     &          + Gu_p(3,1)*dNphi_dxx(jj) + Gu_p(3,2)*dNphi_dyy(jj) 
     &          + Gu_p(3,3)*dNphi_dzz(jj)
          END DO

          ! stiffness matrix 2 K_u_phi_Tmp 60×20
          DO jj=1, nNode
            K_u_phi_Tmp(3*ii-2,jj) = K_u_phi_Tmp(3*ii-2,jj)      
     &     + Gu_p_Tmp(1,1)*dNphi_dxx(jj) + Gu_p_Tmp(1,2)*dNphi_dyy(jj) 
     &          + Gu_p_Tmp(1,3)*dNphi_dzz(jj)
            K_u_phi_Tmp(3*ii-1,jj) = K_u_phi_Tmp(3*ii-1,jj)       
     &     + Gu_p_Tmp(2,1)*dNphi_dxx(jj) + Gu_p_Tmp(2,2)*dNphi_dyy(jj) 
     &          + Gu_p_Tmp(2,3)*dNphi_dzz(jj)
            K_u_phi_Tmp(3*ii,jj)   = K_u_phi_Tmp(3*ii,jj)       
     &     + Gu_p_Tmp(3,1)*dNphi_dxx(jj) + Gu_p_Tmp(3,2)*dNphi_dyy(jj) 
     &          + Gu_p_Tmp(3,3)*dNphi_dzz(jj)
          END DO

          ! stiffness matrix 3 Kup 60×8
          DO jj=1, nNode_p
            Kup(3*ii-2,jj) = Kup(3*ii-2,jj) + b1*Np(jj)
            Kup(3*ii-1,jj) = Kup(3*ii-1,jj) + b2*Np(jj)
            Kup(3*ii,jj)   = Kup(3*ii,jj)   + b3*Np(jj)
          END DO

          ! G_phi^T d, d is the electric permittivity tensor 1×3
          DO pp=1, 3
            G_phi_d(pp) = dNphi_dxx(ii)*dvol*Magneticity(1,pp) 
     &                   + dNphi_dyy(ii)*dvol*Magneticity(2,pp) 
     &                   + dNphi_dzz(ii)*dvol*Magneticity(3,pp) 
          END DO

        ! stiffness matrix 4 K_phi_phi, summation over gauss points 20×20
          DO jj=1, nNode
            K_phi_phi(ii,jj) = K_phi_phi(ii,jj) 
     &                        + G_phi_d(1)*dNphi_dxx(jj)
     &                        + G_phi_d(2)*dNphi_dyy(jj)
     &                        + G_phi_d(3)*dNphi_dzz(jj)
          END DO

          ! residual vector 1, displacement is minus F_u
          ResU(3*ii-2) = ResU(3*ii-2) - b1*CauchySS(1) 
     &                  - b2*CauchySS(2) - b3*CauchySS(3)
          ResU(3*ii-1) = ResU(3*ii-1) - b1*CauchySS(4) 
     &                  - b2*CauchySS(5) - b3*CauchySS(6)
          ResU(3*ii)   = ResU(3*ii)   - b1*CauchySS(7) 
     &                  - b2*CauchySS(8) - b3*CauchySS(9)

        ! pressure Res_phi is minus F_phi
          Res_phi(ii) = Res_phi(ii) 
     &                 - dNphi_dxx(ii)*dvol*Mu0*(MCur(1)-dphi(1))
     &                 - dNphi_dyy(ii)*dvol*Mu0*(MCur(2)-dphi(2))
     &                 - dNphi_dzz(ii)*dvol*Mu0*(MCur(3)-dphi(3))

        END DO ! end of summation over ii=1, nNode --------------------

        ! residual vector 3, pressure ResP
        DO ii=1, nNode_p
          ResP(ii) = ResP(ii) - Np(ii)*(DetA -1.0 -Pressure/kappa)
     &               *dvol0*DetG
        END DO

        ! stiffness matrix 5 Kpp, summation over gauss points 8×8
        DO ii=1, nNode_p
          DO jj=1, nNode_p
            Kpp(ii,jj) = Kpp(ii,jj) - Np(ii)*Np(jj)/kappa*dvol0*DetG
          END DO
        END DO

      END DO ! end of loop over Gauss point gp=1, nGP

      ! stiffness matrix 6 K_phi_u, is the tranpose of K_u_phi 20×60
      DO ii=1, nNode
        DO jj=1, nNode
          K_phi_u(jj,3*ii-2) = K_u_phi_Tmp(3*ii-2,jj)      
          K_phi_u(jj,3*ii-1) = K_u_phi_Tmp(3*ii-1,jj)      
          K_phi_u(jj,3*ii)   = K_u_phi_Tmp(3*ii,jj)  
        END DO
      END DO

      ! stiffness matrix 7 Kpu, is the transpose of Kup 8×60
      DO ii=1, nNode
        DO jj=1, nNode_p
          Kpu(jj,3*ii-2) = Kup(3*ii-2,jj) 
          Kpu(jj,3*ii-1) = Kup(3*ii-1,jj) 
          Kpu(jj,3*ii)   = Kup(3*ii,jj) 
        END DO
      END DO

      ! WRITE(*,*) "#Elem in UEL", JELEM

! assemble Kuu to AMATRIX ---------------------------------------------
      ! WRITE(*,*) "assemble Kuu to AMATRIX"
      ! for row 1 to 8
      DO ii=1, nDofp
        DO jj=1, nDofp
          AMATRX(5*ii-4,5*jj-4) = Kuu(3*ii-2,3*jj-2)
          AMATRX(5*ii-4,5*jj-3) = Kuu(3*ii-2,3*jj-1)
          AMATRX(5*ii-4,5*jj-2) = Kuu(3*ii-2,3*jj)
          AMATRX(5*ii-4,5*jj-1) = K_u_phi(3*ii-2,jj) ! phi term
          AMATRX(5*ii-4,5*jj)   = Kup(3*ii-2,jj)

          AMATRX(5*ii-3,5*jj-4) = Kuu(3*ii-1,3*jj-2)
          AMATRX(5*ii-3,5*jj-3) = Kuu(3*ii-1,3*jj-1)
          AMATRX(5*ii-3,5*jj-2) = Kuu(3*ii-1,3*jj)
          AMATRX(5*ii-3,5*jj-1) = K_u_phi(3*ii-1,jj) ! phi term
          AMATRX(5*ii-3,5*jj)   = Kup(3*ii-1,jj)

          AMATRX(5*ii-2,5*jj-4) = Kuu(3*ii,3*jj-2)
          AMATRX(5*ii-2,5*jj-3) = Kuu(3*ii,3*jj-1)
          AMATRX(5*ii-2,5*jj-2) = Kuu(3*ii,3*jj)
          AMATRX(5*ii-2,5*jj-1) = K_u_phi(3*ii,jj) ! phi term
          AMATRX(5*ii-2,5*jj)   = Kup(3*ii,jj)

          AMATRX(5*ii-1,5*jj-4) = K_phi_u(ii,3*jj-2) ! phi term
          AMATRX(5*ii-1,5*jj-3) = K_phi_u(ii,3*jj-1) ! phi term
          AMATRX(5*ii-1,5*jj-2) = K_phi_u(ii,3*jj) ! phi term
          AMATRX(5*ii-1,5*jj-1) = K_phi_phi(ii,jj) ! phi term
          AMATRX(5*ii-1,5*jj)   = 0

          AMATRX(5*ii,5*jj-4) = Kpu(ii,3*jj-2)
          AMATRX(5*ii,5*jj-3) = Kpu(ii,3*jj-1)
          AMATRX(5*ii,5*jj-2) = Kpu(ii,3*jj)
          AMATRX(5*ii,5*jj-1) = 0
          AMATRX(5*ii,5*jj)   = Kpp(ii,jj)

!           WRITE(*,*) "Kuu, Kup", Kuu(3*ii-2,3*jj-2), 
!      &    Kuu(3*ii-2,3*jj-1), Kuu(3*ii-2,3*jj), Kup(3*ii-2,jj)

!           WRITE(*,*) "Kuu, Kup", Kuu(3*ii-1,3*jj-2), 
!      &    Kuu(3*ii-1,3*jj-1), Kuu(3*ii-1,3*jj), Kup(3*ii-1,jj)

!           WRITE(*,*) "Kuu, Kup", Kuu(3*ii,3*jj-2), 
!      &    Kuu(3*ii,3*jj-1), Kuu(3*ii,3*jj), Kup(3*ii,jj)

!           WRITE(*,*) "Kpu, Kpp", Kpu(ii,3*jj-2),
!      &   Kpu(ii,3*jj-1), Kpu(ii,3*jj), Kpp(ii,jj)
        END DO 
      END DO

      DO ii=1, nDofp
        DO jj=1+nDofp,nNode
          AMATRX(5*ii-4,nDofp+4*jj-3) = Kuu(3*ii-2,3*jj-2)
          AMATRX(5*ii-4,nDofp+4*jj-2) = Kuu(3*ii-2,3*jj-1)
          AMATRX(5*ii-4,nDofp+4*jj-1) = Kuu(3*ii-2,3*jj)
          AMATRX(5*ii-4,nDofp+4*jj)   = K_u_phi(3*ii-2,jj) ! phi term

          AMATRX(5*ii-3,nDofp+4*jj-3) = Kuu(3*ii-1,3*jj-2)
          AMATRX(5*ii-3,nDofp+4*jj-2) = Kuu(3*ii-1,3*jj-1)
          AMATRX(5*ii-3,nDofp+4*jj-1) = Kuu(3*ii-1,3*jj)
          AMATRX(5*ii-3,nDofp+4*jj)   = K_u_phi(3*ii-1,jj) ! phi term

          AMATRX(5*ii-2,nDofp+4*jj-3) = Kuu(3*ii,3*jj-2)
          AMATRX(5*ii-2,nDofp+4*jj-2) = Kuu(3*ii,3*jj-1)
          AMATRX(5*ii-2,nDofp+4*jj-1) = Kuu(3*ii,3*jj)
          AMATRX(5*ii-2,nDofp+4*jj)   = K_u_phi(3*ii,jj) ! phi term

          AMATRX(5*ii-1,nDofp+4*jj-3) = K_phi_u(ii,3*jj-2) ! phi term
          AMATRX(5*ii-1,nDofp+4*jj-2) = K_phi_u(ii,3*jj-1) ! phi term
          AMATRX(5*ii-1,nDofp+4*jj-1) = K_phi_u(ii,3*jj) ! phi term
          AMATRX(5*ii-1,nDofp+4*jj)   = K_phi_phi(ii,jj) ! phi term

          AMATRX(5*ii,nDofp+4*jj-3) = Kpu(ii,3*jj-2)
          AMATRX(5*ii,nDofp+4*jj-2) = Kpu(ii,3*jj-1)
          AMATRX(5*ii,nDofp+4*jj-1) = Kpu(ii,3*jj)
          AMATRX(5*ii,nDofp+4*jj)   = 0
        END DO 
      END DO

      ! for row 9 to 20
      DO ii=1+nDofp,nNode
        DO jj=1, nDofp
          AMATRX(nDofp+4*ii-3,5*jj-4) = Kuu(3*ii-2,3*jj-2)
          AMATRX(nDofp+4*ii-3,5*jj-3) = Kuu(3*ii-2,3*jj-1)
          AMATRX(nDofp+4*ii-3,5*jj-2) = Kuu(3*ii-2,3*jj)
          AMATRX(nDofp+4*ii-3,5*jj-1) = K_u_phi(3*ii-2,jj) ! phi term
          AMATRX(nDofp+4*ii-3,5*jj)   = Kup(3*ii-2,jj)

          AMATRX(nDofp+4*ii-2,5*jj-4) = Kuu(3*ii-1,3*jj-2)
          AMATRX(nDofp+4*ii-2,5*jj-3) = Kuu(3*ii-1,3*jj-1)
          AMATRX(nDofp+4*ii-2,5*jj-2) = Kuu(3*ii-1,3*jj)
          AMATRX(nDofp+4*ii-2,5*jj-1) = K_u_phi(3*ii-1,jj) ! phi term
          AMATRX(nDofp+4*ii-2,5*jj)   = Kup(3*ii-1,jj)

          AMATRX(nDofp+4*ii-1,5*jj-4) = Kuu(3*ii,3*jj-2)
          AMATRX(nDofp+4*ii-1,5*jj-3) = Kuu(3*ii,3*jj-1)
          AMATRX(nDofp+4*ii-1,5*jj-2) = Kuu(3*ii,3*jj)
          AMATRX(nDofp+4*ii-1,5*jj-1) = K_u_phi(3*ii,jj) ! phi term
          AMATRX(nDofp+4*ii-1,5*jj)   = Kup(3*ii,jj)

          AMATRX(nDofp+4*ii,5*jj-4) = K_phi_u(ii,3*jj-2) ! phi term
          AMATRX(nDofp+4*ii,5*jj-3) = K_phi_u(ii,3*jj-1) ! phi term
          AMATRX(nDofp+4*ii,5*jj-2) = K_phi_u(ii,3*jj) ! phi term
          AMATRX(nDofp+4*ii,5*jj-1) = K_phi_phi(ii,jj) ! phi term
          AMATRX(nDofp+4*ii,5*jj)   = 0
        END DO 
      END DO

      DO ii=1+nDofp,nNode
        DO jj=1+nDofp,nNode
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-3) = Kuu(3*ii-2,3*jj-2)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-2) = Kuu(3*ii-2,3*jj-1)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj-1) = Kuu(3*ii-2,3*jj)
          AMATRX(nDofp+4*ii-3,nDofp+4*jj)   = K_u_phi(3*ii-2,jj) ! phi term

          AMATRX(nDofp+4*ii-2,nDofp+4*jj-3) = Kuu(3*ii-1,3*jj-2)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj-2) = Kuu(3*ii-1,3*jj-1)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj-1) = Kuu(3*ii-1,3*jj)
          AMATRX(nDofp+4*ii-2,nDofp+4*jj)   = K_u_phi(3*ii-1,jj) ! phi term

          AMATRX(nDofp+4*ii-1,nDofp+4*jj-3) = Kuu(3*ii,3*jj-2)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj-2) = Kuu(3*ii,3*jj-1)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj-1) = Kuu(3*ii,3*jj)
          AMATRX(nDofp+4*ii-1,nDofp+4*jj)   = K_u_phi(3*ii,jj) ! phi term

          AMATRX(nDofp+4*ii,nDofp+4*jj-3) = K_phi_u(ii,3*jj-2) ! phi term
          AMATRX(nDofp+4*ii,nDofp+4*jj-2) = K_phi_u(ii,3*jj-1) ! phi term
          AMATRX(nDofp+4*ii,nDofp+4*jj-1) = K_phi_u(ii,3*jj) ! phi term
          AMATRX(nDofp+4*ii,nDofp+4*jj)   = K_phi_phi(ii,jj) ! phi term
        END DO 
      END DO

! assemble ResU and ResP to RHS ---------------------------------------
      ! WRITE(*,*) "ResU and ResP to RHS"
      ! for node couple with pressure, node 1 to 8
      kk = 0
      DO ii=1, nDofp
        RHS(5*ii-4,1) =  ResU(3*ii-2) 
        RHS(5*ii-3,1) =  ResU(3*ii-1) 
        RHS(5*ii-2,1) =  ResU(3*ii) 
        RHS(5*ii-1,1) =  Res_phi(ii) ! phi term
        RHS(5*ii  ,1) =  ResP(ii) 
        kk = kk+5
    !     WRITE(*,*) "#ii ResU ResP ", ii, 
    !  &             ResU(3*ii-2), ResU(3*ii-2), ResU(3*ii), ResP(ii)
      END DO
      ! WRITE(*,*) "kk =", kk

      ! for node only have displacement node 9 to 20
      DO ii=1+nDofp, nNode
        RHS(nDofp+4*ii-3,1) = ResU(3*ii-2)
        RHS(nDofp+4*ii-2,1) = ResU(3*ii-1)
        RHS(nDofp+4*ii-1,1) = ResU(3*ii) 
        RHS(nDofp+4*ii  ,1) = Res_phi(ii) ! phi term
        kk = kk+4
    !     WRITE(*,*) "#ii ResU ResP ", ii, 
    !  &             ResU(3*ii-2), ResU(3*ii-2), ResU(3*ii)
      END DO
      ! WRITE(*,*) "kk =", kk

! write ---------------------------------------------------------------
      ! check displacement 
      ! WRITE(*,*) "displacement"
      ! DO ii=1, nNode
      !   WRITE(*,*) dispElem(3*ii-2),dispElem(3*ii-1),dispElem(3*ii)
      ! END DO

      ! check Kuu 
      ! DO ii=1,nDofDisp
      !   DO jj=1,nDofDisp
      !     WRITE(*,*) "Kuu", ii,jj, Kuu(ii,jj)
      !   END DO
      ! END DO

      ! check Kup
      ! DO ii=1,nDofDisp
      !   DO jj=1,nDofp
      !     WRITE(*,*) "Kup",ii,jj, Kup(ii,jj)
      !   END DO
      ! END DO

      ! check K_u_phi
      ! DO ii=1,nDofDisp
      !   DO jj=1,nDofphi
      !     WRITE(*,*) "K_u_phi", ii,jj, K_u_phi(ii,jj)
      !   END DO
      ! END DO

      ! check K_phi_phi
      ! DO ii=1,nDofphi
      !   DO jj=1,nDofphi
      !     WRITE(*,*) "K_phi_phi", ii,jj, K_phi_phi(ii,jj)
      !   END DO
      ! END DO

      ! check ResU 
      ! DO ii=1,nDofDisp
      !   WRITE(*,*) "ResU", ii, ResU(ii)
      ! END DO

      ! check ResP 
      ! DO ii=1,nDofp
      !   WRITE(*,*) "ResP", ii, ResP(ii)
      ! END DO

      ! check Res_phi 
      ! DO ii=1,nDofphi
      !   WRITE(*,*) "Res_phi", ii, Res_phi(ii)
      ! END DO

      RETURN
      END SUBROUTINE UEL

! ------------------------ MagneticityTensor -------------------------- 
      SUBROUTINE MagneticityTensor(ElecTensor, Mu0)
      DOUBLE PRECISION :: ElecTensor(3,3), Mu0
      INTEGER :: ii, jj

      ! -delta_ij
      DO ii=1, 3
        DO jj=1, 3
          IF (ii==jj) THEN
           ElecTensor(ii,jj) = ElecTensor(ii,jj) 
     &         - Mu0
          END IF
        END DO
      END DO

      END SUBROUTINE MagneticityTensor
! ---------------------- CouplingTensor(u phi) ------------------------ 
      SUBROUTINE CouplingTensor1(CoupTensor, Mu0, dphi, DetF, F, MCur)
      DOUBLE PRECISION :: CoupTensor(9,3), Mu0, dphi(3), DetF, F(3,3),
     &                    MCur(3)
      INTEGER :: ii, jj, kk, ll, nn

      ! delta_ik = M_j
    !   DO ii=1, 3
    !     DO jj=1, 3
    !       DO kk=1, 3
    !         IF (ii==kk) THEN
    !          CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
    !  &          + Mu0*MCur(jj)
    !         END IF
    !       END DO
    !     END DO
    !   END DO

      ! dphi_k delta_ij
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

      ! dphi_j delta_ik
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

      ! dphi_i delta_kj
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

! ---------------------- CouplingTensor(phi u) ------------------------ 
      SUBROUTINE CouplingTensor2(CoupTensor, Mu0, dphi, DetF, F, MCur)
      DOUBLE PRECISION :: CoupTensor(9,3), Mu0, dphi(3), DetF, F(3,3),
     &                    MCur(3)
      INTEGER :: ii, jj, kk, ll, nn

      ! delta_ik = M_j
    !   DO ii=1, 3
    !     DO jj=1, 3
    !       DO kk=1, 3
    !         IF (ii==kk) THEN
    !          CoupTensor(3*(ii-1)+jj,kk) = CoupTensor(3*(ii-1)+jj,kk) 
    !  &          + Mu0*MCur(jj)
    !         END IF
    !       END DO
    !     END DO
    !   END DO

      ! dphi_k delta_ij
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

      ! dphi_j delta_ik
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

      ! dphi_i delta_kj
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

! ------------------------- ElasticityTensor --------------------------
      SUBROUTINE ElasticityTensor(ElsTensor, intPressure, 
     1     Mu, kappa, DetF, Ib, b, Mu0, dphi, MCur)
      
      DOUBLE PRECISION :: ElsTensor(9,9), b(3,3), dphi(3), MCur(3)
      INTEGER :: ii, jj, kk, ll, nn
      DOUBLE PRECISION :: Mu, kappa, Mu0
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi

      ! mechanical terms ----------------------------------------------
        ! b_kl delta_ji
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0/3.0*Mu*DetF**(-5.0/3.0)*b(kk,ll)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_lk delta_ji Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj .AND. ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0/3.0*Mu*DetF**(-5.0/3.0)*(-1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_ik b_jl
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*b(jj,ll)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! b_ij delta_lk
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*(-2.0/3.0)*b(ii,jj)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_jk delta_li Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (jj==kk .AND. ll==ii) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*(1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_ji delta_lk Pressure
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

        ! delta_il delta_jk Pressure
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

      ! coupling terms ------------------------------------------------
      Sq_dphi = dphi(1)*dphi(1) + dphi(2)*dphi(2) + dphi(3)*dphi(3)

        ! dphi_n dphi_n delta_jk delta_il
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
                IF (jj==kk .AND. ii==ll) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               + 0.5*Mu0*Sq_dphi
                END IF
              END DO
            END DO
          END DO
        END DO

        ! dphi_n dphi_n delta_ij delta_kl
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
                IF (ii==jj .AND. ll==kk) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               - 0.5*Mu0*Sq_dphi
                END IF
              END DO
            END DO
          END DO
        END DO

        ! dphi_k dphi_l delta_ij
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

        ! dphi_j dphi_k delta_il
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

        ! dphi_i dphi_j delta_lk
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

        ! dphi_i dphi_l delta_kj
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

        ! dphi_i dphi_k delta_jl
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

!         DO ii=1,9
!           WRITE(*,*) "ElasticityTensor",ii, ElsTensor(ii,1),
!      &    ElsTensor(ii,2),ElsTensor(ii,3),
!      &    ElsTensor(ii,4),ElsTensor(ii,5),
!      &    ElsTensor(ii,6),ElsTensor(ii,7),
!      &    ElsTensor(ii,8),ElsTensor(ii,9)
!         END DO

      END SUBROUTINE ElasticityTensor

! ----------------------- ElasticityTensorGent ------------------------
      SUBROUTINE ElasticityTensorGent(ElsTensor, intPressure, 
     1     Mu, kappa, DetF, Ib, b, Mu0, dphi, Im)
      
      DOUBLE PRECISION :: ElsTensor(9,9), b(3,3), dphi(3)
      INTEGER :: ii, jj, kk, ll, nn
      DOUBLE PRECISION :: Mu, kappa, Mu0, Im
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi, Quan1, Quan2

      Quan1 = (1.0-((DetF**(-2.0/3.0)*Ib)-3.0)/Im)**(-1.0)
      Quan2 = (1.0-((DetF**(-2.0/3.0)*Ib)-3.0)/Im)**(-2.0)

      ! mechanical dev terms ------------------------------------------
        ! b_ij delta_kl
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (kk==ll) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0/3.0*Mu*DetF**(-5.0/3.0)*Quan1*b(ii,jj)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_lk delta_ji Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj .AND. ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            -2.0/3.0*Mu*DetF**(-5.0/3.0)*Quan1*(-1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO

        ! b_kl b_ij
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + 2.0*Mu*(DetF**(-7.0/3.0))/Im*Quan2
     &            *b(kk,ll)*b(ii,jj)
              END DO
            END DO
          END DO
        END DO

        ! b_kl delta_ij Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + 2.0*Mu*(DetF**(-7.0/3.0))/Im*Quan2
     &            *b(kk,ll)*(-1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO


        ! b_ij delta_kl Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (kk==ll) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + 2.0*Mu*(DetF**(-7.0/3.0))/Im*Quan2
     &            *b(ii,jj)*(-1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_ij delta_kl Ib Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (jj==ii .AND. ll==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + 2.0*Mu*(DetF**(-7.0/3.0))/Im*Quan2
     &            *(-1.0/3.0)*Ib*(-1.0/3.0)*Ib
               END IF
              END DO
            END DO
          END DO
        END DO

        ! b_jl delta_ik 
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==kk) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*Quan1*b(jj,ll)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! b_kl delta_ij 
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (ii==jj) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*Quan1*b(kk,ll)*(-2.0/3.0)
               END IF
              END DO
            END DO
          END DO
        END DO

        ! delta_kj delta_il Ib
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
               IF (kk==jj .AND. ii==ll) THEN
                ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &            ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &            + Mu*DetF**(-5.0/3.0)*Quan1*Ib*(1.0/3.0)
               END IF
              END DO
            END DO
          END DO
        END DO

      ! mechanical vol terms ------------------------------------------
        ! delta_ji delta_lk Pressure
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

        ! delta_il delta_jk Pressure
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

      ! coupling terms ------------------------------------------------
      Sq_dphi = dphi(1)*dphi(1) + dphi(2)*dphi(2) + dphi(3)*dphi(3)

        ! dphi_n dphi_n delta_jk delta_il
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
                IF (jj==kk .AND. ii==ll) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               + 0.5*Mu0*Sq_dphi
                END IF
              END DO
            END DO
          END DO
        END DO

        ! dphi_n dphi_n delta_ij delta_lk
        DO ii=1, 3
          DO jj=1, 3
            DO kk=1, 3
              DO ll=1, 3
                IF (ii==jj .AND. ll==kk) THEN
                  ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) = 
     &               ElsTensor(3*(jj-1)+ii,3*(ll-1)+kk) 
     &               - 0.5*Mu0*Sq_dphi
                END IF
              END DO
            END DO
          END DO
        END DO

        ! dphi_k dphi_l delta_ij
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

        ! dphi_j dphi_k delta_il
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

        ! dphi_i dphi_j delta_lk
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

        ! dphi_i dphi_l delta_kj
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

        ! dphi_i dphi_k delta_lj
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

      END SUBROUTINE ElasticityTensorGent

! ------------------------- cal_CauchyStress --------------------------
      SUBROUTINE cal_CauchyStress(CauchySS, intPressure, 
     1     Mu, kappa, DetF, Ib, b, Mu0, dphi, F, MCur, Ha)
      
      DOUBLE PRECISION :: b(3,3), CauchySS(9), dphi(3), F(3,3), MCur(3)
      INTEGER :: ii, jj, kk, ll, PP, QQ
      DOUBLE PRECISION :: Mu, kappa, Mu0, Ha(3)
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi, Elast1

      Sq_dphi = dphi(1)*dphi(1) + dphi(2)*dphi(2) + dphi(3)*dphi(3)
      Elast1 = Mu*DetF**(-5.0/3.0)

  ! Calculate my effective Cauchy stress {11,21,31,12,22,32,13,23,33}
      CauchySS(1) = Elast1*(b(1,1) - Ib/3.0) + intPressure
     &            + Mu0*(dphi(1)*dphi(1)- 0.5*Sq_dphi)
     &            - Mu0*MCur(1)*Ha(1)
      CauchySS(2) = Elast1*b(2,1) 
     &            + Mu0*(dphi(2)*dphi(1))
     &            - Mu0*MCur(2)*Ha(1)
      CauchySS(3) = Elast1*b(3,1)
     &            + Mu0*(dphi(3)*dphi(1))
     &            - Mu0*MCur(3)*Ha(1)

      CauchySS(4) = Elast1*b(1,2)
     &            + Mu0*(dphi(1)*dphi(2))
     &            - Mu0*MCur(1)*Ha(2)
      CauchySS(5) = Elast1*(b(2,2) - Ib/3.0) + intPressure
     &            + Mu0*(dphi(2)*dphi(2)- 0.5*Sq_dphi)
     &            - Mu0*MCur(2)*Ha(2)
      CauchySS(6) = Elast1*b(3,2)
     &            + Mu0*(dphi(3)*dphi(2))
     &            - Mu0*MCur(3)*Ha(2)

      CauchySS(7) = Elast1*b(1,3)
     &            + Mu0*(dphi(1)*dphi(3))
     &            - Mu0*MCur(1)*Ha(3)
      CauchySS(8) = Elast1*b(2,3)
     &            + Mu0*(dphi(2)*dphi(3))
     &            - Mu0*MCur(2)*Ha(3)
      CauchySS(9) = Elast1*(b(3,3) - Ib/3.0) + intPressure
     &            + Mu0*(dphi(3)*dphi(3)- 0.5*Sq_dphi)
     &            - Mu0*MCur(3)*Ha(3)

      ! write Cauchy stress for check
    !   DO pp=1,3
    !   WRITE(*,*) "CauchySS", CauchySS(3*pp-2), CauchySS(3*pp-1), 
    !  &    CauchySS(3*pp)
    !   END DO

      END SUBROUTINE cal_CauchyStress

! ----------------------- cal_CauchyStressGent ------------------------
      SUBROUTINE cal_CauchyStressGent(CauchySS, intPressure, 
     1     Mu, kappa, DetF, Ib, b, Mu0, dphi, Im)
      
      DOUBLE PRECISION :: b(3,3), CauchySS(9), dphi(3)
      INTEGER :: ii, jj, kk, ll, PP, QQ
      DOUBLE PRECISION :: Mu, kappa, Mu0, Im
      DOUBLE PRECISION :: DetF, Ib, intPressure, Sq_dphi, Elast1

      Sq_dphi = dphi(1)*dphi(1) + dphi(2)*dphi(2) + dphi(3)*dphi(3)
      Elast1 = Mu*DetF**(-5.0/3.0)
     &       *((1.0-((DetF**(-2.0/3.0)*Ib)-3.0)/Im)**(-1.0))

  ! Calculate my effective Cauchy stress {11,21,31,12,22,32,13,23,33}
      CauchySS(1) = Elast1*(b(1,1) - Ib/3.0)
     &   + intPressure + Mu0*(dphi(1)*dphi(1) - 0.5*Sq_dphi)
      CauchySS(2) = Elast1*b(1,2) 
     &   + Mu0*(dphi(1)*dphi(2))
      CauchySS(3) = Elast1*b(1,3)
     &   + Mu0*(dphi(1)*dphi(3))

      CauchySS(4) = Elast1*b(2,1)
     &   + Mu0*(dphi(2)*dphi(1))
      CauchySS(5) = Elast1*(b(2,2) - Ib/3.0)
     &   + intPressure + Mu0*(dphi(2)*dphi(2) - 0.5*Sq_dphi)
      CauchySS(6) = Elast1*b(2,3)
     &   + Mu0*(dphi(2)*dphi(3))

      CauchySS(7) = Elast1*b(3,1)
     &   + Mu0*(dphi(3)*dphi(1))
      CauchySS(8) = Elast1*b(3,2)
     &   + Mu0*(dphi(3)*dphi(2))
      CauchySS(9) = Elast1*(b(3,3) - Ib/3.0)
     &   + intPressure + Mu0*(dphi(3)*dphi(3) - 0.5*Sq_dphi)

      END SUBROUTINE cal_CauchyStressGent

! ----------------------------- uvarm --------------------------------
      subroutine uvarm(uvar,direct,t,time,dtime,cmname,orname,
     1 nuvarm,noel,npt,layer,kspt,kstep,kinc,ndi,nshr,coord,
     2 jmac,jmatyp,matlayo,laccfla)

      ! this subroutine is used to transfer SDV's from the UEL
      !  onto the dummy mesh for viewing.  Note that an offset of
      !  ElemOffset is used between the real mesh and the dummy mesh.

      !  this will need to be modified.

      use global
      include 'aba_param.inc'

      character*80 cmname,orname
      character*3 flgray(15)
      dimension uvar(nuvarm),direct(3,3),t(3,3),time(2)
      dimension array(15),jarray(15),jmac(*),jmatyp(*),coord(*)

C     The dimensions of the variables FLGRAY, ARRAY and JARRAY
C     must be set equal to or greater than 15.

      NELEMAN = NOEL-ELEMOFFSET

      do i=1,nUvrm
         uvar(i) = GloUVar(NELEMAN,NPT,i)
      end do

      return
      end subroutine uvarm

! ------------------------ GetGaussPoints1D --------------------------
      SUBROUTINE GetGaussPoints1D(nGP, GaussPoints, GaussWeights)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nGP
      DOUBLE PRECISION :: GaussPoints(nGP), GaussWeights(nGP)

      SELECT CASE(nGP)
        CASE(1)  ! 1 Point quadrature rule
          GaussPoints(1) = 0.0
          GaussWeights(1) = 2.0
        CASE(2)  ! 2 Point quadrature rule
          GaussPoints(1) = -0.577350269189626
          GaussWeights(1) = 1.0
          GaussPoints(2) =  0.577350269189626
          GaussWeights(2) = 1.0
        CASE(3)  ! 3 Point quadrature rule
          GaussPoints(1) = -0.774596669241483
          GaussWeights(1) = 0.555555555555556
          GaussPoints(2) =  0.0
          GaussWeights(2) = 0.888888888888889
          GaussPoints(3) =  0.774596669241483
          GaussWeights(3) = 0.555555555555556
        CASE(4)  ! 4 Point quadrature rule
          GaussPoints(1) = -0.861136311594953
          GaussWeights(1) = 0.347854845137454
          GaussPoints(2) = -0.339981043584856
          GaussWeights(2) = 0.652145154862546
          GaussPoints(3) =  0.339981043584856
          GaussWeights(3) = 0.652145154862546
          GaussPoints(4) =  0.861136311594953
          GaussWeights(4) = 0.347854845137454
        CASE(5) ! 5 Point quadrature rule
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
        CASE(6) ! 6 Point quadrature rule
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
        CASE(7) ! 7 Point quadrature rule
          GaussPoints(1) = -0.9491079123427585245261897
          GaussWeights(1) = 0.1294849661688696932706114
          GaussPoints(2) = -0.7415311855993944398638648
          GaussWeights(2) = 0.2797053914892766679014678
          GaussPoints(3) = -0.4058451513773971669066064
          GaussWeights(3) = 0.3818300505051189449503698
          GaussPoints(4) =  0.0			        
          GaussWeights(4) = 0.4179591836734693877551020
          GaussPoints(5) =  0.4058451513773971669066064
          GaussWeights(5) = 0.3818300505051189449503698
          GaussPoints(6) =  0.7415311855993944398638648
          GaussWeights(6) = 0.2797053914892766679014678
          GaussPoints(7) =  0.9491079123427585245261897
          GaussWeights(7) = 0.1294849661688696932706114
        CASE(8) ! 8 Point quadrature rule
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
          WRITE(*,*) " invalid value of 'nGP' ! "
      END SELECT

      END SUBROUTINE GetGaussPoints1D

! ----------------------- getGaussPointsHexa -------------------------
      SUBROUTINE getGaussPointsHexa(nGP, gpts1, gpts2, gpts3, gwts)
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: nGP
      DOUBLE PRECISION :: gpts1(nGP), gpts2(nGP), gpts3(nGP), 
     & gwts(nGP)

      INTEGER ::  nn, ii, jj, kk, ll

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gpoints1, 
     &  gweights1
      ! WRITE(*,*) "Enter getGaussPointsHexa..."
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

! ---------------------- LagrangeBasisFunsTet -----------------------
      SUBROUTINE LagrangeBasisFunsTet(
     1  degree, 
     2  xi1, xi2, xi3,
     3  N, dN_dxi1, dN_dxi2, dN_dxi3)

      INTEGER :: degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi3
      DOUBLE PRECISION :: fact1, fact2, val1, val2, val3, val4
  
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

          dN_dxi1(1) =  1.0
          dN_dxi1(2) =  0.0
          dN_dxi1(4) =  0.0
          dN_dxi1(3) = -1.0

          dN_dxi2(1) =  0.0
          dN_dxi2(2) =  1.0
          dN_dxi2(4) =  0.0
          dN_dxi2(3) = -1.0

          dN_dxi3(1) =  0.0
          dN_dxi3(2) =  0.0
          dN_dxi3(4) =  1.0
          dN_dxi3(3) = -1.0

      CASE DEFAULT

          WRITE(*,*) "no basis functions for degree", degree

      END SELECT

      END SUBROUTINE LagrangeBasisFunsTet

! --------------------- LagrangeBasisFunsPrism -----------------------

      SUBROUTINE LagrangeBasisFunsPrism(degree, xi1, xi2, xi4, 
     &  N, dN_dxi1, dN_dxi2, dN_dxi4)
      !xi1, xi2 and xi3 are the parametric coordinates for the triangle, and
      !xi4 is the parametric coordinate in the direction normal to the triangle

      IMPLICIT NONE

      INTEGER :: degree
      DOUBLE PRECISION, INTENT(IN) :: xi1, xi2, xi4
      DOUBLE PRECISION, DIMENSION(:) :: N, dN_dxi1, dN_dxi2, dN_dxi4
      DOUBLE PRECISION :: xi3

      xi3 = 1.0- xi1 - xi2

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

          WRITE(*,*) "no basis functions defined for this degree IN LagrangeBasisFunsPrism"

      END SELECT

      END SUBROUTINE LagrangeBasisFunsPrism

! ---------------------- LagrangeBasisFunsHex -----------------------

      SUBROUTINE  LagrangeBasisFunsHex(nNode, degree, xi1, xi2, xi3, 
     &  N, dN_dxi1, dN_dxi2, dN_dxi3)

      IMPLICIT NONE

      INTEGER :: nNode, degree
      DOUBLE PRECISION :: xi1, xi2, xi3
      DOUBLE PRECISION :: N(nNode), dN_dxi1(nNode), dN_dxi2(nNode),
     &  dN_dxi3(nNode)

      DOUBLE PRECISION ::  v11, v12, v21, v22, v31, v32

      SELECT CASE (degree)
        case (0)

          N(1) = 1.0

          dN_dxi1(1) = 0.0
          dN_dxi2(1) = 0.0
          dN_dxi3(1) = 0.0

        case (1) ! linear hexahedron element

          v11 = 1.0 - xi1
          v12 = 1.0 + xi1
          v21 = 1.0 - xi2
          v22 = 1.0 + xi2
          v31 = 1.0 - xi3
          v32 = 1.0 + xi3
          ! change 3 <--> 4 and 7 <--> 8 for orientation
          N(1) = 0.125*v11*v21*v31
          N(2) = 0.125*v12*v21*v31
          N(4) = 0.125*v11*v22*v31
          N(3) = 0.125*v12*v22*v31
          N(5) = 0.125*v11*v21*v32
          N(6) = 0.125*v12*v21*v32
          N(8) = 0.125*v11*v22*v32
          N(7) = 0.125*v12*v22*v32

          dN_dxi1(1) = -0.125*v21*v31
          dN_dxi1(2) =  0.125*v21*v31
          dN_dxi1(4) = -0.125*v22*v31
          dN_dxi1(3) =  0.125*v22*v31
          dN_dxi1(5) = -0.125*v21*v32
          dN_dxi1(6) =  0.125*v21*v32
          dN_dxi1(8) = -0.125*v22*v32
          dN_dxi1(7) =  0.125*v22*v32

          dN_dxi2(1) = -0.125*v11*v31
          dN_dxi2(2) = -0.125*v12*v31
          dN_dxi2(4) =  0.125*v11*v31
          dN_dxi2(3) =  0.125*v12*v31
          dN_dxi2(5) = -0.125*v11*v32
          dN_dxi2(6) = -0.125*v12*v32
          dN_dxi2(8) =  0.125*v11*v32
          dN_dxi2(7) =  0.125*v12*v32

          dN_dxi3(1) = -0.125*v11*v21
          dN_dxi3(2) = -0.125*v12*v21
          dN_dxi3(4) = -0.125*v11*v22
          dN_dxi3(3) = -0.125*v12*v22
          dN_dxi3(5) =  0.125*v11*v21
          dN_dxi3(6) =  0.125*v12*v21
          dN_dxi3(8) =  0.125*v11*v22
          dN_dxi3(7) =  0.125*v12*v22

        case (2) ! serendipity hexahedron element

          N(1) =  1.0/8.0*(-1+xi1)*(-1+xi2)*(-1+xi3)*(2+xi1+xi2+xi3)
          N(2) = -1.0/8.0*(1+xi1)*(-1+xi2)*(-1+xi3)*(2-xi1+xi2+xi3)
          N(3) = -1.0/8.0*(1+xi1)*(1+xi2)*(-1+xi3)*(-2+xi1+xi2-xi3)
          N(4) = -1.0/8.0*(-1+xi1)*(1+xi2)*(-1+xi3)*(2+xi1-xi2+xi3)

          N(5) = -1.0/8.0*(-1+xi1)*(-1+xi2)*(2+xi1+xi2-xi3)*(1+xi3)
          N(6) = -1.0/8.0*(1+xi1)*(-1+xi2)*(1+xi3)*(-2+xi1-xi2+xi3)
          N(7) =  1.0/8.0*(1+xi1)*(1+xi2)*(1+xi3)*(-2+xi1+xi2+xi3)
          N(8) =  1.0/8.0*(-1+xi1)*(1+xi2)*(2+xi1-xi2-xi3)*(1+xi3)

          N(9)  = -1.0/4.0*(-1+xi1**2)*(-1+xi2)*(-1+xi3)
          N(10) =  1.0/4.0*(1+xi1)*(-1+xi2**2)*(-1+xi3)
          N(11) =  1.0/4.0*(-1+xi1**2)*(1+xi2)*(-1+xi3)
          N(12) = -1.0/4.0*(-1+xi1)*(-1+xi2**2)*(-1+xi3)
          
          N(13) =  1.0/4.0*(-1+xi1**2)*(-1+xi2)*(1+xi3)
          N(14) = -1.0/4.0*(1+xi1)*(-1+xi2**2)*(1+xi3)
          N(15) = -1.0/4.0*(-1+xi1**2)*(1+xi2)*(1+xi3)
          N(16) =  1.0/4.0*(-1+xi1)*(-1+xi2**2)*(1+xi3)
          
          N(17) = -1.0/4.0*(-1+xi1)*(-1+xi2)*(-1+xi3**2)
          N(18) =  1.0/4.0*(1+xi1)*(-1+xi2)*(-1+xi3**2)
          N(19) = -1.0/4.0*(1+xi1)*(1+xi2)*(-1+xi3**2)
          N(20) =  1.0/4.0*(-1+xi1)*(1+xi2)*(-1+xi3**2)

          ! dN_dxi1
          dN_dxi1(1) =  1.0/8.0*(-1+xi2)*(-1+xi3)*(1+2*xi1+xi2+xi3)
          dN_dxi1(2) = -1.0/8.0*(-1+xi2)*(-1+xi3)*(1-2*xi1+xi2+xi3)
          dN_dxi1(3) =  1.0/8.0*(1+xi2)*(-1+xi3)*(1-2*xi1-xi2+xi3)
          dN_dxi1(4) =  1.0/8.0*(1+xi2)*(-1-2*xi1+xi2-xi3)*(-1+xi3)

          dN_dxi1(5) = -1.0/8.0*(-1+xi2)*(1+2*xi1+xi2-xi3)*(1+xi3)
          dN_dxi1(6) =  1.0/8.0*(-1+xi2)*(1-2*xi1+xi2-xi3)*(1+xi3)
          dN_dxi1(7) =  1.0/8.0*(1+xi2)*(1+xi3)*(-1+2*xi1+xi2+xi3)
          dN_dxi1(8) = -1.0/8.0*(1+xi2)*(1+xi3)*(-1-2*xi1+xi2+xi3)

          dN_dxi1(9)  = -1.0/2.0*xi1*(-1+xi2)*(-1+xi3)
          dN_dxi1(10) =  1.0/4.0*(-1+xi2**2)*(-1+xi3)
          dN_dxi1(11) =  1.0/2.0*xi1*(1+xi2)*(-1+xi3)
          dN_dxi1(12) = -1.0/4.0*(-1+xi2**2)*(-1+xi3)

          dN_dxi1(13) =  1.0/2.0*xi1*(-1+xi2)*(1+xi3)
          dN_dxi1(14) = -1.0/4.0*(-1+xi2**2)*(1+xi3)
          dN_dxi1(15) = -1.0/2.0*xi1*(1+xi2)*(1+xi3)
          dN_dxi1(16) =  1.0/4.0*(-1+xi2**2)*(1+xi3)

          dN_dxi1(17) = -1.0/4.0*(-1+xi2)*(-1+xi3**2)
          dN_dxi1(18) =  1.0/4.0*(-1+xi2)*(-1+xi3**2)
          dN_dxi1(19) = -1.0/4.0*(1+xi2)*(-1+xi3**2)
          dN_dxi1(20) =  1.0/4.0*(1+xi2)*(-1+xi3**2)

          ! dN_dxi2
          dN_dxi2(1) =  1.0/8.0*(-1+xi1)*(-1+xi3)*(1+xi1+2*xi2+xi3)
          dN_dxi2(2) =  1.0/8.0*(1+xi1)*(-1+xi1-2*xi2-xi3)*(-1+xi3)
          dN_dxi2(3) =  1.0/8.0*(1+xi1)*(-1+xi3)*(1-xi1-2*xi2+xi3)
          dN_dxi2(4) = -1.0/8.0*(-1+xi1)*(-1+xi3)*(1+xi1-2*xi2+xi3)

          dN_dxi2(5) = -1.0/8.0*(-1+xi1)*(1+xi1+2*xi2-xi3)*(1+xi3)
          dN_dxi2(6) = -1.0/8.0*(1+xi1)*(1+xi3)*(-1+xi1-2*xi2+xi3)
          dN_dxi2(7) =  1.0/8.0*(1+xi1)*(1+xi3)*(-1+xi1+2*xi2+xi3)
          dN_dxi2(8) =  1.0/8.0*(-1+xi1)*(1+xi1-2*xi2-xi3)*(1+xi3)

          dN_dxi2(9)  = -1.0/4.0*(-1+xi1**2)*(-1+xi3)
          dN_dxi2(10) =  1.0/2.0*(1+xi1)*xi2*(-1+xi3)
          dN_dxi2(11) =  1.0/4.0*(-1+xi1**2)*(-1+xi3)
          dN_dxi2(12) = -1.0/2.0*(-1+xi1)*xi2*(-1+xi3)

          dN_dxi2(13) =  1.0/4.0*(-1+xi1**2)*(1+xi3)
          dN_dxi2(14) = -1.0/2.0*(1+xi1)*xi2*(1+xi3)
          dN_dxi2(15) = -1.0/4.0*(-1+xi1**2)*(1+xi3)
          dN_dxi2(16) =  1.0/2.0*(-1+xi1)*xi2*(1+xi3)

          dN_dxi2(17) = -1.0/4.0*(-1+xi1)*(-1+xi3**2)
          dN_dxi2(18) =  1.0/4.0*(1+xi1)*(-1+xi3**2)
          dN_dxi2(19) = -1.0/4.0*(1+xi1)*(-1+xi3**2)
          dN_dxi2(20) =  1.0/4.0*(-1+xi1)*(-1+xi3**2)

          ! dN_dxi3
          dN_dxi3(1) =  1.0/8.0*(-1+xi1)*(-1+xi2)*(1+xi1+xi2+2*xi3)
          dN_dxi3(2) =  1.0/8.0*(1+xi1)*(-1+xi2)*(-1+xi1-xi2-2*xi3)
          dN_dxi3(3) = -1.0/8.0*(1+xi1)*(1+xi2)*(-1+xi1+xi2-2*xi3)
          dN_dxi3(4) = -1.0/8.0*(-1+xi1)*(1+xi2)*(1+xi1-xi2+2*xi3)

          dN_dxi3(5) = -1.0/8.0*(-1+xi1)*(-1+xi2)*(1+xi1+xi2-2*xi3)
          dN_dxi3(6) =  1.0/8.0*(1+xi1)*(-1+xi2)*(1-xi1+xi2-2*xi3)
          dN_dxi3(7) =  1.0/8.0*(1+xi1)*(1+xi2)*(-1+xi1+xi2+2*xi3)
          dN_dxi3(8) =  1.0/8.0*(-1+xi1)*(1+xi2)*(1+xi1-xi2-2*xi3)

          dN_dxi3(9)  = -1.0/4.0*(-1+xi1**2)*(-1+xi2)
          dN_dxi3(10) =  1.0/4.0*(1+xi1)*(-1+xi2**2)
          dN_dxi3(11) =  1.0/4.0*(-1+xi1**2)*(1+xi2)
          dN_dxi3(12) = -1.0/4.0*(-1+xi1)*(-1+xi2**2)

          dN_dxi3(13) =  1.0/4.0*(-1+xi1**2)*(-1+xi2)
          dN_dxi3(14) = -1.0/4.0*(1+xi1)*(-1+xi2**2)
          dN_dxi3(15) = -1.0/4.0*(-1+xi1**2)*(1+xi2)
          dN_dxi3(16) =  1.0/4.0*(-1+xi1)*(-1+xi2**2)

          dN_dxi3(17) = -1.0/2.0*(-1+xi1)*(-1+xi2)*xi3
          dN_dxi3(18) =  1.0/2.0*(1+xi1)*(-1+xi2)*xi3
          dN_dxi3(19) = -1.0/2.0*(1+xi1)*(1+xi2)*xi3
          dN_dxi3(20) =  1.0/2.0*(-1+xi1)*(1+xi2)*xi3

        CASE DEFAULT
          WRITE(*,*) "degree is", degree
          WRITE(*,*) "no basis functions defined for this degree IN LagrangeBasisFunsHex"

      END SELECT

      END SUBROUTINE LagrangeBasisFunsHex

! --------------------- computeBasisFunctions3D ----------------------
      SUBROUTINE computeBasisFunctions3D(
     1  nNode, ETYPE, degree, 
     2  param, 
     3  xNode, yNode, zNode,
     4  N,
     5  dN_dx, dN_dy, dN_dz,
     6  Jac)

      INTEGER :: nNode, ETYPE, degree
      DOUBLE PRECISION :: param(3), Jac
      DOUBLE PRECISION :: xNode(nNode), yNode(nNode), zNode(nNode)
      DOUBLE PRECISION :: N(nNode), dN_dx(nNode), dN_dy(nNode),
     & dN_dz(nNode)

      INTEGER :: ii=0, jj=0, kk=0, count=0

      DOUBLE PRECISION :: dN_du1(nNode), dN_du2(nNode), dN_du3(nNode)
      DOUBLE PRECISION :: xx=0.0, yy=0.0, zz=0.0, detinv
      DOUBLE PRECISION :: B(3,3), Binv(3,3)

!         DO ii=1,8
!         WRITE(*,*) "xNode, yNode, zNode", xNode(ii), yNode(ii), 
!      &    zNode(ii)
!         END DO

      IF(ETYPE == 4) THEN !tetrahedral elements
        call LagrangeBasisFunsTet(
     1    degree, 
     2    param(1), param(2), param(3),
     3    N, dN_du1, dN_du2, dN_du3)
      ELSE IF(ETYPE == 6) THEN !penta elements
        call LagrangeBasisFunsPrism(
     1    degree,
     2    param(1), param(2), param(3),
     3    N, dN_du1, dN_du2, dN_du3)
      ELSE  !! hex elements
        ! WRITE(*,*) "Enter LagrangeBasisFunsHex..."
        call LagrangeBasisFunsHex(nNode, degree, 
     &    param(1), param(2), param(3), N, dN_du1, dN_du2, dN_du3)
      END IF

      !Gradient of mapping from parameter space to physical space
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

      ! WRITE(*,*) "transformation matrix B = ", B(1,1), B(1,2), B(1,3)
      ! WRITE(*,*) "transformation matrix B = ", B(2,1), B(2,2), B(2,3)
      ! WRITE(*,*) "transformation matrix B = ", B(3,1), B(3,2), B(3,3)

      ! determinant of the matrix
      Jac = -B(1,3)*B(2,2)*B(3,1) + B(1,2)*B(2,3)*B(3,1)
     1      +B(1,3)*B(2,1)*B(3,2) - B(1,1)*B(2,3)*B(3,2)
     2      -B(1,2)*B(2,1)*B(3,3) + B(1,1)*B(2,2)*B(3,3)

      ! inverse determinant of the matrix
      detinv = 1.0/Jac

      ! Calculate the inverse of the matrix
      Binv(1,1) = +detinv * (B(2,2)*B(3,3) - B(2,3)*B(3,2))
      Binv(2,1) = -detinv * (B(2,1)*B(3,3) - B(2,3)*B(3,1))
      Binv(3,1) = +detinv * (B(2,1)*B(3,2) - B(2,2)*B(3,1))
      Binv(1,2) = -detinv * (B(1,2)*B(3,3) - B(1,3)*B(3,2))
      Binv(2,2) = +detinv * (B(1,1)*B(3,3) - B(1,3)*B(3,1))
      Binv(3,2) = -detinv * (B(1,1)*B(3,2) - B(1,2)*B(3,1))
      Binv(1,3) = +detinv * (B(1,2)*B(2,3) - B(1,3)*B(2,2))
      Binv(2,3) = -detinv * (B(1,1)*B(2,3) - B(1,3)*B(2,1))
      Binv(3,3) = +detinv * (B(1,1)*B(2,2) - B(1,2)*B(2,1))

      !!Compute derivatives of basis functions w.r.t physical coordinates
      DO ii=1, nNode
        dN_dx(ii) = dN_du1(ii) * Binv(1,1) + dN_du2(ii) * Binv(1,2)
     &    + dN_du3(ii) * Binv(1,3)
        dN_dy(ii) = dN_du1(ii) * Binv(2,1) + dN_du2(ii) * Binv(2,2)
     &    + dN_du3(ii) * Binv(2,3)
        dN_dz(ii) = dN_du1(ii) * Binv(3,1) + dN_du2(ii) * Binv(3,2)
     &    + dN_du3(ii) * Binv(3,3)
      END DO

      END SUBROUTINE computeBasisFunctions3D

! ----------------------------- GrowthFun -----------------------------
      SUBROUTINE GrowthFun(GF_G, 
     1     GF_xCoord,GF_yCoord,GF_zCoord, GF_time, TotalT, gFactor)

      DOUBLE PRECISION :: GF_xCoord, GF_yCoord, GF_zCoord
      DOUBLE PRECISION :: GF_G(3,3)
      DOUBLE PRECISION :: Lambda1, Lambda2, Lambda3
      DOUBLE PRECISION :: Lambda1z0, Lambda1z1, Lambda2z0, Lambda2z1
      DOUBLE PRECISION :: TotalT, GF_time, Pi, gFactor
      INTEGER :: gp

      Pi=3.14159265359
!       Lambda1z0 = Pi*(2.5-1.5*Cos(2*Pi*GF_xCoord))**0.5
!       Lambda1z1 = 4*Pi/(5.0-3.0*Cos(2*Pi*GF_xCoord))
!       Lambda2z0 = 2*Pi*Sin(Pi*GF_xCoord)
!       Lambda2z1 = (4.0*Pi*2.0**0.5*Sin(Pi*GF_xCoord))/
!      &  ((5.0-3.0*Cos(2.0*Pi*GF_xCoord))**0.5)

      ! isotropic growth
      Lambda1z0 = gFactor
      Lambda1z1 = 0.0
      Lambda2z0 = gFactor
      Lambda2z1 = 0.0

      Lambda1 = Lambda1z0 + GF_zCoord*Lambda1z1
      Lambda2 = Lambda2z0 + GF_zCoord*Lambda2z1
      Lambda3 = gFactor

      GF_G(1,1) = 1.0+(Lambda1-1.0) * (GF_time)/TotalT
      GF_G(1,2) = (0.0) * (GF_time)/TotalT
      GF_G(1,3) = (0.0) * (GF_time)/TotalT

      GF_G(2,1) = (0.0) * (GF_time)/TotalT
      GF_G(2,2) = 1.0+(Lambda2-1.0) * (GF_time)/TotalT
      GF_G(2,3) = (0.0) * (GF_time)/TotalT

      GF_G(3,1) = (0.0) * (GF_time)/TotalT
      GF_G(3,2) = (0.0) * (GF_time)/TotalT
      GF_G(3,3) = 1.0+(Lambda3-1.0) * (GF_time)/TotalT

    !    WRITE(*,*) "G",GF_G(1,1),GF_G(1,2),GF_G(1,3),
    !   &  GF_G(2,1),GF_G(2,2),GF_G(2,3),
    !   &  GF_G(3,1),GF_G(3,2),GF_G(3,3)

      END SUBROUTINE GrowthFun


! ---------------------------- DetMatrix ------------------------------
      SUBROUTINE DetMatrix(M, DetM)

      DOUBLE PRECISION :: M(3,3), DetM

      DetM = -M(1,3)*M(2,2)*M(3,1) + M(1,2)*M(2,3)*M(3,1)
     &       +M(1,3)*M(2,1)*M(3,2) - M(1,1)*M(2,3)*M(3,2)
     &       -M(1,2)*M(2,1)*M(3,3) + M(1,1)*M(2,2)*M(3,3)

      END SUBROUTINE DetMatrix

! ---------------------------- InvMatrix ------------------------------
      SUBROUTINE InvMatrix(MInv, M, DetM)

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

! ------------------------------ CalAe --------------------------------
      SUBROUTINE CalAe(Ae, F, G)

      DOUBLE PRECISION :: Ae(3,3), F(3,3), G(3,3)

        Ae(1,1) = (F(1,3)*G(2,2)*G(3,1)-F(1,2)*G(2,3)*G(3,1)-F(1,3)*
     &    G(2,1)*G(3,2)+F(1,1)*G(2,3)*G(3,2)+F(1,2)*G(2,1)*G(3,3)-
     &    F(1,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(1,2) = (F(1,3)*G(1,2)*G(3,1)-F(1,2)*G(1,3)*G(3,1)-F(1,3)*
     &    G(1,1)*G(3,2)+F(1,1)*G(1,3)*G(3,2)+F(1,2)*G(1,1)*G(3,3)-
     &    F(1,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(1,3) = (F(1,3)*G(1,2)*G(2,1)-F(1,2)*G(1,3)*G(2,1)-F(1,3)*
     &    G(1,1)*G(2,2)+F(1,1)*G(1,3)*G(2,2)+F(1,2)*G(1,1)*G(2,3)-
     &    F(1,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

        Ae(2,1) = (F(2,3)*G(2,2)*G(3,1)-F(2,2)*G(2,3)*G(3,1)-F(2,3)*
     &    G(2,1)*G(3,2)+F(2,1)*G(2,3)*G(3,2)+F(2,2)*G(2,1)*G(3,3)-
     &    F(2,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(2,2) = (F(2,3)*G(1,2)*G(3,1)-F(2,2)*G(1,3)*G(3,1)-F(2,3)*
     &    G(1,1)*G(3,2)+F(2,1)*G(1,3)*G(3,2)+F(2,2)*G(1,1)*G(3,3)-
     &    F(2,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(2,3) = (F(2,3)*G(1,2)*G(2,1)-F(2,2)*G(1,3)*G(2,1)-F(2,3)*
     &    G(1,1)*G(2,2)+F(2,1)*G(1,3)*G(2,2)+F(2,2)*G(1,1)*G(2,3)-
     &    F(2,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

        Ae(3,1) = (F(3,3)*G(2,2)*G(3,1)-F(3,2)*G(2,3)*G(3,1)-F(3,3)*
     &    G(2,1)*G(3,2)+F(3,1)*G(2,3)*G(3,2)+F(3,2)*G(2,1)*G(3,3)-
     &    F(3,1)*G(2,2)*G(3,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))
        Ae(3,2) = (F(3,3)*G(1,2)*G(3,1)-F(3,2)*G(1,3)*G(3,1)-F(3,3)*
     &    G(1,1)*G(3,2)+F(3,1)*G(1,3)*G(3,2)+F(3,2)*G(1,1)*G(3,3)-
     &    F(3,1)*G(1,2)*G(3,3))/(-G(1,3)*G(2,2)*G(3,1)+G(1,2)*G(2,3)
     &    *G(3,1)+G(1,3)*G(2,1)*G(3,2)-G(1,1)*G(2,3)*G(3,2)-G(1,2)
     &    *G(2,1)*G(3,3)+G(1,1)*G(2,2)*G(3,3))
        Ae(3,3) = (F(3,3)*G(1,2)*G(2,1)-F(3,2)*G(1,3)*G(2,1)-F(3,3)*
     &    G(1,1)*G(2,2)+F(3,1)*G(1,3)*G(2,2)+F(3,2)*G(1,1)*G(2,3)-
     &    F(3,1)*G(1,2)*G(2,3))/(G(1,3)*G(2,2)*G(3,1)-G(1,2)*G(2,3)*
     &    G(3,1)-G(1,3)*G(2,1)*G(3,2)+G(1,1)*G(2,3)*G(3,2)+G(1,2)*
     &    G(2,1)*G(3,3)-G(1,1)*G(2,2)*G(3,3))

      END SUBROUTINE CalAe