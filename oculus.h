#ifdef HAVE_CONFIGF
#include "config.f"
#else
#define HAVE_NAG
#endif

module oculus
  
  implicit none

! constants; 25 Mar 15;

  REAL, parameter :: zero       =   0.0
  REAL, parameter :: one        =   1.0
  REAL, parameter :: two        =   2.0
  REAL, parameter :: three      =   3.0
  REAL, parameter :: four       =   4.0
  REAL, parameter :: five       =   5.0
  REAL, parameter :: six        =   6.0
  REAL, parameter :: seven      =   7.0
  REAL, parameter :: eight      =   8.0
  REAL, parameter :: nine       =   9.0
  REAL, parameter :: ten        =  10.0

  REAL, parameter :: pi2        =   6.28318530717958623
  REAL, parameter :: goldenmean =   1.618033988749895 ! golden mean = ( one + sqrt(five) ) / two ;

  REAL, parameter :: small      = 1.0e-12 ! this should really be machine precision; 25 Mar 15;

  REAL, parameter :: half       = one / two
  
! structures; 25 Mar 15;
  
  type magneticaxis
     INTEGER              :: Nfp
     INTEGER              :: Ntor
     REAL                 :: R, Z
     INTEGER              :: maxits, its
     REAL                 :: odetol, tol, error, residue, iota
     REAL                 :: tangent(1:2,1:2), wr(1:2), wi(1:2), vr(1:2,1:2), vi(1:2,1:2)
     REAL   , allocatable :: Ri(:), Zi(:), Rnc(:), Zns(:), Rns(:), Znc(:)
  end type magneticaxis

  type(magneticaxis)      :: axis


  type homoclinictangle 
     INTEGER              :: Nfp, Nhomoclinic
     REAL                 :: R, Z, dU, dS
     INTEGER              :: its
     REAL                 :: odetol, tol, error, residue
     REAL                 :: tangent(1:2,1:2), wr(1:2), wi(1:2), vr(1:2,1:2), vi(1:2,1:2)
  end type homoclinictangle

  type(homoclinictangle)  :: tangle


! type mappeddata
!    INTEGER              :: Nfp, itangent
!    REAL                 :: Ri, Zi, Ro, Zo, odetol, tangent(1:2,1:2)
! end type mappeddata

! type(mappeddata)        :: rzdata


  type transformdata
     INTEGER              :: Nfp, Ppts
     REAL                 :: Ra, Za, R, Z, iota, odetol
     REAL   , allocatable :: data(:,:)
  end type transformdata

  type(transformdata)     :: transform


  type poincaredata
     INTEGER              :: Nfp, Ppts
     REAL                 :: R, Z, odetol
     REAL   , allocatable :: data(:,:)
  end type poincaredata

  type(poincaredata)      :: poincare


! miscellaneous internal varibles; 25 Mar 15;


  INTEGER                 :: itangent
  
  LOGICAL                 :: Lbfieldok

  INTEGER                 :: izeta, fNtor

  INTEGER                 :: iga00aa, iho00aa, ipp00aa, itr00aa


!latex \section{subroutines}

  
contains

  
#ifdef HAVE_NAG
  subroutine bf00aa( zeta, RZ, BRZ ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
    
    integer, parameter  :: Node = 7
#else    
  subroutine bf00aa(Node, zeta, RZ, BRZ ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
    integer, intent(in) :: Node
#endif
    
    REAL  , intent(in)  :: zeta, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
    
    INTEGER             :: ifail
    REAL                :: RpZ(1:3), dBRpZ(1:3,0:3), teta, TM(1:2,1:2)
    
    BRZ(1:Node) = zero ! default intent out; 25 Mar 15;
    
    RpZ(1:3) = (/ RZ(1), zeta, RZ(2) /) ; teta = RZ(7) ! internal shorthand; 25 Mar 15;
    
    call bfield( RpZ(1:3), itangent, dBRpZ(1:3,0:3), ifail ) ! user supplied routine for magnetic field; 25 Mar 15;
    if( ifail.ne.0 ) then ; Lbfieldok = .false. ; return
    endif
    
    BRZ(1:2) = (/ dBRpZ(1,0), dBRpZ(3,0) /) /  dBRpZ(2,0) ! normalize to toroidal field; 25 Mar 15;
    
    if( itangent.eq.0 ) return

    TM(1,1) = ( dBRpZ(1,1) - BRZ(1) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(1,2) = ( dBRpZ(1,3) - BRZ(1) * dBRpZ(2,3) ) / dBRpZ(2,0)
    TM(2,1) = ( dBRpZ(3,1) - BRZ(2) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(2,2) = ( dBRpZ(3,3) - BRZ(2) * dBRpZ(2,3) ) / dBRpZ(2,0)
     
    BRZ(3) = TM(1,1) * RZ(3) + TM(1,2) * RZ(5) ! tangent map obtained by matrix multiplication;
    BRZ(4) = TM(1,1) * RZ(4) + TM(1,2) * RZ(6)
    BRZ(5) = TM(2,1) * RZ(3) + TM(2,2) * RZ(5)
    BRZ(6) = TM(2,1) * RZ(4) + TM(2,2) * RZ(6)
     
    BRZ(7) = cos(teta) * ( TM(2,1) * cos(teta) + TM(2,2) * sin(teta) ) - sin(teta) * ( TM(1,1) * cos(teta) + TM(1,2) * sin(teta) ) 
     
    return
    
  end subroutine bf00aa


  subroutine bf00ab( zeta, RZ, BRZ ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
    
    INTEGER, parameter  :: Node = 5
    
    REAL  , intent(in)  :: zeta, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
    
    INTEGER             :: afail, bfail
    REAL                :: RpZa(1:3), RpZb(1:3), dBRpZa(1:3,0:3), dBRpZb(1:3,0:3), dR, dZ
    
    BRZ(1:Node) = zero ! default intent out; 25 Mar 15;
    
    RpZa(1:3) = (/ RZ(1), zeta, RZ(2) /)
    RpZb(1:3) = (/ RZ(3), zeta, RZ(4) /)
    
    call bfield( RpZa(1:3), itangent, dBRpZa(1:3,0:3), afail ) ! user supplied routine for magnetic field; 25 Mar 15;
    call bfield( RpZb(1:3), itangent, dBRpZb(1:3,0:3), bfail ) ! user supplied routine for magnetic field; 25 Mar 15;

    if( afail.ne.0 .or. bfail.ne.0 ) then ; Lbfieldok = .false. ; return
    endif
    
    BRZ(1:2) = (/ dBRpZa(1,0), dBRpZa(3,0) /) /  dBRpZa(2,0) ! normalize to toroidal field; 25 Mar 15;
    BRZ(3:4) = (/ dBRpZb(1,0), dBRpZb(3,0) /) /  dBRpZb(2,0) ! normalize to toroidal field; 25 Mar 15;
    
    dR = RpZb(1) - RpZa(1)
    dZ = RpZb(3) - RpZa(3)

    BRZ(5) = ( dR * ( BRZ(4)-BRZ(2) ) - dZ * ( BRZ(3)-BRZ(1) ) ) / ( dR*dR + dZ*dZ )

    return
    
  end subroutine bf00ab


!latex \subsection{ga00aa : find the magnetic axis}
!latex \bi
!latex \item[1.] Iterative fieldline tracing methods are used to find the magnetic axis,
!latex           defined as the magnetic fieldline that closes on itself after a toroidal distance of $\Delta \p = 2\pi$/\verb+Nfp+, i.e. ${\bf x}(\Delta\p)={\bf x}(0)$,
!latex           where \verb+Nfp+ is the field periodicity.
!latex \item[ *] The fieldline mapping is defined by integrating along the magnetic field, and is constructed numerically in cylindrical coordinates by integrating the o.d.e.'s
!latex           \be \frac{dR(\p)}{d\p} & = & \frac{B^R(R,\p,Z)}{B^\p(R,\p,Z)} \equiv \dot R(R,\p,Z), \label{eq:BR} \\
!latex               \frac{dZ(\p)}{d\p} & = & \frac{B^Z(R,\p,Z)}{B^\p(R,\p,Z)} \equiv \dot Z(R,\p,Z), \label{eq:BZ}
!latex           \ee
!latex           from an initial, user-supplied starting point, $(R_0,0,Z_0)$.
!latex           The toroidal angle, $\p$, is used as the integration parameter, and so $B^\p$ cannot be zero.
!latex           Upon request, this routine will be modified in order to follow field lines in regions where $B^\p=0$.
!latex \item[ *] A Newton-iterative method is used to find the zero of
!latex           \be {\bf f}\left(\begin{array}{c}R_0\\Z_0\end{array}\right) \equiv \left(\begin{array}{c}R_1-R_0\\Z_1-Z_0\end{array}\right)
!latex           \ee
!latex           where $R_1 \equiv R(\Delta\p)$ and $Z_1 \equiv Z(\Delta\p)$.
!latex \item[ *] Given an initial guess, ${\bf x}\equiv(R_0, Z_0)^T$,
!latex           a better guess for the location of the axis, $(R_0,Z_0)^T+(\delta R,\delta Z)^T$, is given by the linear approximation
!latex           \be {\bf f}\left(\begin{array}{c}R_0+\delta R_0\\Z_0+\delta Z_0\end{array}\right) = {\bf f}\left(\begin{array}{c}R_0\\Z_0\end{array}\right)
!latex               + \underbrace{
!latex               \left(\begin{array}{lcl}\partial_{R_0}R_1 -1&, & \partial_{Z_0}R_1 \\ \partial_{R_0}Z_1&, & \partial_{Z_0}Z_1-1\end{array}\right)}_{\nabla {\bf f}}  \cdot
!latex               \left(\begin{array}{c}\delta R_0\\ \delta Z_0\end{array}\right) + {\cal O}(\delta^2) = 0,
!latex           \ee
!latex           and the correction is given by $\delta{\bf x} = - (\nabla {\bf f})^{-1} \cdot {\bf f}({\bf x})$.
!latex \item[ *] The derivatives, $\partial_{R_0}R_1$, $\partial_{Z_0}R_1$, etc. are determined by fieldline integration,
!latex           \be \frac{d}{d\p} \left( \begin{array}{cc}\partial_{R_0}R(\p), & \partial_{Z_0}R(\p) \\ \partial_{R_0}Z(\p), & \partial_{Z_0}Z(\p) \end{array}\right) = 
!latex                             \left( \begin{array}{cc}\partial_{R  }\dot R, & \partial_{Z  }\dot R\\ \partial_{R  }\dot Z, & \partial_{Z  }\dot Z \end{array}\right) \cdot
!latex                             \left( \begin{array}{cc}\partial_{R_0}R(\p), & \partial_{Z_0}R(\p) \\ \partial_{R_0}Z(\p), & \partial_{Z_0}Z(\p) \end{array}\right),
!latex           \ee
!latex           from an initial starting point being the identity matrix,
!latex           \be               \left( \begin{array}{cc}\partial_{R_0}R( 0), & \partial_{Z_0}R( 0) \\ \partial_{R_0}Z( 0), & \partial_{Z_0}Z( 0) \end{array}\right) = 
!latex                             \left( \begin{array}{cc} 1, & 0 \\ 0,& 1                                                                            \end{array}\right).
!latex           \ee
!l!tex           i.e. $\partial_{R_0}R_1(0)=1$, $\partial_{Z_0}R_1(0)=0$, $\partial_{R_0}Z_1(0)=0$ and $\partial_{Z_0}Z_1(0)=1$.
!latex \item[ *] The above definition of the magnetic axis does not have a unique solution:
!latex           an $\iotabar=1/1$ fieldline also satisfies this definition, as does the $\iotabar=2/1$, $3/1$, etc., as also does the ``X'' point at the separatrix.
!latex           Furthermore, during a sawteeth cycle, the $\iotabar=1/1$ fieldline and the original magnetic axis swap places.
!latex           If there is a continuous family of ``magnetic axes", e.g. there exists an intact $q=1$ surface,
!latex           then $\nabla {\bf f}$ will not be invertible (unless singular value decomposition methods are used).
!latex           Thus, this routine should be used with care.
!latex           Which closed field line that \verb+ga00aa+ locates is determined by the initial guess provided.
!l!tex           Only the ``true'' magnetic axis should be located with \verb+ga00aa+.
!l!tex           Frequently, the coordinate harmonics of the magnetic axis returned by \verb+ga00aa+ will be used as the coordinate axis of toroidal ``chaotic'' coordinates.
!l!tex           Other routines will be provided for locating the $X$ point etc.
!l!tex \item[ *] The definition of the ``true'' magnetic axis may be made precise by examining the Fourier harmonic content of the closed fieldline,
!l!tex           and defining the true magnetic axis may be that closed fieldline with minimal Fourier content.
!latex \item[ *] The returned information includes: 
!latex \bi
!latex \item[i.] the Fourier representation of $R(\p)$ and $Z(\p)$;
!latex \item[ii.] the tangent-mapping near the axis, which allows the rotational-transform on axis to be determined;
!latex \item[iii.] Greene's residue \cite{Greene_79} calculated at the magnetic axis (determines stability).
!latex \ei
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : magneticaxis, ga00aa+ \\ \\
!latex           \verb+type(magneticaxis) :: axis+ \\ \\ in their source that calls \verb+ga00aa+,
!latex           where \verb+axis+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+axis+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+axis%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+axis%Ntor             : integer ;+
!latex \bi
!latex \item[i.] the desired Fourier resolution of the magnetic axis, 
!latex \item[ii.] if it is not required to have a Fourier decomposition of the magnetic axis, or the magnetic field is axisymmetric, choose \verb+Ntor=0+;
!latex \ei
!latex \item[  ] \verb+axis%R               : real    ;+
!latex \bi
!latex \item[i.] guess for the $R$ location of the magnetic axis on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+axis%Z               : real    ;+
!latex \bi
!latex \item[i.] guess for the $Z$ location of the magnetic axis on the $\phi=0$ plane;
!latex \ei
!latex \item[  ] \verb+axis%maxits           : integer ;+
!latex \bi
!latex \item[i.] max. iterations allowed in search;
!latex \item[ii.] e.g. \verb+maxits=16+;
!latex \ei
!latex \item[  ] \verb+axis%tol              : real    ;+
!latex \bi
!latex \item[i.] required accuracy to which the position of the magnetic axis on the $\p=0$ plane is required
!latex \item[ii.] e.g. \verb+tol=1.0e-06+;
!latex \ei
!latex \item[  ] \verb+axis%odetol           : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+ifail                      : integer ;+
!latex \bi
!latex \item[i.] error control flag:
!latex \item[ii.] usually, the user should choose \verb+ifail = 0+, which will result in some information being printed to screen;
!latex \item[iii.] ``quiet'' mode is \verb-ifail = +1-;
!latex \item[iv.] ``noisy'' mode is \verb+ifail = -1+, which will show additional information that may be helpful for debugging.
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call ga00aa( axis, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+axis%R               : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+axis%Z               : real    ;+
!latex \bi
!latex \item[i.] updated;
!latex \ei
!latex \item[  ] \verb+axis%tangent(1:2,1:2) : real    ;+ 
!latex \bi
!latex \item[i.]  the tangent mapping at axis;
!latex \item[ii.] if the eigenvalues of the tangent map are imaginary, e.g. $\lambda\equiv\alpha+\beta i$, 
!latex            then the rotational-transform on axis satisfies $\tan(||\iotabar||)=\beta/\alpha$, where $||\iotabar||\equiv\iotabar \mod 2\pi$.
!latex \item[iii.] if the eigenvalues of the tangent map are real, then the eigenvalues give the direction of the stable and unstable manifolds.
!latex \ei
!latex \item[  ] \verb+axis%wr(1:2) : real    ;+ 
!latex \item[  ] \verb+axis%wi(1:2) : real    ;+ 
!latex \item[  ] \verb+axis%vr(1:2,1:2) : real    ;+ 
!latex \item[  ] \verb+axis%vi(1:2,1:2) : real    ;+ 
!latex \bi
!latex \item[i.]  the eigenvalues and eigenvectors of the tangent mapping at axis;
!latex \ei
!latex \item[  ] \verb+axis%iota         : real    ;+
!latex \bi
!latex \item[i.] rotational-transform on axis
!latex \item[ii.] will only be meaningful if axis is stable, which is indicated by both sign of residue and eigenvalues and eigenvectors;
!latex \ei
!latex \item[  ] \verb+axis%Ri(0:4*Ntor)  : real    ;+
!latex \bi
!latex \item[i.] the magnetic axis, $R(i\Delta\varphi)$, for $i=0,4*$\verb+Ntor+, where $\Delta\varphi=\Delta\p/(4*$\verb+Ntor+$)$;
!latex \item[ii.] \verb+Ri+ is allocated internally; if on input \verb+Ri+ is already allocated it will first be deallocated; similarly for \verb+Zi+
!latex \ei
!latex \item[  ] \verb+axis%Zi(0:4*Ntor)  : real    ;+
!latex \bi
!latex \item[i.] the magnetic axis, $Z(i\Delta\varphi)$, for $i=0,4*$\verb+Ntor+, where $\Delta\varphi=\Delta\p/(4*$\verb+Ntor+$)$;
!latex \ei
!latex \item[  ] \verb+axis%Rnc(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $R(\phi)=\sum_n [R_{n,c} \cos(-n\phi)+R_{n,s} \sin(-n\phi)]$;
!latex \item[ii.] \verb+Rnc+ is allocated internally; if on input \verb+Rnc+ is already allocated it will first be deallocated;
!latex            similarly for \verb+Zns+, \verb+Rns+ and \verb+Znc+.
!latex \ei
!latex \item[  ] \verb+axis%Zns(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $Z(\phi)=\sum_n [Z_{n,c} \cos(-n\phi)+Z_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%Rns(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $R(\phi)=\sum_n [R_{n,c} \cos(-n\phi)+R_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%Znc(0:Ntor)      : real    ;+
!latex \bi
!latex \item[i.] the Fourier harmonics, $Z(\phi)=\sum_n [Z_{n,c} \cos(-n\phi)+Z_{n,s} \sin(-n\phi)]$;
!latex \ei
!latex \item[  ] \verb+axis%error        : real    ;+
!latex \bi
!latex \item[i.] the error, $\sqrt{\Delta R^2 + \Delta Z^2}$, where $\Delta R \equiv R(\Delta\p)-R_0$ and $\Delta Z \equiv Z(\Delta\p)-Z_0$.
!latex \ei
!latex \item[  ] \verb+axis%its          : integer ;+
!latex \bi
!latex \item[i.] the number of iterations required;
!latex \item[ii.] if a good initial guess is given, this number should be small, as Newton methods should converge rapidly; 
!latex            however, if there are multiple magnetic axes (as during a sawtooth event) then the Newton method may encounter problems;
!latex            also, numerical errors in the magnetic field (perhaps $\nabla\cdot{\bf B}$ is not exactly zero) can cause the fieldline integration to be inaccurate, 
!latex            and so it may be difficult to find the solution to the desired accuracy.
!latex \item[iii.] presently, there is no limitation on the maximum number of iterations allowed, but this could easily be changed if required;
!latex \ei
!latex \item[  ] \verb+axis%residue      : real    ;+
!latex \bi
!latex \item[i.] Greene's residue of the magnetic axis; [Greene, 1979] \cite{Greene_79};
!latex \ei
!latex \item[  ] \verb+ifail                      : integer ;+
!latex \bi
!latex \item[i.] on normal execution \verb+ifail=0+;
!latex \item[ii.] more detailed error flags will be provided on request;
!latex \ei
!latex \item[6.] The NAG routine \verb+C05PBF+ is used for the nonlinear root find, and \verb+tol+ is given directly to \verb+C05PBF+.\\
!latex           The NAG routine \verb+D02BJF+ is used for the o.d.e. integration, and \verb+odetol+ is supplied directly to \verb+D02BJF+.
!latex \ei
  

  subroutine ga00aa( laxis, ifail )
#ifndef HAVE_NAG
    use fft2d_mod, only: fft2d
#endif
    
    implicit none
    

    type(magneticaxis)   :: laxis
    INTEGER              :: ifail
    

    integer, parameter   :: NN = 2, Ldfjac = NN, Lworka = NN * ( 3 * NN + 13 ) / 2
    integer, parameter   :: Nmat = 2, Lda = Nmat, Ldvr = Nmat, Ldvi = Nmat, Lworkb = 4 * Nmat
    integer              :: ic05pbf, if02ebf, ic06eaf, iev, astat, iRZ, iflag, ii
    real                 :: tol, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), worka(1:Lworka), tdot
    real                 :: workb(1:Lworkb) ! eigenvalues/vectors of tangent;
    character            :: job
#ifndef HAVE_NAG
    integer              :: lphi,isign
    logical              :: dealflag
    real, allocatable    :: rcoef(:,:), eigv(:,:)
    complex, allocatable :: fcoef(:,:)
    interface
      SUBROUTINE HYBRJ1(FCN,N,X,FVEC,FJAC,LDFJAC,TOL,INFO,WA,LWA)
      INTEGER N,LDFJAC,INFO,LWA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N),WA(LWA)
      EXTERNAL FCN
      END SUBROUTINE HYBRJ1
    end interface
    interface
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
      END SUBROUTINE DGEEV
    end interface
#endif

    
    iga00aa = ifail


! begin : some cursory preliminary checks on input variables; 02 Mar 15;
    
    if( laxis%Nfp.le.0 .or. laxis%Ntor.lt.0 .or. laxis%R.le.zero .or. laxis%tol.le.zero .or. laxis%odetol.le.zero ) then
     ic05pbf = 1
     goto 9000
    endif
    
!  end  : some cursory preliminary checks on input variables; 02 Mar 15;
    
    
! begin : allocate R and Z, etc., which will contain the magnetic axis as a function of toroidal angle; 02 Mar 15;
    
    axis%Nfp = laxis%Nfp ; axis%Ntor = laxis%Ntor ; axis%R = laxis%R 
    axis%Z = laxis%Z ; axis%tol = laxis%tol ; axis%odetol = laxis%odetol 
    axis%maxits = laxis%maxits
    
    fNtor = 4 * max(axis%Ntor,1) ! shorthand; 02 Mar 15;
    
    RALLOCATE(axis%Ri,(0:fNtor))
    RALLOCATE(axis%Zi,(0:fNtor))
    
    if( allocated(laxis%Ri) ) deallocate(laxis%Ri)
    if( allocated(laxis%Zi) ) deallocate(laxis%Zi)
    
    RALLOCATE(laxis%Ri,(0:fNtor))
    RALLOCATE(laxis%Zi,(0:fNtor))
    
    if( allocated(laxis%Rnc) ) deallocate(laxis%Rnc)
    if( allocated(laxis%Zns) ) deallocate(laxis%Zns)
    if( allocated(laxis%Rns) ) deallocate(laxis%Rns)
    if( allocated(laxis%Znc) ) deallocate(laxis%Znc)
    
    RALLOCATE(laxis%Rnc,(0:axis%Ntor))
    RALLOCATE(laxis%Zns,(0:axis%Ntor))
    RALLOCATE(laxis%Rns,(0:axis%Ntor))
    RALLOCATE(laxis%Znc,(0:axis%Ntor))
    
!  end  : allocate R and Z, which will contain the magnetic axis as a function of toroidal angle; 02 Mar 15;
    
    
! begin : call NAG nonlinear root finder; 02 Mar 15;
    
    ic05pbf = 1 ; tol = axis%tol ; axis%its = 0 ; RZ(1:2) = (/ axis%R, axis%Z /)
    
#ifdef HAVE_NAG
    call C05PBF( ga00ab, NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, tol, worka(1:Lworka), Lworka, ic05pbf )
#else 
    call HYBRJ1( ga00ab, NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, tol, ic05pbf, worka(1:Lworka), Lworka )
    if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif
    
!  end  : call NAG nonlinear root finder; 02 Mar 15;
    
    laxis%Ri(0:fNtor) = axis%Ri(0:fNtor)
    laxis%Zi(0:fNtor) = axis%Zi(0:fNtor)
    
    laxis%its = axis%its ; laxis%iota = zero
     
! begin : if root finder was successful, construct relevant output quantities ; 02 Mar 15;
    
    if( ic05pbf.ne.0 ) then
     
    write(0,'("ga00aa : ic05pbf="i2" ; error locating axis; its="i3" ;")') ic05pbf, axis%its
     
     laxis%Rnc(0:laxis%Ntor) = zero
     laxis%Zns(0:laxis%Ntor) = zero
     laxis%Rns(0:laxis%Ntor) = zero
     laxis%Znc(0:laxis%Ntor) = zero
     
     laxis%residue = zero ; laxis%tangent(1,1) = one ; laxis%tangent(1,2) = zero 
     laxis%tangent(2,1) = zero ; laxis%tangent(2,2) = one

     laxis%error = sqrt(sum(FRZ(1:NN)**2))
     
    else
     
     iflag = 2 ; call ga00ab( NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, iflag ) ! make sure axis and tangent map are constructed at solution; was 21 Mar 15;

     laxis%error = sqrt(sum(FRZ(1:NN)**2))
     
     laxis%R = RZ(1) ; laxis%Z = RZ(2) ; laxis%residue = axis%residue
     
     if( laxis%residue.gt.zero ) laxis%iota = axis%iota

     laxis%tangent(1:2,1:2) = axis%tangent(1:2,1:2)
     
     ! Fourier routines
#ifdef HAVE_NAG
    ! c06eaf is NAG routine for Fourier modes
     ic06eaf = 1 ; call C06EAF( axis%Ri(0:fNtor-1), fNtor, ic06eaf ) ! this is destructive; 02 Mar 15;
     ic06eaf = 1 ; call C06EAF( axis%Zi(0:fNtor-1), fNtor, ic06eaf ) ! this is destructive; 02 Mar 15;
#else
!     lphi=log(real(fNtor))/log(2.) ! JRK - round off error here?
!     isign=1
!     dealflag=.false.
!     allocate(rcoef(1,1:fNtor),fcoef(1,1:fNtor/2+1))
!     rcoef(1,:)=axis%Ri
!     CALL fft2d(rcoef,fcoef,1,lphi,isign,dealflag) 
!     !axis%Ri(0:fNtor/2)=real(fcoef(1:fNtor/2+1))
!     ! Indexing is tricky here needs to be reversed
!     !do ii=fNtor-1,fNtor/2,-1
!     !  axis%Ri(ii)=IMAG(fcoef(fNtor+1-ii))
!     !enddo
!     rcoef(1,:)=axis%Zi
!     CALL fft2d(rcoef,fcoef,1,lphi,isign,dealflag) 
!     !TODO - unpack the fcoef array into axis%Zi
!     deallocate(rcoef,fcoef)
#endif

     
     laxis%Rnc(0:axis%Ntor) = axis%Ri(      0:      axis%Ntor   ) * one / sqrt( one * fNtor )
     laxis%Rns(1:axis%Ntor) = axis%Ri(fNtor-1:fNtor-axis%Ntor:-1) * two / sqrt( one * fNtor )
     laxis%Znc(0:axis%Ntor) = axis%Zi(      0:      axis%Ntor   ) * one / sqrt( one * fNtor )
     laxis%Zns(1:axis%Ntor) = axis%Zi(fNtor-1:fNtor-axis%Ntor:-1) * two / sqrt( one * fNtor )
     
    endif ! end of if( c05pbf.eq.0 ) ; 02 Mar 15;
    
!  end  : if root finder was successful, construct relevant output quantities ; 02 Mar 15;

    DEALLOCATE(axis%Ri)
    DEALLOCATE(axis%Zi)
9000 continue

    job = 'V' ; if02ebf = 1 ! construct eigenvalues etc of tangent map at magnetic axis; 02 Mar 15;
#ifdef HAVE_NAG
    call F02EBF( job, Nmat, axis%tangent(1:Lda,1:Nmat), Lda, &
                 laxis%wr(1:Nmat), laxis%wi(1:Nmat), laxis%vr(1:Ldvr,1:Nmat), Ldvr, laxis%vi(1:Ldvi,1:Nmat), Ldvi, &
                 workb(1:Lworkb), Lworkb, if02ebf )
#else
#ifdef HAVE_LAPACK
    allocate(eigv(1:Lda,1:Nmat))
    call DGEEV( 'N', job, Nmat, axis%tangent(1:Lda,1:Nmat), Lda, laxis%wr(1:Nmat), laxis%wi(1:Nmat), &
                eigv, Ldvi, eigv, Ldvi, &
                workb(1:Lworkb), Lworkb, if02ebf )
    laxis%vr(1,:)=eigv(:,1)
    laxis%vr(2,:)=eigv(:,2)
    laxis%vi=0. ! Assume real eigenvectors
    deallocate(eigv)
#else
    write(0,'(" oculus : ", 10x ," : requires LAPACK")')
    stop 'oculus requires lapack'
#endif
#endif
     
     if( laxis%residue.gt.zero ) then
      tdot = atan2(laxis%wi(1),laxis%wr(1)) / (pi2/laxis%Nfp)
      laxis%iota = tdot
     !if    ( laxis%iota - tdot.lt.-one/laxis%Nfp ) then ; laxis%iota = tdot - one / laxis%Nfp
     !elseif( laxis%iota - tdot.gt. one/laxis%Nfp ) then ; laxis%iota = tdot + one / laxis%Nfp
     !else                                               ; laxis%iota = tdot
     !endif
     endif

9999 continue

    ifail = ic05pbf

    return
    
!1000 format("ga00aa : ic05pbf="i2" ; ( R, Z ) = ("es13.05" ,"es13.05" ) ; error=",es10.3" ; residue="es13.5" ; its="i3" ;")
!1010 format("ga00aa :         "2x" ; tangent(1:2,1:2) = "4es13.5," ;")
!1020 format("ga00aa : if02ebf="i2" ; ":"evalue="f8.3" +"f8.3" i ; evector=("f8.3" +"f8.3" i,"f8.3" +"f8.3" i ) ;")
    
  end subroutine ga00aa

    
  subroutine ga00ab( NN, oRZ, FRZ, dFRZ, Ldfjac, iflag )
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif
    
    implicit none
    
    INTEGER, intent(in)    :: NN
    REAL                   :: oRZ(1:NN), FRZ(1:NN)
    
    INTEGER, intent(in)    :: Ldfjac
    REAL                   :: dFRZ(1:Ldfjac,1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    integer, parameter     :: Node = 7, Lworka = 20 * Node
    
    integer                :: id02bjf
    real                   :: dRZ(1:Node), phistart, phiend, tol
#ifdef HAVE_NAG
    real                   :: worka(1:Lworka)
#else
    integer, parameter     :: liw=20, lrw=20+16*Node
    integer                :: iwork(liw)
    real                   :: rwork(lrw)
    integer                :: iopt,istate,itask,itol,mf,ii
    real                   :: atol,rtol,dphi,phic,phie
#endif
    
    CHARACTER              :: relabs
    
#ifdef HAVE_NAG
    external               :: D02BJW
#endif
    
    tol = axis%odetol ; relabs = 'D'
    
    phistart = zero ; phiend = phistart + ( pi2 / axis%Nfp ) ! integration endpoints ; 05 Mar 14;
    
    dRZ(1:Node) = (/ oRZ(1), oRZ(2),  one, zero, zero, one,  zero /) ! initial guess, intialize tangent map integration; intialize transform; 31 Jul 13;
    
    izeta = 0 ! this counter is incremented in ga00ac; 31 Jul 13;
    
    select case( iflag )
    case( 1 ) ; itangent = 0 ! derivatives (i.e. tangent map) is not required; 21 Mar 15; this is passed through to user-supplied bfield;
    case( 2 ) ; itangent = 1 ! derivatives (i.e. tangent map) is     required; 21 Mar 15; this is passed through to user-supplied bfield;
    end select
    
    id02bjf = 1
    
    Lbfieldok = .true.

#ifdef HAVE_NAG
    ! d02bjf is the ODE integrator
    call D02BJF( phistart, phiend, Node, dRZ(1:Node), bf00aa, tol, relabs, ga00ac, D02BJW, worka(1:Lworka), id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=tol
!    atol :  absolute tolerance
     atol=tol
!    initializations for loop
     dphi=(phiend-phistart)/fNtor
     phic=phistart
     phie=phic
     call ga00ac(phic,dRZ)
     do ii=1,fNtor
       phie=phie+dphi
       call lsode(bf00aa,(/Node/),dRZ,phic,phie, &
                  itol,(/rtol/),(/atol/),itask,istate,iopt, &
                  rwork,lrw,iwork,liw,mf=mf)
       call ga00ac(phie,dRZ)
     enddo
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
    
    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
  
    if( id02bjf.ne.0 ) then
     
    write(0,'("ga00ab : id02bjf="i2" ; error integrating along field;")') id02bjf
     
     axis%tangent(1,1) =  one   ; axis%tangent(1,2) = zero ! supply dummy values; 21 Mar 15;
     axis%tangent(2,1) = zero   ; axis%tangent(2,2) =  one
     
     axis%residue = zero ! supply dummy values; 21 Mar 15;
     
     axis%iota = zero
     
     iflag = -1 ! tell NAG that an error has occured;  5 Jun 13;
     
    else
     
     select case( iflag ) 
     case( 1 ) 
       FRZ(1  ) = dRZ(1) - oRZ(1)                            ! must return function;  5 Jun 13;
       FRZ(2  ) = dRZ(2) - oRZ(2)
     case( 2 )
       dFRZ(1,1) = dRZ(3) - one    ; dFRZ(1,2) = dRZ(4)       ! must return Jacobian;  5 Jun 13;
       dFRZ(2,1) = dRZ(5)          ; dFRZ(2,2) = dRZ(6) - one
     end select
     
     axis%tangent(1,1) = dRZ(3) ; axis%tangent(1,2) = dRZ(4)
     axis%tangent(2,1) = dRZ(5) ; axis%tangent(2,2) = dRZ(6) 

     axis%residue = ( two - ( dRZ(3) + dRZ(6) ) ) / four
    !determinant = dRZ(3)*dRZ(6) - dRZ(4)*dRZ(5) ! this should be equal to unity at fixed points; 21 Mar 15;
     
     axis%iota = dRZ(7) / ( pi2 / axis%Nfp )

     if( iflag.eq.1 ) then
     if( iga00aa.lt.0 ) write(0,'("ga00ab : id02bjf="i2" ; its="i3" ; xx="2es13.5" ; FF="2es13.5" ;")') &
                              id02bjf, axis%its, oRZ(1:2), FRZ(1:2) ! screen output; 25 Mar 15;
      axis%its = axis%its + 1
     endif

     if( axis%its.ge.axis%maxits ) iflag = -1

    endif
    
    return
    
  end subroutine ga00ab

  
  subroutine ga00ac( zeta, RZ )
    
    implicit none
    
    INTEGER, parameter :: Node = 7
    
    REAL               :: zeta, RZ(1:Node)
    
    axis%Ri(izeta) = RZ(1) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
    axis%Zi(izeta) = RZ(2) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
    
    izeta = izeta + 1
    
    zeta = izeta * ( pi2/axis%Nfp ) / fNtor
    
    return
    
  end subroutine ga00ac


!latex \subsection{ho00aa : find the unstable manifold}
!latex \bi
!latex \item[1.] Documentation under construction; contact \verb+shudson@pppl.gov+ if you wish to use this routine..
!latex \ei

  subroutine ho00aa( ltangle, ifail )
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif

    implicit none

    type(homoclinictangle) :: ltangle
    INTEGER                :: ifail
    

    integer, parameter     :: NN = 2, Ldfjac = NN, Lworka = NN * ( 3 * NN + 13 ) / 2
    integer, parameter     :: Nmat = 2, Lda = Nmat, Ldvr = Nmat, Ldvi = Nmat, Lworkb = 4 * Nmat
    integer                :: ic05pbf, if02ebf, iev, astat, iRZ, iflag, ipip, ic05nbf, id02bjf, ii
    real                   :: RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), worka(1:Lworka) ! nonlinear root finder; 02 Mar 15;
    real                   :: tempmatrix(1:Lda,1:Nmat), workb(1:Lworkb) ! eigenvalues/vectors of tangent;
    REAL                   :: XX(1:2), FF(1:2)
    CHARACTER              :: job
#ifndef HAVE_NAG
    real, allocatable      :: eigv(:,:)
#endif
#ifdef UNDERCONSTRUCTION
    INTEGER, parameter     :: Node = 7, Lworkc = 20 * Node
    REAL                   :: lRZ(1:Node,1:2), phistart, phiend
    CHARACTER              :: relabs
#ifdef HAVE_NAG
    real                   :: workc(1:Lworkc)
    external               :: D02BJX, D02BJW
#else
    integer, parameter     :: liw=20, lrw=20+16*Node
    integer                :: iwork(liw)
    real                   :: rwork(lrw)
    integer                :: iopt,istate,itask,itol,mf
    real                   :: atol,rtol

    interface
      SUBROUTINE HYBRJ1(FCN,N,X,FVEC,FJAC,LDFJAC,TOL,INFO,WA,LWA)
      INTEGER N,LDFJAC,INFO,LWA
      DOUBLE PRECISION TOL
      DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N),WA(LWA)
      EXTERNAL FCN
      END SUBROUTINE HYBRJ1
    end interface
    interface
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, &
                        LDVR, WORK, LWORK, INFO )
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ), &
                         WI( * ), WORK( * ), WR( * )
      END SUBROUTINE DGEEV
    end interface
#endif
#endif
    
    if( ltangle%Nfp.le.0 .or. ltangle%R.le.zero .or. ltangle%odetol.le.zero .or. ltangle%tol.le.zero ) then
     ic05pbf = -1
     goto 9999
    endif
    
    
    iho00aa = ifail
    
    tangle%Nfp = ltangle%Nfp ; tangle%odetol = ltangle%odetol

    
    ic05pbf = 1 ; tangle%its = 0 ; RZ(1:2) = (/ ltangle%R, ltangle%Z /)
    
#ifdef HAVE_NAG
    call C05PBF( ho00ab, NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, ltangle%tol, worka(1:Lworka), Lworka, ic05pbf )
#else 
    call HYBRJ1( ho00ab, NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, ltangle%tol, ic05pbf, worka(1:Lworka), Lworka )
    if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif
    
    ltangle%its = tangle%its ; ltangle%error = sqrt(sum(FRZ(1:NN)**2))
    
    
    if( ic05pbf.ne.0 ) then
     
     write(0,'("ho00aa : ic05pbf="i2" ; error locating X point ; its="i3" ;")') ic05pbf, tangle%its
     
     ltangle%residue = zero ; ltangle%tangent(1,1) = one ; ltangle%tangent(1,2) = zero 
     ltangle%tangent(2,1) = zero ; ltangle%tangent(2,2) = one
     
    else
     
     iflag = 2 ; call ho00ab( NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, iflag )
     
     ltangle%R = RZ(1) ; ltangle%Z = RZ(2) 
     ltangle%tangent(1:Lda,1:Nmat) = tangle%tangent(1:Lda,1:Nmat) 
     ltangle%residue = tangle%residue
     
    endif
    
    
    job = 'V' ; if02ebf = 1 ! construct eigenvalues etc of tangent map; 02 Mar 15;
#ifdef HAVE_NAG
    call F02EBF( job, Nmat, tangle%tangent(1:Lda,1:Nmat), Lda, &
                 ltangle%wr(1:Nmat), ltangle%wi(1:Nmat), ltangle%vr(1:Ldvr,1:Nmat), Ldvr, ltangle%vi(1:Ldvi,1:Nmat), Ldvi, &
                 workb(1:Lworkb), Lworkb, if02ebf )
#else
#ifdef HAVE_LAPACK
     allocate(eigv(1:Lda,1:Nmat))
     call DGEEV( 'N', job, Nmat, tangle%tangent(1:Lda,1:Nmat), Lda, ltangle%wr(1:Nmat), ltangle%wi(1:Nmat), &
                 eigv, Ldvr, eigv, Ldvi, workb(1:Lworkb), Lworkb, if02ebf )
     ltangle%vr(1,:)=eigv(:,1)
     ltangle%vr(2,:)=eigv(:,2)
     ltangle%vi=0. ! Assume real eigenvectors
     deallocate(eigv)
#else
     write(0,'(" oculus : ", 10x ," : requires LAPACK")')
     stop 'oculus requires lapack'
#endif
#endif
    
    
    if( ltangle%Nhomoclinic.le.0 ) goto 9999


#ifdef UNDERCONSTRUCTION
    
    FATALERROR(ho00aa, tangle%dU.lt.small, initial guess needs to be provided)
    
    FATALERROR(ho00aa, abs(tangle%wi(1))+abs(tangle%wi(1)).gt.small, X point is not unstable)
    
    FATALERROR(ho00aa, tangle%wr(1).lt.one, unstable direction needs to be swapped)
    
    write(0,1030) tangle%Nhomoclinic, tangle%R, tangle%Z, tangle%dU, tangle%dS, tangle%tol, tangle%odetol
    

    
    lRZ(1:2,2) = (/ tangle%R, tangle%Z /)
    
    lRZ(1:2,1) = (/ tangle%R, tangle%Z /) + tangle%dU * tangle%vr(1:2,1)
    lRZ(3:7,1) = (/ one, zero, zero, one,  zero /)
    
    phistart = zero
    
    tangle%Nhomoclinic = -1
    
    do ; tangle%Nhomoclinic = tangle%Nhomoclinic + 1
     phiend = phistart + ( pi2 / tangle%Nfp)
     id02bjf = 0 ; relabs = 'D'; Lbfieldok = .true.
#ifdef HAVE_NAG
     call D02BJF( phistart, phiend, Node, lRZ(1:Node,1), bf00aa, tol, relabs, D02BJX, D02BJW, workc(1:Lworkc), id02bjf )
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=tol
!    atol :  absolute tolerance
     atol=tol
!    do integration
     call lsode(bf00aa,(/Node/),lRZ,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,mf=mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif

     write(0,'("ho00aa : Nhomoclinic="i3" ; RZ="2f15.10" ; Lbfieldok="L2" ;")') tangle%Nhomoclinic, lRZ(1:2,1), Lbfieldok
     if( (lRZ(1,1)-tangle%R)**2 + (lRZ(2,1)-tangle%Z)**2 .lt. (lRZ(1,2)-tangle%R)**2 + (lRZ(2,2)-tangle%Z)**2 ) exit
     lRZ(1:2,2) = lRZ(1:2,1)
    enddo
    
    write(0,'("ho00aa : Nhomoclinic="i3" ;")') tangle%Nhomoclinic
    
    
    do ipip = 1, 2 ! loop over primary intersection points; 19 Nov 14;
     
     select case( ipip )
     case( 1 ) ; XX(1:NN) = (/ tangle%dU, tangle%dS /)
     case( 2 ) ; XX(1:NN) = (/ tangle%dU, tangle%dS /) * sqrt( tangle%wr(1:2) )
     end select
     
     ic05nbf = 1 ; FF(1:2) = one
     
#ifdef HAVE_NAG    
     call C05NBF( ho00ac, NN, XX(1:NN), FF(1:NN), tangle%tol, worka(1:Lworka), Lworka, ic05nbf ) ! root finder; should use derivatives; 25 Mar 15;
#else 
     call HYBRJ1( ho00ac, NN, XX(1:NN), FF(1:NN), tangle%tol, ic05nbf, worka(1:Lworka), Lworka )
     if (ic05pbf==1) ic05pbf=0 ! because the MINPACK error code is stupid
#endif
     
     select case( ic05nbf )
     case( 0 ) ; write(0,'("ho00aa : " 10x " : ( dU, dS )="2es23.15" ; ic05nbf="i3" ; success          ; FF="2es13.5" ;")') XX(1:NN), ic05nbf, FF(1:NN)
     case( 1 ) ; write(0,'("ho00aa : " 10x " : ( dU, dS )="2es23.15" ; ic05nbf="i3" ; input error      ; FF="2es13.5" ;")') XX(1:NN), ic05nbf, FF(1:NN)
     case( 2 ) ; write(0,'("ho00aa : " 10x " : ( dU, dS )="2es23.15" ; ic05nbf="i3" ; consider restart ; FF="2es13.5" ;")') XX(1:NN), ic05nbf, FF(1:NN)
     case( 3 ) ; write(0,'("ho00aa : " 10x " : ( dU, dS )="2es23.15" ; ic05nbf="i3" ; xtol too small   ; FF="2es13.5" ;")') XX(1:NN), ic05nbf, FF(1:NN)
     case( 4 ) ; write(0,'("ho00aa : " 10x " : ( dU, dS )="2es23.15" ; ic05nbf="i3" ; bad progress     ; FF="2es13.5" ;")') XX(1:NN), ic05nbf, FF(1:NN)
     case default
     end select
     
     if( ipip.eq.1 .and. sqrt(sum(FF(1:NN)**2)).lt.sqrt(tangle%tol) ) then ; tangle%dU = XX(1) ; tangle%dS = XX(2) ! update input guesses; 25 Mar 15;
     endif
     
    enddo ! end of do ipip; 19 Nov 14;
    
#endif
    
    
9999 continue
    
    ifail = ic05pbf
    
    return   
    
!1000 format("ho00aa : ic05pbf="i2" ; ( R, Z ) = ("es13.05" ,"es13.05" ) ; error=",es10.3" ; residue="es13.5" ; its="i3" ;")
!1010 format("ho00aa :         "2x" ; tangent(1:2,1:2) = "4es13.5," ;")
!1020 format("ho00aa : if02ebf="i2" ; ":"evalue="f8.3" +"f8.3" i ; evector=("f8.3" +"f8.3" i,"f8.3" +"f8.3" i ) ;")
!1030 format("ho00aa : " 10x " : Nhomoclinic="i4" ; R="f9.5", Z="f9.5" ; dU="es12.5", dS="es12.5" ; xtol="es12.5" ; odetol="es12.5" ;")
    
  end subroutine ho00aa
  

  subroutine ho00ab( NN, oRZ, FRZ, dFRZ, Ldfjac, iflag )
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif
    
    implicit none
    
    INTEGER, intent(in)    :: NN
    REAL                   :: oRZ(1:NN), FRZ(1:NN)
    
    INTEGER, intent(in)    :: Ldfjac
    REAL                   :: dFRZ(1:Ldfjac,1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    INTEGER, parameter     :: Node = 7, Lworka = 20 * Node
    
    INTEGER                :: id02bjf
    REAL                   :: dRZ(1:Node), phistart, phiend, tol
    CHARACTER              :: relabs
#ifdef HAVE_NAG
    real                   :: worka(1:Lworka)
    external               :: D02BJX, D02BJW
#else
    integer, parameter     :: liw=20, lrw=20+16*Node
    integer                :: iwork(liw)
    real                   :: rwork(lrw)
    integer                :: iopt,istate,itask,itol,mf
    real                   :: atol,rtol
#endif
    
    tol = tangle%odetol ; relabs = 'D'
    
    phistart = zero ; phiend = phistart + ( pi2 / tangle%Nfp ) ! integration endpoints ; 05 Mar 14;
    
    dRZ(1:Node) = (/ oRZ(1), oRZ(2),  one, zero, zero, one,  zero /) ! initial guess, intialize tangent map integration; 31 Jul 13;
    
    select case( iflag )
    case( 1 ) ; itangent = 0 ! derivatives (i.e. tangent map) is not required; 21 Mar 15; this is passed through to user-supplied bfield;
    case( 2 ) ; itangent = 1 ! derivatives (i.e. tangent map) is     required; 21 Mar 15; this is passed through to user-supplied bfield;
    end select
    
    id02bjf = 1
    
    Lbfieldok = .true.

#ifdef HAVE_NAG
    call D02BJF( phistart, phiend, Node, dRZ(1:Node), bf00aa, tol, relabs, D02BJX, D02BJW, worka(1:Lworka), id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=tol
!    atol :  absolute tolerance
     atol=tol
!    do integration
     call lsode(bf00aa,(/Node/),dRZ,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,mf=mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
    
    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
  
    if( id02bjf.ne.0 ) then

    !if( iho00aa.lt.0 ) write(0,'("ho00ab : id02bjf="i2" ; its="i3" ; x="2es13.5" ;":" F="2es13.5" ;")') id02bjf, tangle%its, oRZ(1:2)
     
     write(0,'("ho00ab : id02bjf="i2" ; error integrating along field ;")') id02bjf
     
     tangle%tangent(1,1) =  one ; tangle%tangent(1,2) = zero ! supply dummy values; 21 Mar 15;
     tangle%tangent(2,1) = zero ; tangle%tangent(2,2) =  one

     tangle%residue = zero ! supply dummy values; 21 Mar 15;
     
     iflag = -1 ! tell NAG that an error has occured;  5 Jun 13;
     
    else
     
     select case( iflag ) 
     case( 1 ) 
       FRZ(1  ) = dRZ(1) - oRZ(1)                            ! must return function;  5 Jun 13;
       FRZ(2  ) = dRZ(2) - oRZ(2)
     case( 2 ) 
       dFRZ(1,1) = dRZ(3) - one    ; dFRZ(1,2) = dRZ(4)       ! must return Jacobian;  5 Jun 13;
       dFRZ(2,1) = dRZ(5)          ; dFRZ(2,2) = dRZ(6) - one
     end select
     
     tangle%tangent(1,1) = dRZ(3) ; tangle%tangent(1,2) = dRZ(4)
     tangle%tangent(2,1) = dRZ(5) ; tangle%tangent(2,2) = dRZ(6) 

     tangle%residue = ( two - ( dRZ(3) + dRZ(6) ) ) / four

    !tangle%iota = zero

    !determinant = dRZ(3)*dRZ(6) - dRZ(4)*dRZ(5) ! this should be equal to unity at fixed points; 21 Mar 15;
     
     if( iflag.eq.1 ) then
     !if( iho00aa.lt.0 ) write(0,'("ho00ab : id02bjf="i2" ; its="i3" ; xx="2es13.5" ; FF="2es13.5" ;")') id02bjf, tangle%its, oRZ(1:2), FRZ(1:2)
      tangle%its = tangle%its + 1
     endif
     
    endif
    
    return
    
  end subroutine ho00ab
  
  
  subroutine ho00ac( NN, XX, FF, iflag ) ! should use NAG routine that returns derivatives; 25 Mar 15;
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif
    
    implicit none

    INTEGER            :: NN, iflag
    REAL               :: XX(1:NN), FF(1:NN)
    
    INTEGER, parameter :: Node = 7, Lwork= 20 * Node
    INTEGER            :: ius, ihomo, id02bjf
    REAL               :: lRZ(1:Node,1:2), tol, phistart, phiend
    CHARACTER          :: relabs
#ifdef HAVE_NAG
    real               :: work(1:Lwork)
    external           :: D02BJX, D02BJW
#else
    integer, parameter :: liw=20, lrw=20+16*Node
    integer            :: iwork(liw)
    real               :: rwork(lrw)
    integer            :: iopt,istate,itask,itol,mf
    real               :: atol,rtol
#endif
    
    FF(1:NN) = zero ; tol = tangle%odetol / ten ; relabs = 'D'
    
    do ius = 1, 2

     lRZ(1:2,ius) = (/ tangle%R, tangle%Z /) + XX(ius) * tangle%vr(1:2,ius) ! displace along unstable direction; ! 19 Nov 14;
     
     lRZ(3:7,ius) = (/ one, zero, zero, one,  zero /)
     
     ihomo = 0 ; phistart = zero
     
     FATALERROR(ho00ac, tangle%Nhomoclinic.eq.0, under construction)
     
     do ihomo = 1, tangle%Nhomoclinic
      
      if( ius.eq.1 ) phiend = phistart + ( pi2 / tangle%Nfp ) ! forwards or backwards; 19 Nov 14;
      if( ius.eq.2 ) phiend = phistart - ( pi2 / tangle%Nfp ) ! forwards or backwards; 19 Nov 14;
      
      Lbfieldok = .true.

      id02bjf = 1 
#ifdef HAVE_NAG
      call D02BJF(phistart,phiend,Node,lRZ(1:Node,ius),bf00aa,tol,relabs,D02BJX,D02BJW,work(1:Lwork),id02bjf)
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=tol
!    atol :  absolute tolerance
     atol=tol
!    do integration
     call lsode(bf00aa,(/Node/),lRZ,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,mf=mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
      
      if( .not.Lbfieldok ) id02bjf = -2

      select case( id02bjf )
      case( -2 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; Lbfieldok fail ;              ")') id02bjf, lRZ(1:2,ius)
      case( -1 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; user termination ;            ")') id02bjf, lRZ(1:2,ius)
      case(  0 ) ! write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; success ;                     ")') id02bjf, lRZ(1:2,ius)
      case(  1 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; input error ;                 ")') id02bjf, lRZ(1:2,ius)
      case(  2 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; no further progress possible ;")') id02bjf, lRZ(1:2,ius)
      case(  3 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; tol too small ;               ")') id02bjf, lRZ(1:2,ius)
      case(  4 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; xsol not reset ;              ")') id02bjf, lRZ(1:2,ius)
      case(  5 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; xsol not reset ;              ")') id02bjf, lRZ(1:2,ius)
      case(  6 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; function did not change sign ;")') id02bjf, lRZ(1:2,ius)
      case(  7 ) ; write(0,'("ho00ac : id02bjf="i2" : (R,Z)=("2f9.5" ) ; serious error ;               ")') id02bjf, lRZ(1:2,ius)
      case default
      end select
      
     enddo ! end of do ihomo; 19 Nov 14;
     
    enddo ! end of do ius; 19 Nov 14;
    
    FF(1:NN) = lRZ(1:NN,1) - lRZ(1:NN,2)
    
   !write(0,1000) id02bjf, XX(1:2), lRZ(1:2,-1), lRZ(1:2,1), FF(1:2)
 
   !pause
    
    return
    
1000 format("ho00ac : id02bjf="i2" ; xx = ("es22.15","es22.15" ) ; (R,Z)=(" &
            es23.15" ,"es23.15" ) = ("es23.15" ,"es23.15" ) ; ":"ff="2es13.5" ;")
    
  end subroutine ho00ac
  
!latex \subsection{pm00aa : find all fixed points}
!latex \bi
!latex \item[1.] Documentation under construction; contact \verb+shudson@pppl.gov+ if you wish to use this routine..
!latex \ei
    
!  subroutine pm00aa( lrzdata, ifail )
!    
!    implicit none
!    
!    type(mappeddata)       :: lrzdata
!    INTEGER                :: ifail
!    
!    INTEGER, parameter     :: Node = 7, Lworka = 20 * Node
!    
!    INTEGER                :: id02bjf
!    REAL                   :: dRZ(1:Node), phistart, phiend, worka(1:Lworka)
!    CHARACTER              :: relabs
!    
!    external               :: D02BJX, D02BJW
!    
!    relabs = 'D' ; phistart = zero ; phiend = phistart + ( pi2 / lrzdata%Nfp ) ! integration endpoints ; 05 Mar 14;
!    
!    dRZ(1:Node) = (/ lrzdata%Ri, lrzdata%Zi,  one, zero, zero, one,  zero /) ! initial guess, intialize tangent map integration; intialize transform; 31 Jul 13;
!    
!    itangent = lrzdata%itangent
!    
!    id02bjf = 1
!    
!    Lbfieldok = .true.
!    
!    call D02BJF( phistart, phiend, Node, dRZ(1:Node), bf00aa, lrzdata%odetol, relabs, D02BJX, D02BJW, worka(1:Lworka), id02bjf ) ! NAG ode integration;
!    
!    if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
!    
!    if( id02bjf.ne.0 ) then
!     
!     write(0,'("pm00aa : id02bjf="i2" ; error integrating along field;")') id02bjf
!     
!     lrzdata%R = lrzdata%Ri
!     lrzdata%Z = lrzdata%Zi
!     
!     lrzdata%tangent(1,1) =  one   ; lrzdata%tangent(1,2) = zero ! supply dummy values; 21 Mar 15;
!     lrzdata%tangent(2,1) = zero   ; lrzdata%tangent(2,2) =  one
!     
!    else
!     
!     lrzdata%R = dRZ(1)
!     lrzdata%Z = dRZ(2)
!     
!     lrzdata%tangent(1,1) = dRZ(3) ; lrzdata%tangent(1,2) = dRZ(4)
!     lrzdata%tangent(2,1) = dRZ(5) ; lrzdata%tangent(2,2) = dRZ(6) 
!     
!    endif
!    
!    ifail = id02bjf
!    
!    return
!    
!  end subroutine pm00aa



  
!latex \subsection{tr00aa : measure rotational-transform}
!latex \bi
!latex \item[1.] Fieldline tracing methods are used to determine the relative rotational-transform of one fieldline about a given ``reference'' fieldline,
!latex           will usually be a magnetic axis.
!latex \item[ *] The equations governing the fieldlines are the same as that given in \Eqn{BR} and \Eqn{BZ}.
!latex \item[ *] The user must supply two starting points, $(R_a,Z_a)$ and $(R,Z)$.
!latex           The location of the magnetic axis, $(R_a,Z_a)$, can be obtained from a previous call to \verb+ga00aa+;
!latex           however, the values of $(R_a,Z_a)$ and $(R,Z)$ are completely arbitrary.
!latex           What is really measured by this routine is average ``linking'' of one fieldline about another.
!latex \item[ *] A poloidal angle, $\t(\phi)$, is introduced as 
!latex           \be \tan\t(\phi)=\frac{\delta Z(\phi)}{\delta R(\phi)}, \ee
!latex           where $\delta R(\phi) = R(\phi) - R_a(\phi)$ and $\delta Z(\phi) = Z(\phi) - Z_a(\phi)$.
!latex \item[ *] This angle varies with $\phi$ according to 
!latex           \be \frac{d\t}{d\phi} = \frac{\delta R\,(Z^\prime-Z_a^\prime) - \delta Z\,(R^\prime-R_a^\prime)}{\delta R^2 + \delta Z^2}, \label{eq:diotadphi}
!latex           \ee
!latex           where $\prime$ denotes total derivative with respect to $\phi$.
!latex \item[ *] The o.d.e. integration defined in \Eqn{diotadphi} may be initialized with $\t(0)=0$,
!latex           and after a sufficiently large distance, $\Delta\phi$,
!latex           this angle satisfies $\Delta \t \approx \iotabar \Delta\phi$, where $\iotabar$ is the rotational-transform.
!latex \item[ *] (A more accurate calculation of $\iotabar$ is enabled by fitting a straight line to $\t(\phi)$, rather than just subtracting the endpoints . . .)
!latex \item[ *] Formally, the rotational-transform is defined as the limit
!latex           \be \iotabar \equiv \lim_{\Delta \phi \rightarrow \infty} \frac{\Delta \t}{\Delta \phi}.
!latex           \ee
!latex           This limit only converges on regular fieldlines. 
!latex           For irregular, or chaotic, fieldlines, this limit does not converge and the rotational-transform is not defined.
!latex \item[ *] Note that \Eqn{diotadphi} requires knowledge of $R_a(\phi)$ and $Z_a(\phi)$, and $R(\phi)$ and $Z(\phi)$, and these are obtained by integrating
!latex           \Eqn{BR} and \Eqn{BZ}. So, in total there are $5$ coupled o.d.e.s.
!latex \item[2.] \underline{\bf The user must include} \\ \\ 
!latex           \verb+use oculus, only : transformdata, tr00aa+ \\ \\
!latex           \verb+type(transformdata)    :: transform+ \\ \\ in their source that calls \verb+tr00aa+,
!latex           where \verb+transform+ is a derived type (i.e. structure) that contains both the required input and output information.
!latex           The variable name, \verb+transform+, is arbitrary.
!latex \item[3.] \underline{\bf Required inputs}
!latex \item[  ] \verb+transform%Nfp              : integer ;+
!latex \bi
!latex \item[i.] the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+;
!latex \ei
!latex \item[  ] \verb+transform%Pts              : integer ;+
!latex \bi
!latex \item[i.] the number of toroidal transits that the fieldlines will be followed, e.g. \verb+Ppts=100+;
!latex \ei
!latex \item[  ] \verb+transform%odetol           : real    ;+
!latex \bi
!latex \item[i.] o.d.e. integration tolerance; 
!latex \item[ii.] e.g. \verb+odetol=1.0e-08+;
!latex \ei
!latex \item[  ] \verb+transform%Ra               : real    ;+
!latex \item[  ] \verb+transform%Za               : real    ;+
!latex \item[  ] \verb+transform%R                : real    ;+
!latex \item[  ] \verb+transform%Z                : real    ;+
!latex \bi
!latex \item[i.] starting points for fieldline integrations;
!latex \item[ii.] usually, $(R_a,Z_a)$ will be the location of the magnetic axis on the $\phi=0$ plane, and $(R,Z)$ is arbitrary.
!latex \ei
!latex \item[4.] \underline{\bf Execution}
!latex \item[  ] \verb+call tr00aa( transform, ifail )+
!latex \item[5.] \underline{\bf Outputs}
!latex \item[  ] \verb+transform%iota        : real    ;+
!latex \bi
!latex \item[i.] the ``rotational-transform'' of the fieldline starting at $(R,Z)$ relative to the fieldline starting at $(R_a,Z_a)$.
!latex \ei
!latex \item[  ] \verb+transform%data(1:2,0:Ppts)  : real ;+
!latex \bi
!latex \item [i.] the \Poincare plot information;
!latex \ei
!latex \ei
  
  subroutine tr00aa( ltransform, ifail )
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif
    
    implicit none
    
    type(transformdata) :: ltransform
    INTEGER             :: ifail
    
    INTEGER, parameter  :: Node = 5, Lworka = 20 * Node
    INTEGER             :: id02bjf, ipoint, astat
    REAL                :: RZRZt(1:Node), phistart, phiend
    CHARACTER           :: relabs
#ifdef HAVE_NAG
    real                :: worka(1:Lworka)
    external            :: D02BJX, D02BJW
#else
    integer, parameter  :: liw=20, lrw=20+16*Node
    integer             :: iwork(liw)
    real                :: rwork(lrw)
    integer             :: iopt,istate,itask,itol,mf
    real                :: atol,rtol
#endif
    
    itr00aa = ifail
    
    if( ltransform%Ppts.le.0 ) then
     ifail = -1
     write(0,'("tr00aa : ifail  ="i2" : Ppts="i6" ;")') ifail, ltransform%Ppts
     return
    endif
    
    if( ltransform%Nfp.le.0 ) then
     ifail = -1
     write(0,'("tr00aa : ifail  ="i2" : Nfp="i3" ;")') ifail, ltransform%Nfp
     return
    endif
    
    if( ltransform%odetol.le.0 ) then
     ifail = -1
     write(0,'("tr00aa : ifail  ="i2" : odetol="es13.5" ;")') ifail, ltransform%odetol
     return
    endif
    
    relabs = 'D' ; itangent = 0
    
    if( allocated(ltransform%data) ) deallocate(ltransform%data)
    
    RALLOCATE(ltransform%data,(1:2,0:ltransform%Ppts)) ! for block writing to file; 25 Mar 15;
    
    ipoint = 0
    
    RZRZt(1:5) = (/ ltransform%Ra, ltransform%Za, ltransform%R, ltransform%Z, zero /)
    
    ltransform%data(1:2,ipoint) = RZRZt(3:4)
    
    do ipoint = 1, ltransform%Ppts
     
     phistart = zero ; phiend = ( pi2 / ltransform%Nfp )
     
     id02bjf = 1
     
     Lbfieldok = .true.
     
#ifdef HAVE_NAG
     call D02BJF(phistart,phiend,Node,RZRZt(1:Node),bf00ab,ltransform%odetol,relabs,D02BJX,D02BJW,worka(1:Lworka),id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=ltransform%odetol
!    atol :  absolute tolerance
     atol=ltransform%odetol
!    do integration
     call lsode(bf00aa,(/Node/),RZRZt,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,mf=mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
     
     if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
     
     if( id02bjf.ne.0 ) exit
     
     ltransform%data(1:2,ipoint) = RZRZt(3:4)
     
    enddo
    
    ltransform%iota = RZRZt(5) / ( ltransform%Ppts * pi2 / ltransform%Nfp )

    ifail = id02bjf
    
    return
    
  end subroutine tr00aa



  
!latex \subsection{pp00aa : construct \Poincare plot (in cylindrical coordinates)}
!latex \bi
!latex \item[1.] Documentation under construction; contact \verb+shudson@pppl.gov+ if you wish to use this routine..
!latex \ei

  subroutine pp00aa( lpoincare, ifail )
#ifndef HAVE_NAG
    use lsode_mod, only: lsode 
#endif
    
    implicit none
    
    type(poincaredata) :: lpoincare
    INTEGER            :: ifail
    
    integer, parameter :: Node = 7, Lworka = 20 * Node
    integer            :: id02bjf, ipoint, astat
    real               :: RZ(1:Node), phistart, phiend
    character          :: relabs
#ifdef HAVE_NAG
    real               :: worka(1:Lworka)
    external           :: D02BJX, D02BJW
#else
    integer, parameter :: liw=20, lrw=20+16*Node
    integer            :: iwork(liw)
    real               :: rwork(lrw)
    integer            :: iopt,istate,itask,itol,mf
    real               :: atol,rtol
#endif

    ipp00aa = ifail
    
    relabs = 'D' ; itangent = 0
    
    if( allocated(lpoincare%data) ) deallocate(lpoincare%data)
    
    RALLOCATE(lpoincare%data,(1:2,0:lpoincare%Ppts))
    
    ipoint = 0
    
    RZ(1:7) = (/ lpoincare%R, lpoincare%Z, one, zero, zero, one, zero /)
    
    lpoincare%data(1:2,ipoint) = RZ(1:2)
    
    do ipoint = 1, lpoincare%Ppts
     
     phistart = zero ; phiend = ( pi2 / lpoincare%Nfp )
     
     id02bjf = 1
     
     Lbfieldok = .true.
     
#ifdef HAVE_NAG
     call D02BJF( phistart, phiend, Node, RZ(1:Node), bf00aa, lpoincare%odetol, relabs, D02BJX, D02BJW, worka(1:Lworka), id02bjf ) ! NAG ode integration;
#else
!    set the integrator parameters.
     rwork=0.
     iwork=0
!    istate=1 :  indicates the first lsode call
     istate=1
!    itask=4 :  normal integration with limited over-shoot (set by
!               rwork(1) in the i_lsode loop
     itask=1
!    iopt=1 :  optional inputs are used
     iopt=1
!    rwork(6) :  set maximum lsode-internal step size
     iwork(6)=400000
!    rwork(7) :  set minimum lsode-internal step size
!    iwork(6) :  set maximum lsode-internal steps
!    iwork(7) :  set maximum lsode-internal error messages printed
!    mf=10 :  non-stiff Adams method of integration
     mf=10
!    itol=1 :  indicates absolute tolerance is just a scalar
     itol=1
!    rtol :  relative tolerance
     rtol=lpoincare%odetol
!    atol :  absolute tolerance
     atol=lpoincare%odetol
!    do integration
     call lsode(bf00aa,(/Node/),RZ,phistart,phiend, &
                itol,(/rtol/),(/atol/),itask,istate,iopt, &
                rwork,lrw,iwork,liw,mf=mf)
     id02bjf=istate
     if (istate>0) id02bjf=0
#endif
     
     if( .not.Lbfieldok ) id02bjf = -1 ! override error flag; 05 Mar 14; an error has occured in user-supplied bfield;
     
     if( id02bjf.ne.0 ) exit
     
     lpoincare%data(1:2,ipoint) = RZ(1:2)
     
    enddo
    
    ifail = id02bjf
    
    return
    
  end subroutine pp00aa
    
   
   
end module oculus
