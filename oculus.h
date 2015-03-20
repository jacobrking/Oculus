
!latex \section{the user must provide the magnetic field in cylindrical coordinates}

!latex \bi
!latex \item[] \verb+subroutine bfield( RpZ, itangent, BRpZ, ifail )+\\
!latex \item[1.] \verb+RpZ(1:3)+ is real: contains the $R$, $\phi$ and $Z$ coordinates at which the field, and possibly the derivatives, are required;
!latex \item[2.] \verb+itangent+ is integer: if \verb+itangent=0+ only the magnetic field is required; if \verb+itangent=1+ then the derivatives are also required;
!latex \item[3.] \verb+BRpZ(1:3,0:3)+ is real: on output {\em must always} contain 
!latex           $B^R   \equiv{\bf B}\cdot\nabla R$,
!latex           $B^\phi\equiv{\bf B}\cdot\nabla \phi$, and 
!latex           $B^Z   \equiv{\bf B}\cdot\nabla Z$;\\
!latex           if \verb+itangent=1+, the derivatives of the field must also be provided. \\
!latex           The required format is shown:
!latex \bi
!latex \item[] \verb+BRpZ(1,0)+ $ =               B^R,   \;\;$
!latex         \verb+BRpZ(1,1)+ $ = \partial_R    B^R,   \;\;$
!latex         \verb+BRpZ(1,2)+ $ = \partial_\phi B^R,   \;\;$
!latex         \verb+BRpZ(1,3)+ $ = \partial_Z    B^R,   \;\;$
!latex \item[] \verb+BRpZ(2,0)+ $ =               B^\phi,\;\;$
!latex         \verb+BRpZ(2,1)+ $ = \partial_R    B^\phi,\;\;$
!latex         \verb+BRpZ(2,2)+ $ = \partial_\phi B^\phi,\;\;$
!latex         \verb+BRpZ(2,3)+ $ = \partial_Z    B^\phi,\;\;$
!latex \item[] \verb+BRpZ(3,0)+ $ =               B^Z,   \;\;$
!latex         \verb+BRpZ(3,1)+ $ = \partial_R    B^Z,   \;\;$
!latex         \verb+BRpZ(3,2)+ $ = \partial_\phi B^Z,   \;\;$
!latex         \verb+BRpZ(3,3)+ $ = \partial_Z    B^Z    \;\;$
!latex \ei
!latex \item[!!] Note that $B^\phi = {\bf B} \cdot \hat \phi / R$, and $\partial_R B^\phi = ( \partial_R {\bf B} \cdot \hat \phi - B^\phi ) / R$.
!latex \item[4.] \verb+ifail+ is integer: returns an error flag; if calculation of the field was successful please return \verb+ifail=0+;
!latex \ei

!latex \section{macro expansion and compilation}

!latex \bi
!latex \item[1.] The \verb+oculus.h+ file is converted to \verb+oculus.F90+ via \verb+>m4 -P oculus.macros oculus.h > oculus.F90+
!latex \item[2.] On compilation, it is required to convert single precision to double precision.
!latex \item[3.] Presently, the NAG library is required, but replacement routines will be implemented as time allows.
!latex \item[4.] Presently, the \verb+oculus.h+ and \verb+oculus.macros+ files are available at \verb+http://w3.pppl.gov/~shudson/Oculus/oculus.tar+.
!latex           At some time in the future, the routines will be kept under version control.
!latex \item[5.] Please inform \verb+shudson@pppl.gov+ of any errors; and suggestions, requests and contributions are very welcome!
!latex \ei

module oculus
  
  implicit none
  
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
  
! REAL, parameter :: pi         =  pi2 / two
! REAL, parameter :: mu0        =  2.0E-07 * pi2

  type coordinates
     INTEGER              :: Nfp
     INTEGER              :: Mpol, Ntor
     REAL                 :: Ro, Zo
     INTEGER              :: izeta, itangent
     REAL                 :: odetol, axistol, axiserror
     REAL   , allocatable :: Raxis(:,:), Zaxis(:,:)
  end type coordinates

  type(coordinates)       :: toroidalcoordinates
  
  LOGICAL              :: Lbfieldfail     
  INTEGER              :: izeta
  
contains

!latex \section{available subroutines}

!latex \newpage

!latex \subsection{ga00aa(toroidalcoordinates,ifail) : find the magnetic axis}
!latex {\bf required inputs}
!latex \bi 
!latex \item[  ] The user must include \verb+"use oculus, only : toroidalcoordinates, ga00aa"+ in their source that calls \verb+ga00aa+.
!latex           The following must be set:
!latex \item[1.] \verb+toroidalcoordinates%Nfp     : integer + the toroidal periodicity of the magnetic field, e.g. \verb+Nfp=1+;
!latex \item[2.] \verb+toroidalcoordinates%Ntor    : integer + the toroidal Fourier resolution that the magnetic axis, e.g. +Ntor=0+;
!latex \item[3.] \verb+toroidalcoordinates%Ro      : real    + an initial guess for the $R$ location of the magnetic axis on the $\phi=0$ plane;
!latex \item[4.] \verb+toroidalcoordinates%Zo      : real    + an initial guess for the $Z$ location of the magnetic axis on the $\phi=0$ plane;
!latex \item[5.] \verb+toroidalcoordinates%axistol : real    + required accuracy, e.g. \verb+axistol=1.0E-06+;
!latex \item[6.] \verb+toroidalcoordinates%odetol  : real    + magnetic fieldline o.d.e. integration tolerance; e.g. \verb+odetol=1.0e-08+;
!latex \item[7.] \verb+ifail                       : integer + error control flag, e.g. \verb+ifail=0+; quiet mode is \verb+ifail=1+;
!latex \ei
!latex {\bf outputs}
!latex \bi
!latex \item[] \verb+toroidalcoordinates%Ro        : real   + updated;
!latex \item[] \verb+toroidalcoordinates%Zo        : real   + updated;
!latex \item[] \verb+toroidalcoordinates%tangent   : real   + the tangent mapping at axis; under construction;
!latex \item[] \verb+toroidalcoordinates%Rm        : real   + the Fourier harmonics, $R(\phi)=\sum_m R_m \cos(m\phi)$; under construction;
!latex \item[] \verb+toroidalcoordinates%Zm        : real   + the Fourier harmonics, $Z(\phi)=\sum_m Z_m \sin(m\phi)$; under construction;
!latex \item[] \verb+toroidalcoordinates%axiserror : real   + the error;
!latex \item[] \verb+ifail                         : integer+ on normal execution \verb+ifail=0+; can also check that \verb+axiserror+ is small;
!latex \ei
!latex {\bf method}
!latex \bi 
!latex \item[] From the supplied $(R,Z)$, fieldline tracing methods are used to find the fieldline that closes on itself after a toroidal distance of $2\pi$/\verb+Nfp+.
!latex \item[] The NAG routine \verb+C05PBF+ is used for the nonlinear root find, and the NAG routine \verb+D02BJF+ is used for the o.d.e. integration.
!latex \item[] The \verb+ifail+ flag is passed directly to \verb+C05PBF+: please see the NAG documentation for how this flag is treated.
!latex \ei
  
  subroutine ga00aa( toroidalcoordinates, ifail )
    
    implicit none
    
    type(coordinates)  :: toroidalcoordinates
    
    INTEGER            :: ifail
    
    INTEGER, parameter :: NN = 2, Ldfjac = NN, Lrwork = NN * ( 3 * NN + 13 ) / 2
    INTEGER            :: ic05pbf
    REAL               :: tol, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), rwork(1:Lrwork)
    CHARACTER          :: message*21

    if( toroidalcoordinates%Ro      .le. zero .or. &
        toroidalcoordinates%axistol .le. zero ) then
     ic05pbf = 1
     goto 9000
    endif

    ic05pbf = ifail ; tol = toroidalcoordinates%axistol
    
    RZ(1:2) = (/ toroidalcoordinates%Ro, toroidalcoordinates%Zo /)
    
    call C05PBF( ga00ab, NN, RZ(1:NN), FRZ(1:NN), dFRZ(1:Ldfjac,1:NN), Ldfjac, tol, rwork(1:Lrwork), Lrwork, ic05pbf )
    
    select case( ic05pbf )
    case( 0 )    ; toroidalcoordinates%Ro = RZ(1) ; toroidalcoordinates%Zo = RZ(2) ; toroidalcoordinates%axiserror = sqrt(sum(FRZ(1:NN)**2))
    case default
    end select
    
9000 continue
    
    if( ifail.le.0 ) then ! shall print messages to screen; 02 Mar 15;
     
     select case( ic05pbf ) !01233456789012345678901
     case( 0 )    ; message = "success ;         "
     case( 1 )    ; message = "input error ;     "
     case default ; message = "NAG C05PBF error ;"
     end select
     
     write(*,1000) ic05pbf, toroidalcoordinates%Ro, toroidalcoordinates%Zo, toroidalcoordinates%axistol, toroidalcoordinates%axiserror, message
     
    endif
    
    ifail = ic05pbf
    
    return
    
1000 format("ga00aa : ifail="i3" ; ( Ro, Zo ) = ("es23.15" ,"es23.15" ) ; axistol="es10.3" ; axiserror=",es10.3" ; "a21)
    
  end subroutine ga00aa

    
  subroutine ga00ab( NN, oRZ, FRZ, dFRZ, Ldfjac, iflag )
    
    implicit none
    
    INTEGER, intent(in)    :: NN
    REAL                   :: oRZ(1:NN), FRZ(1:NN)
    
    INTEGER, intent(in)    :: Ldfjac
    REAL                   :: dFRZ(1:Ldfjac,1:NN)
    
    INTEGER, intent(inout) :: iflag
    
    INTEGER, parameter     :: Node = 6, Lrwork = 20 * Node
    
    INTEGER                :: id02bjf
    REAL                   :: dRZ(1:Node), zetastart, zetaend, tol, rwork(1:Lrwork), determinant, residue
    CHARACTER              :: relabs
    
    external               :: D02BJW
    
    tol = toroidalcoordinates%odetol ; relabs = 'D'
    
    zetastart = zero ; zetaend = zetastart + ( pi2 / toroidalcoordinates%Nfp ) ! integration endpoints ; 05 Mar 14;
    
    dRZ(1:Node) = (/ oRZ(1), oRZ(2), one, zero, zero, one /) ! initial guess, intialize tangent map integration; 31 Jul 13;
    
    izeta = 0 ! this counter is incremented in ga00ad; 31 Jul 13;
    
    id02bjf = 1

    select case( iflag )
    case( 0 ) ; toroidalcoordinates%itangent = 0 ! derivatives are not required; 02 Mar 15;
    case( 1 ) ; toroidalcoordinates%itangent = 1
    end select
    
    call D02BJF( zetastart, zetaend, Node, dRZ(1:Node), bf00aa, tol, relabs, ga00ac, D02BJW, rwork(1:Lrwork), id02bjf ) ! NAG ode integration;
  
    select case( id02bjf )
    case( -1 ) ! write(*,'("ga00ac : " 10x " : id02bjf="i3" : user termination ; RZ="2es23.15" ;            ")') id02bjf, oRZ(1:2)
    case(  0 ) ! write(*,'("ga00ac : " 10x " : id02bjf="i3" : success ; RZ="2es23.15" ;                     ")') id02bjf, oRZ(1:2)
    case(  1 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : input error ;                                 ")') id02bjf! oRZ(1:2)
    case(  2 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : no further progress possible ; RZ="2es23.15" ;")') id02bjf, oRZ(1:2)
    case(  3 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : tol too small ;                               ")') id02bjf! oRZ(1:2)
    case(  4 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : xsol not reset ;                              ")') id02bjf! oRZ(1:2)
    case(  5 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : xsol not reset ;                              ")') id02bjf! oRZ(1:2)
    case(  6 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : function did not change sign ;                ")') id02bjf! oRZ(1:2)
    case(  7 ) ; write(*,'("ga00ac : " 10x " : id02bjf="i3" : serious error ;                               ")') id02bjf! oRZ(1:2)
    end select
    
    if( Lbfieldfail ) id02bjf = -1 ! override error flag; 05 Mar 14;
    
    if( id02bjf.eq.0 ) then ! o.d.e. integration was successful; 02 Mar 15;
     
    !determinant = dRZ(3)*dRZ(6) - dRZ(4)*dRZ(5) ; residue = ( two - (dRZ(3)+dRZ(6)) ) / four
     
     select case( iflag ) 
     case( 1 ) ;  FRZ(1  ) = dRZ(1) - oRZ(1)                            ! must return function;  5 Jun 13;
      ;        ;  FRZ(2  ) = dRZ(2) - oRZ(2)
     case( 2 ) ; dFRZ(1,1) = dRZ(3) - one    ; dFRZ(1,2) = dRZ(4)       ! must return Jacobian;  5 Jun 13;
      ;        ; dFRZ(2,1) = dRZ(5)          ; dFRZ(2,2) = dRZ(6) - one
     end select
     
    !AxisTangent(1,1,iaxis) = dRZ(3) ; AxisTangent(1,2,iaxis) = dRZ(4) ! tangent mapping, returned through global; 05 Mar 14;
    !AxisTangent(2,1,iaxis) = dRZ(5) ; AxisTangent(2,2,iaxis) = dRZ(6) 
     
    else
     
     iflag = -1 ! tell NAG that an error has occured;  5 Jun 13;
     
    endif
    
    return
    
  end subroutine ga00ab

  
  subroutine ga00ac( zeta, RZ )
    
    implicit none
    
    INTEGER, parameter :: Node = 6
    
    REAL               :: zeta, RZ(1:Node)
    
   !Raxis(izeta,iaxis) = RZ(1) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
   !Zaxis(izeta,iaxis) = RZ(2) ! save information along magnetic axis in preparation for Fourier decomposition; 05 Mar 14;
    
    izeta = izeta + 1

    zeta = izeta * ( pi2 / toroidalcoordinates%Nfp ) / ( 4 * max(toroidalcoordinates%Ntor,1) )
    
    return
    
  end subroutine ga00ac

  
  subroutine bf00aa( zeta, RZ, BRZ ) ! the format of this subroutine is constrained by the NAG ode integration routines. 
    
    INTEGER, parameter  :: Node = 6
    
    REAL  , intent(in)  :: zeta, RZ(1:Node)
    REAL  , intent(out) :: BRZ(1:Node)
  
    INTEGER             :: ibfield
    REAL                :: RpZ(1:3), dBRpZ(1:3,0:3)
    
    REAL                :: TM(1:2,1:2)
    
    RpZ(1:3) = (/ RZ(1), zeta, RZ(2) /)
    
    call bfield( RpZ(1:3), toroidalcoordinates%itangent, dBRpZ(1:3,0:3), ibfield )
  
    select case( ibfield ) 
    case( 0 )    ; Lbfieldfail = .false.
    case( 1 )    ; Lbfieldfail = .true.  ; BRZ(1:Node) = zero ; goto 9999
    case default ; write(*,'("bf00aa : illegal ibfield returned from bfield ;")') 
    end select
    
    BRZ(1:2) = (/ dBRpZ(1,0), dBRpZ(3,0) /) /  dBRpZ(2,0)
    
    TM(1,1) = ( dBRpZ(1,1) - BRZ(1) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(1,2) = ( dBRpZ(1,3) - BRZ(1) * dBRpZ(2,3) ) / dBRpZ(2,0)
    TM(2,1) = ( dBRpZ(3,1) - BRZ(2) * dBRpZ(2,1) ) / dBRpZ(2,0)
    TM(2,2) = ( dBRpZ(3,3) - BRZ(2) * dBRpZ(2,3) ) / dBRpZ(2,0)
    
    BRZ(3) = TM(1,1) * RZ(3) + TM(1,2) * RZ(5) ! tangent map obtained by matrix multiplication;
    BRZ(4) = TM(1,1) * RZ(4) + TM(1,2) * RZ(6)
    BRZ(5) = TM(2,1) * RZ(3) + TM(2,2) * RZ(5)
    BRZ(6) = TM(2,1) * RZ(4) + TM(2,2) * RZ(6)
    
9999 continue
    
    return
    
  end subroutine bf00aa

 
end module oculus
