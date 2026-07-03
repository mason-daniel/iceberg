
    module Lib_LogisticFunction
!---^^^^^^^^^^^^^^^^^^^^^^^^^^^^
!*  
!*      Simple implementation of a logistic function
!*          f(x) = 1 / (1 + exp( -(x-x0)/w ))
!*      with width w>0
!*      x<<x0 f(x) = 0
!*      x>>x0 f(x) = 1
!*
!*      note: null constructor sets w=-1, which can be used as a check for no log func
!*
!*      Daniel Mason
!*      (c) UKAEA 2026
!*

        use iso_fortran_env
        use Lib_SafeExp
        implicit none
        private


        public      ::      LogisticFunction_ctor
        public      ::      report
        public      ::      delete

        public      ::      func
        public      ::      getx0,getw
        public      ::      isUnset

        type,public     ::      LogisticFunction
            real(kind=real64)           ::      x0,w
            real(kind=real64)           ::      iw
        end type


        interface       LogisticFunction_ctor
            module procedure    LogisticFunction_null
            module procedure    LogisticFunction_ctor0
        end interface

        interface       report
            module procedure    report0
        end interface

        interface       delete
            module procedure    delete0
        end interface

        
        interface       func
            module procedure    func0
        end interface

        interface       getx0
            module procedure    getx00
        end interface

        interface       getw
            module procedure    getw0
        end interface
        
 
 

    contains
!---^^^^^^^^

        function LogisticFunction_null() result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction)      ::      this
            this%x0 = 0
            this%w = -1
            this%iw = -1
            return
        end function LogisticFunction_null

        function LogisticFunction_ctor0(x0,w) result(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            real(kind=real64),intent(in)    ::      x0,w
            type(LogisticFunction)      ::      this
            this%x0 = x0
            this%w = w
            this%iw = 1/this%w
            return
        end function LogisticFunction_ctor0


        subroutine report0(this,u,o)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(in)       ::      this
            integer,intent(in),optional             ::      u,o
            integer         ::      uu,oo
            uu = 6 ; if (present(u)) uu = u 
            oo = 0 ; if (present(o)) oo = o 
            if (isUnset( this )) then
                write(unit=uu,fmt='(3(a,f16.6))') repeat(" ",oo)//"LogisticFunction unset"
            else
                write(unit=uu,fmt='(3(a,f16.6))') repeat(" ",oo)//"LogisticFunction [x0,w = ",this%x0,",",this%w,"]"
            end if
            return
        end subroutine report0
    
        subroutine delete0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(inout)    ::      this
            this = LogisticFunction_null() 
            return
        end subroutine delete0


    !---

        pure logical function isUnset( this )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(in)       ::      this
            isUnset = (this%iw < 0)
            return
        end function isUnset

        elemental real(kind=real64) function func0(this,x)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(in)       ::      this
            real(kind=real64),intent(in)            ::      x
            real(kind=real64)       ::      dd

            if (isUnset( this )) then
                func0 = 1
            else
                dd = (x-this%x0)*this%iw
                func0 = 1/(1+safeExp(-dd))
            end if

            return
        end function func0

        pure real(kind=real64) function getx00(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(in)       ::      this
            getx00 = this%x0
            return
        end function getx00

        pure real(kind=real64) function getw0(this)
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
            type(LogisticFunction),intent(in)       ::      this
            getw0 = this%w
            return
        end function getw0

    end module Lib_LogisticFunction