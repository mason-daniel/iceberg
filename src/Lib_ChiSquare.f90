
    module Lib_ChiSquare
!---^^^^^^^^^^^^^^^^^^^^
!*      A simple implementation of a chi-square test
!*      usage:
!*          alpha = chi_square_test( observed_dat,expected_dat )
!*              given the observed data, return the alpha score 
!*              ie the prob that a random sample would give a higher chi-square value 
!*              assuming 1 d histogram data so ndof = size(dat) - 1
!*
!*      Daniel Mason
!*      (c) UKAEA June 2024



        use lib_linint
        use iso_fortran_env

        implicit none
        private

        include 'chi_square.h'

        public          ::      chi_square
        public          ::      chi_square_test
        public          ::      alpha_score

        interface   chi_square
            module procedure        chi_square0
        end interface

        interface   chi_square_test
            module procedure        chi_square_test0
        end interface

    contains
!---^^^^^^^^

        pure real(kind=real64) function chi_square0( observed_dat,expected_dat )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the observed data, return the chi_square score 
    !*      assumes 1 d histogram data so ndof = size(dat) - 1
    !*      use Yates's correction
    !*          X^2  = sum_i ( |O_i - E_i| - 0.5 )^2/E_i 
            integer,dimension(:),intent(in)                 ::      observed_dat
            real(kind=real64),dimension(:),intent(in)       ::      expected_dat
    
            real(kind=real64)               ::      sum_observed

            integer                         ::      ii,nn
            real(kind=real64)               ::      xx 
            real(kind=real64),dimension(size(observed_dat))     ::      normalised_expected_dat

            

            nn = size(observed_dat)

            sum_observed = sum(observed_dat)

            normalised_expected_dat = max(1.0d-16,expected_dat)     


            if (sum_observed == 0) then
                chi_square0 = 0.25d0 * sum_observed      
            else
                chi_square0 = 0.0d0
                do ii = 1,nn
                    xx = observed_dat(ii) - normalised_expected_dat(ii)              
                    chi_square0 = chi_square0 + ( abs(xx)-0.5d0 )**2 / normalised_expected_dat(ii)   
                end do
            end if
            return
        end function chi_square0

        pure real(kind=real64) function chi_square_test0( observed_dat,expected_dat,ndof )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given the observed data, return the alpha score ( prob that a random sample would give a higher chi-square value )
    !*      assumes 1 d histogram data  
            integer,dimension(:),intent(in)                 ::      observed_dat
            real(kind=real64),dimension(:),intent(in)       ::      expected_dat
            integer,intent(in),optional                     ::      ndof
    
            real(kind=real64)               ::      chi_square
            integer                         ::      df            

            df = size(observed_dat)         !   assume sum of expected distribution is not sum of observed distribution
            if (present(ndof)) df = ndof    !   unless you've done something tricksy to the expected
            
            chi_square = chi_square0( observed_dat,expected_dat )
            chi_square_test0 = alpha_score( chi_square, df )
            return
        end function chi_square_test0


        pure real(kind=real64) function alpha_score( chi_square, df )
    !---^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
    !*      given a value of the chi-square statistic and the number of degrees of freedom,
    !*      compute the alpha- score, ie the probability that we might randomly find a value of chi-square larger than 
    !*      the value observed
            real(kind=real64),intent(in)            ::      chi_square
            integer,intent(in)                      ::      df

            real(kind=real64),dimension(size(CHI_SQ_TABLE,dim=1))    ::      right_tail_prob

            integer             ::      ii,jj,nAlpha,nDof
            real(kind=real64)   ::      aa

            nAlpha = size(ALPHA)
            nDof = size(DOF)
 

        !---    find the right tail probability for the required number of dof

            if (df>DOF(nDof)) then

                !   use a simple extrapolation to high dof.
                do jj = 1,nAlpha
                    if (ALPHA(jj)>0.5d0) then
                        right_tail_prob(jj) = (-13.382*ALPHA(jj)*ALPHA(jj) + 24.208*ALPHA(jj) - 10.065)*df   &
                                            + (-435.71*ALPHA(jj)*ALPHA(jj) + 787.72*ALPHA(jj) - 362.32)
                    else
                        right_tail_prob(jj) = (-0.039*log(ALPHA(jj)) + 1.0313)*df   &
                                            + (-3.427*log(ALPHA(jj)) - 0.5687)
                    end if
                end do

            else     

                do ii = 1,nDof-1
                    
                    if (df > DOF(ii+1)) then
                        cycle     
                    else if (df == DOF(ii+1)) then
                        !   quite a likely scenario - DOFs are enumerated 1,2,3..30.
                        right_tail_prob(:) = CHI_SQ_TABLE(:,ii+1)
                        exit
                    else 
                        !   I know DOF(ii) < df < DOF(ii+1) , and DOF is an ordered list
                        aa = real(df-DOF(ii),kind=real64) / (DOF(ii+1)-DOF(ii))             !   fraction of way across interval
                        
                        do jj = 1,nAlpha
                            right_tail_prob(jj) = linint( CHI_SQ_TABLE(jj,:),ii-1,aa )      !   linint expects i=0:N-1
                        end do
                        exit
                    end if
                end do
            end if
           

        !---    now find the value of alpha corresponding to the chi-squared stat
            if (chi_square < right_tail_prob(1)) then
                alpha_score = ALPHA(1)
            else if (chi_square > right_tail_prob(nAlpha)) then
                alpha_score = ALPHA(nAlpha)
            else
                do jj = 1,nAlpha-1
                    if (chi_square > right_tail_prob(jj+1)) then
                        cycle
                    else                         
                        aa = (chi_square-right_tail_prob(jj)) / (right_tail_prob(jj+1)-right_tail_prob(jj))             !   fraction of way across interval
                        alpha_score = ALPHA(jj) + (ALPHA(jj+1)-ALPHA(jj))*aa
                        exit
                    end if
                end do
            end if
            
            return
        end function alpha_score
 
    end module Lib_ChiSquare



