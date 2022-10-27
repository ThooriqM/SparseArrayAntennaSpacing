%To plot radiation pattern of Dolph-Chebyshev array
clc;
close all;
clear all;
f=94e8; % frequency
lam=3e8/f;% wavelength
d=lam/4;% separation between the elements (for d=lam/2, EFA gives bothe side
% pattern, for
% d<lam/2; one part of EFA lowers
theta=0:0.01:2*pi;
N=input('Enter the Number of Elements = ');% elements in the array
%===================== to find coeff of Chebyshev Array ===================
R0=input('Enter the side lobe level in dB = ');
R0=10^(R0/20);
Z0=cosh((1/(N-1))*acosh(R0));
aa=mod(N,2); % check N-> even or odd, if N is even then aa=0 else aa=1
if aa==0
    M=N/2;% when N is even refer balanis page 340 eq no. 6-77a and 6-77b
    for n=1:M
        amp=0;
        for q=n:M
            f1=((-1)^(M-q))*((Z0)^((2*q)-1));
            f2=factorial(q+M-2)*((2*M)-1);
            f3=factorial(q-n)*factorial(q+n-1)*factorial(M-q);
            amp=amp+f1*f2/f3;
        end
        amp1(n)=amp;
    end
else
    M=(N-1)/2 % when N is odd the equation for amp coeff is different
    for n=1:M+1
        for q=n:M+1
            f1=((-1)^(M-q+1))*((Z0)^(2*(q-1)));
            f2=factorial(q+M-2)*(2*M);
            if n==1
                epsilon=2;
            else
                epsilon=1;
            end
            f3=epsilon*factorial(q-n)*factorial(q+n-2)*factorial(M-q+1);
            amp=f1*f2/f3;
        end
        amp1(n)=amp;
    end
end

    amp2=amp1/amp1(1) ; %Normalised amplitude
%==========================================================================
%========================to find the array factor==========================
af=0;
u=pi*d*cos(theta)/lam;
if aa==0 % if N is even
    M=N/2;
    for n=1:M
        af=af+(amp2(n)*cos((2*n-1)*u));
    end
else % if N is odd
    M=(N-1)/2;
        for n=1:M
        af=af+(amp2(n)*cos(2*(n-1)*u));
    end
end
%========================to find the array factor==========================
polar(theta,af);
spacing = (1./amp2).*0.5;

s1=title('RADIATION PATTERN OF DOLPH CHEBESHEV ARRAY');
xlabel(['Having the Wavelength 1 meter'])