
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Math 4220 Fall 2018 
% Assignment 6 Solutions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% It is generally recommended to clear all previously defined
% variables from the workspace.
clear; 
clc;

format compact;

% Toggle between problems

index=1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 2

if index == 1
    
  % Interval  
  a=0;
  b=1;
  
  for isw=1:2
      
      disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
      if isw == 1
          disp('Results for Part (a)');
          I=pi;
      else
          disp('Results for Part (b)');
          I=2/3;
      end
      fprintf('I_trap   e_trap     I_mid     e_mid     I_simp     e_simp \n');
      
      for i=1:5
          
          % Number of Subintervals
          r=2^i;
          h=(b-a)/r;
          
          % Composite Trapezoidal Method
          x=a:h:b;
          y=funex2(isw,x);
          I_trap(i)=(h/2)*(y(1)+2*sum(y(2:end-1))+y(end));
          e_trap(i)=I-I_trap(i);
          
          % Composite Midpoint Method
          x=a+h/2:h:b-h/2;
          y=funex2(isw,x);
          I_mid(i)=h*sum(y);
          e_mid(i)=I-I_mid(i);
          
          % Composite Simpson's Method
          x=a:h:b;
          y=funex2(isw,x);
          I_simp(i)=(h/3)*(y(1)+4*sum(y(2:2:end-1))+2*sum(y(3:2:end-2))+y(end));
          e_simp(i)=I-I_simp(i);
          
          fprintf('%6.4f   %6.2e   %6.4f   %6.2e   %6.4f   %6.2e\n',I_trap(i),e_trap(i),I_mid(i),e_mid(i),I_simp(i),e_simp(i))
          
      end
      
  end
  
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 3

if index == 2
    
    N=10;
    
    c=10*ones(N+1,N+1);
    g=11*ones(N+1,N);
    
    for j=0:N
        
        if j == 0
            c(1,1)=1;
        elseif j == 1
            c(2,1)=1;
            c(2,2)=0;
            g(2,1)=roots(c(2,1:2));
        else
            c(j+1,1:j+1) = (2*j-1)/(j)*[c(j,1:j) 0]-(j-1)/(j)*[0 0 c(j-1,1:j-1)];
            g(j+1,1:j) = roots(c(j+1,1:j+1));
        end
        
    end
       
    % The Gaussian points for n, n=0:9, are now in the first 
    % n+1 entries of the (n+2)nd row of g
    
    for j=1:N
        plot([-0.99 0.99],[j,j],'k')
        hold on 
        plot(g(j+1,1:j),j*ones(1,j),'o')
    end
    axis([-1 1 0 N+1])
    hold off
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problem 5

if index == 3
    
    format long
    
    func=@(x) 4./(1+x.^2);
    a=0;
    b=1;
    k=5;
    Rtab=romberg(func,a,b,k);
    disp(Rtab);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function Definitions

% Problem 2
function f=funex2(isw,x)
    
    if (isw==1)
        f=4./(1+x.^2);
    else
        f=sqrt(x);
    end

end

% Problem 5
function Rtab=romberg(fname,a,b,k)
    % romberg(fname,a,b,k) computes the k by k Romberg table to provide 
    % approximations to the definite integral of fname over the interval [a,b]
    
    % Zero the k by k Rtab (romberg table)
    Rtab=zeros(k);
    h=b-a;
    
    % Compute Rtab11
    Rtab(1,1)=h*(fname(a)+fname(b))./2.0;
    
    % Compute the remainder of the Romberg table
    for j=2:k
        % Compute the new Composite Trapezoidal approx with h/2
        Rtab(j,1)=Rtab(j-1,1) ;
        for i=1:2^(j-2)
            Rtab(j,1)=Rtab(j,1)+h*fname(a+(i-0.5).*h);
        end
        Rtab(j,1)=Rtab(j,1)./2.0;
        
        % Compute the remaining columns of the table's jth row.
        for r=2:j
            Rtab(j,r)=Rtab(j,r-1)+(Rtab(j,r-1)-Rtab(j-1,r-1))./(4.0.^(r-1)-1);
        end
        h=h/2;
    end
end