function Jac=NextGenJ(X,P,S,freq)

% Jacobian for next-generation NMM

g=X(1); x=X(3); y=X(4); A=X(5); tau=X(7);

Jac=[0, 1, 0, 0, 0, 0, 0;
      -P.a^2, -2*P.a, P.a^2*P.k*dfx(x,y), P.a^2*P.k*dfy(x,y), 0 , 0, 0;
      dxdg(x,y,P.v), 0, drfdx(x,y,P.n+A,P.d)+drgdx(x,y,g,P.v), drfdy(x,y,P.n+A,P.d)+drgdy(x,y,g,P.v), dxdA(x,y), 0, 0;
      dydg(x,y,P.v), 0, difdx(x,y,P.n+A,P.d)+digdx(x,y,g,P.v), difdy(x,y,P.n+A,P.d)+digdy(x,y,g,P.v), dydA(x,y), 0, 0;
      0, 0, 0, 0, 0, 1, 0;
      0, 0, 0, 0, -P.ad^2, -2*P.ad, P.ad^2*(dI(tau,S,freq));
      0, 0, 0, 0, 0, 0, 0];
  
end
  
    function X = dfx(x,y)
        X=-(1/pi)*2*(x.^2+2*x-y.^2+1)./((1+2*x+x.^2+y.^2).^2);
    end

    function X = dfy(x,y)
        X=-(1/pi)*4*y.*(1+x)./((1+2*x+x.^2+y.^2).^2);
    end

    function X=dxdg(x,y,v)
        X=-v*y*(x+1)-(x^2-y^2-1)/2;
    end

    function X=dydg(x,y,v)
        X=(x^2-y^2+2*x+1)*v/2-x*y;
    end

    function X=drfdx(x,y,n,d)
        X=y-d*(x+1)-y*n;
    end

    function X=drgdx(x,y,g,v)
        X=-y*v*g-x*g;
    end

    function X=drfdy(x,y,n,d)
        X=x-1+d*y-n*(x+1);
    end

    function X=drgdy(x,y,g,v)
        X=-(x+1)*v*g+y*g;
    end

    function X=difdx(x,y,n,d)
        X=-x+1-d*y+n*(x+1);
    end

    function X=digdx(x,y,g,v)
        X=(x+1)*v*g-y*g;
    end

    function X=difdy(x,y,n,d)
        X=y-d*(x+1)-n*y;
    end

    function X=digdy(x,y,g,v)
        X=-y*v*g-x*g;
    end

    function X=dxdA(x,y)
        X=-x*y-y;
    end

    function X=dydA(x,y)
        X=(x^2-y^2+2*x+1)/2;
    end

    function X=dI(t,S,f)
        X=f*pi*S*cos(f*2*pi*t);
    end