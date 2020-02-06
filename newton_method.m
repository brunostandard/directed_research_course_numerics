function y = newton_method(f,df,y0,k);
  u(:,1)  = y0;
  tol = 0.001;
  for j = 1:size(y0,2)
    u = [u , zeros(size(y0,1),1) ];
    u(:,j+1) = u(:,j) - f(u(:,j))./df(u(:,j));
    if abs(u(:,j)-u(:,j)) < tol
      break;
    end
  end
  y = u(:,end-1);
end
