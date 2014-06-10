function x = revised_simplex(A,b,c)
% x = revised_simplex(A,b,c)
% Implement the revised_simplex algorithm for a problem of the form
% minimize  c'x
% such that Ax = b 
% and        x>=0 , b>=0
    c = c(:);
    b = b(:);
    [M,N] = size(A);
    % First we need to solve the augmented LP to get our first BFS
    % minimize 1'y
    sprintf('Solution to part 1: \r\n')
    [x,bv] = find_bfs(A,b)
    dv = setdiff([1:N],bv);

    lambda = linsolve(A(:,bv).',c(bv));
    r = c(dv).'-lambda.'*A(:,dv);

    while ~isempty(find(r<0));
        J = find(r<0);
        JIDX = J(1);
        J = dv(JIDX); % basis to enter 
        ratios = b./A(:,J);
        if isempty(find(ratios>=0))
            error('Initial feasible solution could not be found. Phase I system is unbounded from below.');
        else
            min_ratios=min(ratios(find(ratios>=0)));
            K=find(ratios==min_ratios);
            KIDX = K(1);
            K=bv(KIDX); % basis to leave
            bv(KIDX) = J;
            dv(JIDX) = K;

            b = linsolve(A(:,bv),b);
            lambda = linsolve(A(:,bv).',c(bv));
            r = c(dv).'-lambda.'*A(:,dv);

        end

    end

    x = zeros(1,N);
    x = x(:);
    x(bv) = b;

end


function [x,bv] = find_bfs(A,b)
  b = b(:)
  [M,N] = size(A);
  B = eye(M);
  A=[A,B];
  bv = N+(1:M);
  c = zeros(1,M+N);
  c = c(:);
  c(bv) = 1;
  dv = setdiff([1:N+M],bv);
  lambda = linsolve(A(:,bv).',c(bv));
  r = c(dv).'-lambda.'*A(:,dv)
  while ~isempty(find(r<0));
      J = find(r<0);
      JIDX = J(1);
      J = dv(JIDX); % basis to enter 
      ratios = b./A(:,J);
      if isempty(find(ratios>=0))
          error('Initial feasible solution could not be found. Phase I system is unbounded from below.');
      else
          min_ratios=min(ratios(find(ratios>=0)));
          K=find(ratios==min_ratios);
          KIDX = K(1);
          K=bv(KIDX); % basis to leave
          sprintf('basis to enter: %d basis to leave: %d\r\n',J,K)
          bv(KIDX) = J;
          dv(JIDX) = K;

          b = linsolve(A(:,bv),b);
          lambda = linsolve(A(:,bv).',c(bv));

          r = c(dv).'-lambda.'*A(:,dv)
          
      end

  end
  x = zeros(1,N);
  x = x(:);
  x(bv) = b;
        
end
