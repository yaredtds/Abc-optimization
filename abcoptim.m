function [vect, cost] = abcoptim(costf, vlen, pop, iter, maxtr)

  % usage:
  %   [vect, cost] = fminabc(costf, vlen, pop, iter, maxtr)
  %    vect - return value
  %   costf - cost function
  %    iter - total iteration
  %   maxtr - maximum number of trials
  %     pop - population size
  %    vlen - length of a vector
  pop = pop + 1;
  [vmax, ~] = costf(ones(1, vlen).*1e3);
  [vmin, ~] = costf(-ones(1, vlen).*1e3);
  vmid = (vmax + vmin)./2;
  vstd = (vmax - vmin)./4;
  state = ['e', 'o', 's'];
  
  %hist.c = []; % store cost value per iteration
  %hist.i = []; % store iteration
  
  b = [];
  b(1).v = vmid + randn(1, vlen).*vstd;
  [b(1).v, b(1).c] = costf(b(1).v);
  b(1).s = state(randi([1,3]));
  b(1).t = 0;
  
  
  for i=2:pop
      b(i).v = vmid + randn(1, vlen).*vstd;
      [b(i).v, b(i).c] = costf(b(i).v);
      b(i).s = state(randi([1,2]));
      b(i).t = 0;
      
      m = i-1;
      tmp = b(i);
      while (m > 0 && tmp.c < b(m).c)
          b(m+1) = b(m);
          m = m - 1;
      end
      b(m+1) = tmp;
      
  end
  % best of current solution
  
  % exploit + exploring phase
  for j=1:iter
      for k=2:pop
          % Desired cost value
          dc = b(k).c + rand()*(b(1).c-b(k).c);
          
          % Modified binary search
          lh = 1; rh = k;
          while lh <= rh && b(lh).c <= dc
              m = floor((lh + rh)/2);
              if b(m).c == dc
                  break;
              elseif b(m).c < dc
                  lh = m + 1;
              else
                  rh = m - 1;
              end %@if
          end %@while
          
          pb = b(randi([1 m]));
          
          switch b(k).s
              case 'e'
                  % vs will be zero when b choose itself
                  % updating whole element increase/decrease cost
                  % *1. update element of maximum deviation
                  % *2. sequentially update the whole vect element
                  %# -> increasing vect element dont minimize cost, try
                  % increasing else try average value
                  % -> At this point b(k) has two choice, either to choose
                  % neighbour soln which is better than it or take any pose
                  % better than current nhbr solution
                  for n=1:vlen
                      si = randn(1);
                      wh1 = exp(pb.c - b(k).c);
                      va = (b(k).v(n) * wh1 + pb.v(n))/(1 + wh1);
                      vs = (b(k).v(n) - pb.v(n))./4;
                      
                      vp = b(k).v;
                      vn = b(k).v;
                      dv = si.*vs;
                      
                      vp(n) = va + dv; [vp, cp] = costf(vp);
                      vn(n) = va - dv; [vn, cn] = costf(vn);
                      
                      if (cp < cn && cp < b(k).c)
                          b(k).c = cp; b(k).v = vp; b(k).t = 0; b(k).s = 'o';
                      elseif (cn < b(k).c)
                          b(k).c = cn; b(k).v = vn; b(k).t = 0; b(k).s = 'o';
                      else
                          %b(k) = pb;
                          b(k).t = b(k).t + 1;
                          if b(k).t > maxtr
                              b(k).s = 's';
                          end
                      end
                      % update upper & lower limits                      
                  end
                  
              case 'o'
                  for n=1:vlen
                      si = randn(1,2);
                      wh = exp([b(1).c-b(k).c  b(1).c-pb.c  1]);
                      va = ([b(k).v(n) pb.v(n) b(1).v(n)]*wh')/(1 + wh(1) + wh(2));
                      vs = [b(k).v(n) - b(1).v(n)  pb.v(n) - b(1).v(n)]./4;
                      
                      vp = b(k).v;
                      vn = b(k).v;
                      dv = si.*vs;
                      
                      vp(n) = va + dv(1) + dv(2); [vp, cp] = costf(vp); %++ +- -+ --
                      vn(n) = va - dv(1) - dv(2); [vn, cn] = costf(vn);
                      if (cp < cn) && (cp < b(k).c)
                          b(k).c = cp; b(k).v = vp; b(k).t = 0; b(k).s = 'e';
                      elseif (cn < b(k).c)
                          b(k).c = cn; b(k).v = vn; b(k).t = 0; b(k).s = 'e';
                      else   
                          %b(k) = pb;
                          b(k).t = b(k).t + 1;
                          if b(k).t > maxtr
                              b(k).s = 's';
                          end
                      end
                  end
                  
              case 's'
                  b(k).v = vmid + randn(1, vlen).*vstd;
                  [b(k).v  b(k).c] = costf(b(k).v);
                  b(k).t = 0; b(k).s = 'e';
          end
          
          % Update global optima
          % buble sort the location
          m = k-1; tmp = b(k);
          while (m>0 && tmp.c < b(m).c)
              b(m+1) = b(m);
              m = m - 1;
          end
          b(m+1) = tmp;
          
      end
      %hist.c(j)= b(1).c;
      %hist.i(j) = j;
  end % End Main Loop
  vect = b(1).v; cost=b(1).c;
  %range = sprintf('%c%d:%c%d',xlsc,2,xlsc,j+1);
  %xlswrite('SimResult.xlsx',[hist.c'],'Sheet-1',range);
end % End of fminabc
