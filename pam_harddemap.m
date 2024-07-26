function b = pam_harddemap(y, ldM)
  M = pow2(ldM);
  y = y * sqrt(2*(M^2-1)/3);    % scale to {-(M-1), ...,-1,1,3, ..., M-1 }

  br = y > 0; 
  y = abs(y);
  if ldM == 2
    b = [br; y < 2];
  elseif ldM == 3
    b = [br; y<4; abs(y-4)<2];
  elseif ldM == 4
      b = [ br; y<8; abs(y-8)<4; (abs(y-4)<2) | (abs(y-12)<2) ];
  else
    error('unknown ldM')
  end
end
