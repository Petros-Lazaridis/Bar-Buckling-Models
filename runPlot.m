%Created by Petros Lazaridis (2020)
clear;

L = [ 112 84 154];
D = [ 14 14 14];
eps_y = [0.002, 0.003, 0.004];
eps_sh = [1, 1, 1].*eps_y;
eps_u = [60 38.3 30].*eps_y;
fy = [400, 540, 800];
fu = [600, 637.2, 1000];
eps = [zeros(3,1), linspace(eps_y, eps_u, 100)];
failStress = 0.85*fy
E_sh = zeros(length(fy),1)

epsilon85fy = zeros(length(L), length(fy))
for steel = 1:length(fy)
 for LDind = 1:length(D)
   
  [sigmaT(:,LDind,steel), sigmaDM(:,LDind,steel), f_sc(:,LDind,steel)] =...
  BucklingModels(eps(steel,:), L(LDind), D(LDind), eps_y(steel), eps_sh(steel),...
  eps_u(steel), fy(steel), fu(steel) );
  [maxf_sc, indMax] = max(f_sc(:,LDind,steel))
  epsilon85fy(LDind,steel) = interp1(f_sc(indMax:end,LDind,steel), eps(steel, indMax:end), failStress(steel), 'linear')
  
 end
 yield_point = find( sigmaT(:, LDind,steel) == fy(steel)) %find(eps(:,steel) == eps_y(steel))
 yield_point = yield_point(end)
 E_T(steel,:) = (sigmaT(3:end,LDind,steel) - sigmaT(1:end-2,LDind,steel))./(eps(steel,3:end) - eps(steel,1:end-2))'
 E_sh(steel) = E_T(steel, yield_point+1)
end

eps = eps*1000
epsilon85fy = epsilon85fy*1000
eps_y = eps_y*1000
eps_u = eps_u*1000
figure
set(gcf, 'Units','centimeters', 'Position',[0 0 15 18.5])
for cSS = 1:steel
  
  subplot(steel,2, (cSS-1)*2+1)
  plot(eps(cSS,:), sigmaDM(:,:,cSS), 'linewidth', 1.8, 'color', [1 0 0])
  %text(eps(cSS,end)*ones(length(L),1), sigmaDM(end,:,cSS), 'DM', 'Fontsize', 10.3, 'Fontname', 'Times')
  hold
  plot(eps(cSS,:), f_sc(:,:,cSS), 'linewidth', 1.8, 'color', [0 1 0])
  %text(eps(cSS,end)*ones(length(L),1), f_sc(end,:,cSS), 'RDM', 'Fontsize', 10.3, 'Fontname', 'Times')
  plot(eps(cSS,:), sigmaT(:,1,cSS), 'linewidth', 1.8, 'color', [0 0 0 ])
  plot([0 max(epsilon85fy(:,cSS))], failStress(cSS)*ones(1,2), 'linestyle', '--','linewidth', 0.7, 'color', [0 0 0 ])
  plot(epsilon85fy(:,cSS).*ones(3,2),[0 ;failStress(cSS)],'linestyle', '--','linewidth', 0.7, 'color', [0 0 0 ])
  text(epsilon85fy(:,cSS)'*1.05,0.2*failStress(cSS)*[1, .6, .5], num2str(epsilon85fy(:,cSS), '%2.1f'), 'color', [0 0 0 ], 'Fontsize', 8.3, 'Fontname', 'Times' )
  text(eps(cSS,50)*ones(length(L),1), (sigmaDM(50,:,cSS) + f_sc(50,:,cSS))/2,...
  [num2str(L'./D')], 'Fontsize', 10.3, 'Fontname', 'Times')
  text(1.4*eps(yield_point), 1.2*sigmaT(yield_point,1,cSS), ['E_{sh} = ', num2str(E_sh(cSS), '%4.2f'),' MPa'],...
  'Fontsize', 10.3, 'Fontname', 'Times')
  text(-0.12*eps_u(cSS), failStress(cSS), num2str(failStress(cSS)), 'Fontsize', 10.5, 'Fontname', 'Times')
  grid
  %title(['f_y=', num2str(fy(cSS)), ' MPa',' ', 'f_u=', num2str(fu(cSS)), ' MPa',' ',...
  %'\epsilon_y=', num2str(eps_y(cSS)), ' ', '\epsilon_{sh}=', num2str(eps_sh(cSS)),' ','\epsilon_u=', num2str(eps_u(cSS))])
  title(['f_y=', num2str(fy(cSS)), ' MPa',' ', 'f_u=', num2str(fu(cSS)), ' MPa',' ',...
  '\epsilon_y=', num2str(eps_y(cSS)), '^o/_{oo}', ' ','\epsilon_u=', num2str(eps_u(cSS)), '^o/_{oo}'])
  set(gca, 'TickLength', [0 0], 'Fontsize', 10.5, 'Fontname', 'Times', 'xlabel',...
  'Stress (f_s)', 'ylabel', ...
  'Strain (\epsilon_s)', 'linewidth', 1.5, 'box', 'off');
  xlabel('Strain (\epsilon_s) ^o/_{oo}')
  ylabel('Stress (f_s) MPa')
  
 %------------------------------------------------------------------------------- 
  subplot(steel,2, 2*cSS)
  p1 = plot(eps(cSS,:)/eps_y(cSS), sigmaDM(:,:,cSS)/fy(cSS), 'linewidth', 1.8, 'color', [1 0 0])
  %text(eps(cSS,end)*ones(length(L),1)/eps_y(cSS), sigmaDM(end,:,cSS)/fy(cSS), 'DM', 'Fontsize', 10.3, 'Fontname', 'Times')
  hold
  p2 = plot(eps(cSS,:)/eps_y(cSS), f_sc(:,:,cSS)/fy(cSS), 'linewidth', 1.8, 'color', [0 1 0])
  
  %text(eps(cSS,end)*ones(length(L),1)/eps_y(cSS), f_sc(end,:,cSS)/fy(cSS), 'RDM', 'Fontsize', 10.3, 'Fontname', 'Times')
  plot(eps(cSS,:)/eps_y(cSS), sigmaT(:,1,cSS)/fy(cSS), 'linewidth', 1.8, 'color', [0 0 0 ])
  text(-0.12*eps_u(cSS)/eps_y(cSS), failStress(cSS)/fy(cSS), num2str(failStress(cSS)/fy(cSS)), 'Fontsize', 8.5, 'Fontname', 'Times')
  plot([0 max(epsilon85fy(:,cSS))]/eps_y(cSS), failStress(cSS)*ones(1,2)/fy(cSS), 'linestyle', '--','linewidth', 0.7, 'color', [0 0 0 ])
  plot(epsilon85fy(:,cSS).*ones(3,2)/eps_y(cSS),[0 ;failStress(cSS)]/fy(cSS),'linestyle', '--','linewidth', 0.7, 'color', [0 0 0 ])
  text(epsilon85fy(:,cSS)'/eps_y(cSS)*1.05,0.2*failStress(cSS)*[1, .6, .5]/fy(cSS), num2str(epsilon85fy(:,cSS)/eps_y(cSS), '%2.1f'), 'color', [0 0 0 ], 'Fontsize', 8.3, 'Fontname', 'Times' )
  
  text(eps(cSS,50)*ones(length(L),1)/eps_y(cSS), (sigmaDM(50,:,cSS) + f_sc(50,:,cSS))/2/fy(cSS),...
  [num2str(L'./D')], 'Fontsize', 10.3, 'Fontname', 'Times')
  grid
  title(['f_y=', num2str(fy(cSS)/fy(cSS)), ' ', 'f_u=', num2str(fu(cSS)/fy(cSS)), ' ',...
  '\epsilon_y=', num2str(eps_y(cSS)/eps_y(cSS)), ' ','\epsilon_u=', num2str(eps_u(cSS)/eps_y(cSS))])
  set(gca, 'TickLength', [0 0], 'Fontsize', 10.5, 'Fontname', 'Times', 'xlabel',...
  'Stress (f_s)', 'ylabel', 'Strain (\epsilon_s)', 'linewidth', 1.5, 'box', 'off');
  xlabel('Normalized Strain (\epsilon_s/\epsilon_y)')
  ylabel('Normalized Stress (f_s/f_y)')
  if cSS == 1
    p = [p1(1), p2(2)]
    l = legend(p,'DM','RDM')
    set(l, 'Fontname', 'Times', 'box', 'off')
  end

end  


pkg load io
for i = 1:steel
  fN = [num2str(fy(i)), '.xlsx']
  xlswrite(fN, [eps(i,:)', sigmaT(:,1,i)])
endfor  


print('-fillpage','results','-dpdf')  





