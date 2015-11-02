function [CE,FF,J_SC,V_OC,J_maxpower,V_maxpower] = SCperformance0(V,J,isShow)
% V - applied voltage
% J - current density

J_SC = J(1);
% J_SC = fzero(@(j)interp1(J(2:end),V(2:end),j,'PCHIP'),J(1)-1e-6);

ind = find(sign(J(1))~=sign(J)); 
if ~isempty(ind)
    V_OC = fzero(@(v)interp1(V(2:end),J(2:end),v,'PCHIP'),V(ind(1)));
else
    V_OC = fzero(@(v)interp1(V(2:end),J(2:end),v,'PCHIP'),0.6);
end    
V_maxpower = fzero(@(x)interp1((V(1:end-1)+V(2:end))/2,diff(V.*J)./diff(V),x,'PCHIP'),.6);
J_maxpower = interp1(V,J,V_maxpower,'PCHIP');
CE = V_maxpower*J_maxpower;
FF = V_maxpower*J_maxpower/V_OC/J_SC*100;

fprintf('\nJSC\t%7.4g mA/cm2\nVOC\t%7.4g mV\nCE\t%7.4g %%\nFF\t%7.4g %%',J_SC,V_OC/1e-3,CE,FF)

fprintf('\n\nfull precision:\nJSC\t%17.15g mA/cm2\nVOC\t%17.15g mV\nCE\t%17.15g %%\nFF\t%17.15g %%\n',J_SC,V_OC/1e-3,CE,FF)


if isShow
    figure, plot(V,J,'.-')
    hold on
    plot(V,V.*J,'.-','color','red')
    xlabel('Voltage, V')
    ylabel('Density of current, mA/cm{^2}')
    xlim_tmp = xlim;
    ylim_tmp = ylim;
    ylim_tmp(1) = -10;
    plot((V(1:end-1)+V(2:end))/2,diff(V.*J)./diff(V),'.-','color','cyan')
    xlim(xlim_tmp)
    ylim(ylim_tmp)

    plot(xlim,zeros(2,1),'color','green')
    plot(xlim,ones(2,1)*J_SC,'color','green')
    plot(ones(2,1)*V_OC,ylim,'color','green')

    plot(xlim,ones(2,1)*J_maxpower,'color','magenta')
    plot(ones(2,1)*V_maxpower,ylim,'color','magenta')
    string = sprintf('\nJ_{SC} = %g {mA/cm^2}\nV_{OC} = %g mV\n{\\eta} = %g%%\nFF = %g%%',J_SC,V_OC/1e-3,CE,FF);
    annotation('textbox', [.2 .6, .1, .1],'String' ,string)
end
