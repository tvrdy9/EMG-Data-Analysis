function EMdelay = delaySlider(fsamp, iEMG, force, filtCST, filtForce, steadyEMGInd, steadyForceInd, titleString)
    time = 0:1/fsamp:(length(force)-1)/fsamp;
    padSecs = 2;
    if steadyEMGInd(1)-padSecs*fsamp < 1
        leftPad = length(filtCST(1:steadyEMGInd(1)-1))/fsamp;
    else
        leftPad = padSecs;
    end
    if steadyEMGInd(2)+padSecs*fsamp-1 > length(filtCST)
        rightPad = length(filtCST(steadyEMGInd(2):end))/fsamp;
    else
        rightPad = padSecs;
    end
    steadyForce = filtForce(steadyForceInd(1):steadyForceInd(2)-1);
    fluctForce = filtForce(steadyForceInd(1)-leftPad*fsamp:steadyForceInd(2)+rightPad*fsamp-1);
    fluctCST = filtCST(steadyEMGInd(1)-leftPad*fsamp:steadyEMGInd(2)+rightPad*fsamp-1);
    timePad = -leftPad:1/fsamp:10+rightPad-1/fsamp;
    time10 = 0:1/fsamp:10-1/fsamp;
    steadyForceOffset = double(steadyForceInd(1))/fsamp;
    EMdelay = 0;
    
    fig = uifigure("Name",titleString);
    g = uigridlayout(fig);
    g.RowHeight = {'1x','1x','fit'};
    g.ColumnWidth = {'1x'};
    ax1 = uiaxes(g);
    ax2 = uiaxes(g);

    plot(ax1,time,force/max(abs(force)),...
        timePad+steadyForceOffset,force(steadyForceInd(1)-leftPad*fsamp:steadyForceInd(2)+rightPad*fsamp-1)/max(abs(force)),...
        time10+steadyForceOffset,force(steadyForceInd(1):steadyForceInd(2)-1)/max(abs(force)),...
        time,iEMG/max(abs(iEMG)),...
        'LineWidth',2);
    legend(ax1,'Force','Padded Force','Steady 10s Force','iEMG','Location','southwest');
    plot(ax2,timePad,fluctForce/max(abs(fluctForce)),...
        time10,steadyForce/max(abs(fluctForce)),...
        timePad,fluctCST/max(abs(fluctCST)),...
        'LineWidth',2);
    legend(ax2,'Padded Force','Steady 10s Force','Padded CST','Location','northwest');
    ax2.Children(3).Color = "#D95319";
    ax2.Children(2).Color = "#EDB120";
    ax2.Children(1).Color = "#7E2F8E";
    
    sld = uislider(g,...
        "Limits",[-2000 2000],...
        "MajorTicks",[-2000:100:2000],...
        "Value",0);
    
    sld.ValueChangingFcn = @(src,event) updateEMdelay(src,event);
    
    uiwait(fig);
    
    function updateEMdelay(~,event)
        EMdelay = event.Value;
        updatePlots(EMdelay);
    end

    function updatePlots(delay)
        ax1.Children(2).XData = time10 + steadyForceOffset + delay/fsamp;
        ax1.Children(3).XData = timePad + steadyForceOffset + delay/fsamp;
        ax1.Children(4).XData = time + delay/fsamp;
        ax2.Children(2).XData = time10 + delay/fsamp;
        ax2.Children(3).XData = timePad + delay/fsamp;
    end
end